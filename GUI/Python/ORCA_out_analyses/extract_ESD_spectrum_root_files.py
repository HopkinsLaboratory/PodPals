from Python.constants_and_conversions import convert_energy
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import os
import re
from openpyxl import Workbook
from datetime import datetime
import time

from PyQt6.QtWidgets import QApplication

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('QtAgg') #Use matplotlib backend that is compatible w/ PyQt6 to prevent competition w/ GUI event loop

matplotlib.rcParams.update({'font.size': 20})

def extract_ESD_spectrum_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting, plots_per_row):
   
    def process_VGFC_spectrum_file(filename, directory, input_unit, shift_unit, output_unit, shift, x_final):
        '''
        Processes spectrum data to interpolate 'TotalSpectrum' based on energy conversion and shift.

        Parameters:
        - filename: Name of the file containing the data.
        - directory: Directory where the file is located.
        - input_unit, shift_unit, output_unit: Units for the energy conversion process.
        - shift: Value to be added to the energy during the conversion process.
        - x_final: Output interpolation range defined by output_lowervalue, output_uppervalue, output_spacing

        Returns:
        - y_final: Interpolated 'TotalSpectrum' values on the x_final grid after applying the shift or None if an error is encountered during the processing.
        '''

        #Load the data
        #/t spacing is inconsistent in the header for some dumb reason, so we'll use regex to parse by whitespace
        data = pd.read_csv(os.path.join(directory, filename), sep=r'\s+')

        #Check that the dataframe has the expected headers
        required_columns = {'Energy', 'TotalSpectrum', 'IntensityFC', 'IntensityHT'}
        if not required_columns.issubset(data.columns):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The header of {filename} does not contain the required columns: {", ".join(required_columns)}\n. Processing of this file will be skipped.')
            return None
        
        #Convert energy from input unit to shift unit, and apply the shift
        try:
            data['Energy'] = data['Energy'].apply(lambda x: convert_energy(x, input_unit, shift_unit) + shift)
        
        except ValueError as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error in conversion of energy from {filename} from the input unit to the shift unit: {e}')
            return None
        
        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An unexpected error occured when converting the energy from {filename} from the input unit to the shift unit: {e}.\nPlease report this error to the issue section of the GitHub repo.')
            return None
                
        #Convert energy from shift unit to output unit
        try:
            data['Energy'] = data['Energy'].apply(lambda x: convert_energy(x, shift_unit, output_unit))
            
            #Interpolate onto the final output grid
            f_interp = interp1d(data['Energy'], data['TotalSpectrum'], bounds_error = False, fill_value = 0)
            y_final = f_interp(x_final)

            return y_final
        
        except ValueError as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error in conversion of energy from {filename} from the shift unit to the output unit: {e}')
            return None
        
        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An unexpected error occured when converting the energy from {filename} from the shift unit to the output unit: {e}.\nPlease report this error to the issue section of the GitHub repo.')
            return None

    def process_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting, plots_per_row):
        '''
        Processes each .spectrum.root file in the specified directory by converting energy units, applying a shift, and interpolating the total spectrum data. The function groups data according to the base name extracted from associated .spectrum files. It requires both .spectrum and .spectrum.root files to be present in the directory. The processed data can be normalized and plotted based on the provided parameters.

        Parameters:
        - directory (str): Path to the directory containing the .spectrum and .spectrum.root files.
        - input_unit (str): The unit of the input energy values ('nm' for wavelength, 'cm**-1' for wavenumber, or 'eV' for electron volts).
        - shift_unit (str): The unit to which the energy values will be temporarily converted for applying the shift ('nm', 'cm**-1', or 'eV').
        - output_unit (str): The unit of the output energy values after conversion and shifting ('nm', 'cm**-1', or 'eV').
        - shift (float): The value to be added to the energy during conversion. Can be positive or negative.
        - output_lowervalue (float): The lower bound of the output energy range for interpolation.
        - output_uppervalue (float): The upper bound of the output energy range for interpolation.
        - output_spacing (float): The spacing between consecutive energy values in the output interpolated data.
        - output_file_basename (str): The base name for the output Excel file. This should not include a file extension.
        - normalize_output (bool): If True, the combined spectra will be normalized to unity. If False, raw intensity values are used.
        - plotting (bool): If True, generates a plot of the processed spectra. Each root file's spectrum is plotted with slightly transparent dots, and the combined spectrum is also plotted if multiple root files are processed.

        The function saves the processed spectra to an Excel file named after `output_file_basename` with '_roots.xlsx' appended. The Excel file contains multiple sheets, each corresponding to a different base name group of .spectrum.root files.
        '''

        #ensure that no file extension is in the output file basename
        output_file_basename = os.path.splitext(output_file_basename)[0]

        #Write all roots to a unique Excel file name, ensuring that previous versions do not get overwritten
        roots_xlsx = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_wroots.xlsx')
        i = 2 
        while os.path.isfile(roots_xlsx):
            roots_xlsx = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_wroots_v{i}.xlsx')
            i += 1

        #Write all Summed Spectra to a unique Excel file name, ensuring that previous versions do not get overwritten
        total_xlsx = os.path.join(directory, f'{output_file_basename}_TotalSpectrum.xlsx')
        i = 2 
        while os.path.isfile(total_xlsx):
            total_xlsx = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_v{i}.xlsx')
            i += 1

        writer = pd.ExcelWriter(roots_xlsx, engine='openpyxl')

        #Define the final interpolation range
        x_final = np.arange(output_lowervalue, output_uppervalue + output_spacing, output_spacing)

        #Dictionary to write the summed and normalized (if requested) spectra, with first column being the energy range
        summary_data = {}

        #Prepare for plotting roots and add dynamic resizing to the plotting canvas
        if plotting:
            #need a way to know how many plots to make so that they all will fit on the output figure. This can be from the associated .spectrum files
            total_plots = len([name for name in os.listdir(directory) if name.endswith('.spectrum')])

            #plots_per_row = 2 #can be adjusted to user preference
            rows = np.ceil(total_plots / plots_per_row).astype(int)
            cols = plots_per_row
            
            #Adjust figure size dynamically 
            fig, axs = plt.subplots(rows, cols, figsize=((cols * 8), (rows * 4)))
            axs = axs.flatten() if total_plots > 1 else [axs]  #Ensure axs is always a list for consistency

            #Hide any unused subplots if total_plots isn't a perfect multiple of plots_per_row
            for i in range(total_plots, len(axs)):
                axs[i].axis('off')
            
            #Initialize plot counter
            plot_count = 0
                                  
        for filename in os.listdir(directory):
            #get basename of each file from the corresponding .spectrum file
            if filename.endswith('.spectrum'):
                base_name = filename.split('.spectrum')[0]
                sheet_data = {}

                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Processing .spectrum.rootN files from {base_name}.')
                QApplication.processEvents()

                #need regex matching to get a basename. Files named as basename_1 will also accepted basename_12 if using .startswith(basename_1) logic - found that out the hard way :(
                basename_pattern = re.escape(base_name) + r'(?=\.\s*spectrum\.root)'

                #look again through the specified directory, now searching for the .spectrum.rootN files
                for root_file in os.listdir(directory):
                    if re.match(basename_pattern, root_file):
                        root_num = root_file.split('.spectrum.root')[1]

                        #shift/convert the spectral data by the specified amount/units, then interpolate on the final output range and write to excel file - process_VGFC_spectrum_file returns None on an erronous extraction
                        y_final = process_VGFC_spectrum_file(root_file, directory, input_unit, shift_unit, output_unit, shift, x_final)
                        
                        if y_final is not None:
                            sheet_data[f'root{root_num}'] = y_final

                            #Add each root spectrum to a subplot plot as slightly transparent lines; each subplot corresponds to a unique basename defined by the prefix to .spectrum files
                            if plotting:
                                ax = axs[plot_count]
                                ax.plot(x_final, y_final, label=f'Root {root_num}', linestyle='solid', alpha=0.75, linewidth=1.5)
                        
                        #If spectral data could not be extracted from the .spectrum.rootN file, skip the processing and proceed to the next file
                        else:
                            continue      
                                
                #Include x_final as a column and convert to DataFrame
                final_df = pd.DataFrame(sheet_data, index=x_final)
                final_df.reset_index(inplace=True)
                final_df.rename(columns={'index': f'energy / {output_unit}'}, inplace=True)

                #Calculate the sum of intensities for each energy level by dropping the energy column from the summation
                summed_intensity = final_df.drop(f'energy / {output_unit}', axis=1).sum(axis=1)
                final_df[f'Summed Intensity {base_name}'] = summed_intensity

                #Check the max of the summed intensity - if it is below 1000, then the absorption is very weak and the user has probably specificed an incorrect range to output the data and/or has applied an incorrect shift to their data.
                max_intensity = final_df[f'Summed Intensity {base_name}'].max()

                if max_intensity < 800:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The root files that contribute to {filename} are returning a small absorption between {output_lowervalue} - {output_lowervalue} {output_unit}. Check that you have applied a correct shift, and that you are outputting spectra in a region where your analyte actually absorbs.')
                    pass

                if np.isclose(max_intensity, 0, atol=1E-3):
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The root files that contribute to {filename} are returning a zero absorption between {output_lowervalue} - {output_lowervalue} {output_unit}. Check that you have applied a correct shift, and that you are outputting spectra in a region where your analyte actually absorbs.\nProcessing of this file will be skipped.')
                    continue

                #Normalize the summed intensities if requested
                if normalize_output:
                    final_df[f'Normalized Intensity {base_name}'] = final_df[f'Summed Intensity {base_name}'] / max_intensity
                    summary_data[f'{base_name}'] = final_df[f'Normalized Intensity {base_name}'].values

                else:
                    summary_data[f'{base_name}'] = final_df[f'Summed Intensity {base_name}'].values

                #Save the root data to the roots Excel file
                final_df.to_excel(writer, sheet_name=base_name, index=False)

                #Add the summed intensity and labels to the subplot 
                if plotting:
                    ax = axs[plot_count]
                    ax.plot(final_df[f'energy / {output_unit}'], final_df[f'Summed Intensity {base_name}'], label='Summed Intensity', color='k', linestyle='solid', linewidth = 1.5)
                    ax.set_title(base_name)
                    ax.set_xlabel(f'Energy / {output_unit}')
                    ax.set_ylabel('Intensity')

                    plot_count += 1

        writer.close()

        #Prepare and write a summary file containg all the summed intesities / normalized summed intensities, with energies being added to the df as the first column
        final_df_summary = pd.DataFrame(summary_data, index=x_final)
        final_df_summary.reset_index(inplace=True)
        final_df_summary.rename(columns={'index': f'energy / {output_unit}'}, inplace=True)
        final_df_summary.to_excel(total_xlsx, index=False)   

        #Prepare and create the plot of the spectra (if requested)
        if plotting:
           
            '''Plotting the spectrum of each file showing contributin from each root as a separate plot.'''
            
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Plotting data from all .spectrum.rootN files....')
            QApplication.processEvents()
            
            #Because all subplots have the same data labels, we can set up the legend using handles and labels from the first subplot
            if len(axs) > 0:
                handles, labels = axs[0].get_legend_handles_labels()

            plt.tight_layout()
            plt.subplots_adjust(bottom=0.11) #spacing to accomodate the legend

            #Because of the dynamic sizing of the plots, to minimize headaches, there needs to be a consistent legend size. In this case, the legend will always be three rows.
            num_items = len(handles)
            num_cols = np.ceil(num_items / 3)  #This will ensure 3 rows
            
            #Set the legend dynamically below the subplots
            #legend_offset = -1.1306 * exp(-0.8627 * rows) - exp fit to data for rows 1 - 7
            fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, (-1.1306 * np.exp(-0.8627 * rows))), ncol=num_cols, fancybox=True, shadow=True)
       
            #Generate a filename for the plot to be saved as that won't interfere w/ previous versions
            save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_wroots.png')
            
            i = 2
            while os.path.isfile(save_file):
                save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_wroots_v{i}.png')
                i += 1  
           
            #Save the plot
            plt.savefig(save_file, bbox_inches='tight')

            '''Plotting the total spectrum '''
            plt.figure(figsize=(16, 10))
            
            #Generate a filename for the plot to be saved as that won't interfere w/ previous versions
            save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum.png')
            
            i = 2 
            while os.path.isfile(save_file):
                save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_v{i}.png')
                i += 1  

            #Initialize an empty list to store modified filenames as entries for the plot legend
            plt_labels = []            
            
            for key, y_final in summary_data.items():
                
                #Check if the label length exceeds a maximum length so that the legend doesn't get hella cluttered
                if len(key) > 20:
                    match = re.search(r"(\d+)(?!.*\d)", key)  #regex to look for the last number in a filename
                    basename = key.split('_')[0]

                    if match and basename:
                        new_label = f'{basename}_{match.group(1)}' #Use first text preceeded by an underscore and the last number as the legend entry
                    
                    elif match and not basename:
                        new_label = f'File_{match.group(1)}' #If no underscores in fileanme, use the 'File' prefix and the last number
                    
                    else:
                        new_label = key[:int(np.ceil(len(key) * 0.25))]  #Take first 25% of characters in the filename and remove the file extension [:-4]... users can adjust this number if they don't like it
                
                #If a filename (key) is short enough, we can useit as is
                else:
                    new_label = key

                plt_labels.append(new_label)  #Store the (potentially) modified label                
                plt.plot(x_final, y_final, linestyle='solid', alpha=0.75, label=new_label)
            
            #Add plot details
            plt.title('Total Absorbtion Spectrum')
            plt.xlabel(f'Energy / {output_unit}') 

            if normalize_output:
                plt.ylabel('Norm. Intensity')
            else:
                plt.ylabel('Intensity')

            plt.tight_layout()
            plt.subplots_adjust(bottom=0.5) #spacing to accomodate the legend
            
            #Use the labels list for the legend and adjust its placement
            plt.legend(labels=plt_labels, loc='lower center', bbox_to_anchor=(0.5, -0.55), ncol=5, fancybox=True, shadow=True)
            plt.grid(True)
            
            #Save the plot
            plt.savefig(save_file, bbox_inches='tight')

    #Execute the function  
    stime = time.time()
    process_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting, plots_per_row) 

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} .spectrum.root files have been processed in {np.round((time.time() - stime), 2)} seconds.')
    return

if __name__ == '__main__':
    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T10_ORCA_extract_VGFC_Spectra\spectrum_and_root_files'
    input_unit = 'nm' #Can be 'wavelength', 'wavenumber', or 'eV'
    output_unit = 'nm' #Can be 'wavelength', 'wavenumber', or 'eV'
    shift = -0.25 #Must be a number or float 
    shift_unit = 'eV' #Can be 'wavelength', 'wavenumber', or 'eV'
    output_lowervalue = 150 #Must be a number or float
    output_uppervalue = 400 #Must be a number or float
    output_spacing = 1 #Must be a number or float
    output_file_basename = 'Fentanyl_OProt' #basename only! no file extension!

    normalize_output = False #normalize combined spectra to unity
    plotting = True
    plots_per_row = 2

    '''Main section'''
    #check for .spectrum.root and .spectrum files, and proceed with processing if they are found
    contains_spectrum_root = any('.spectrum.root' in file for file in os.listdir(directory))
    ends_with_spectrum = any(file.endswith('.spectrum') for file in os.listdir(directory))

    #If directory contains both .spectrum and .spectrum.root files or just .spectrum.root files, the full analysis can be done using only .spectrum.root files:
    if (contains_spectrum_root and ends_with_spectrum):
        
        stime = time.time()
        
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Both .spectrum.root and .spectrum files were found in the provided directory. All spectra will be processed from the root files. Please ensure that all .spectrum.root files are present, otherwise the analysis will be incomplete.')
        extract_ESD_spectrum_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting, plots_per_row)
        
        #print(f'{datetime.now().strftime("[ %H:%M:%S ]")} .spectrum.root files have been processed in {np.round((time.time() - stime), 2)} seconds.')
    
    #If directory does not both .spectrum files and .spectrum.root, inform the user:
    else: 
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .spectrum or .spectrum.root files were found in {directory}. Please ensure that you have specified the correct directory.')
