from Python.constants_and_conversions import convert_energy
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import os
import re
from openpyxl import Workbook
from datetime import datetime
import time

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('QtAgg') #Use matplotlib backend that is compatible w/ PyQt6 to prevent competition w/ GUI event loop

matplotlib.rcParams.update({'font.size': 18})

def extract_ESD_spectrum_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting):
   
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
        - y_final: Interpolated 'TotalSpectrum' values on the x_final grid after applying the shift.
        '''

        # Load the data
        #/t spacing is inconsistent in the header for some dumb reason, so we'll use regex to parse by whitespace
        data = pd.read_csv(os.path.join(directory, filename), sep=r'\s+')

        # Check that the dataframe has the expected headers
        required_columns = {'Energy', 'TotalSpectrum', 'IntensityFC', 'IntensityHT'}
        if not required_columns.issubset(data.columns):
            raise ValueError(f'{filename} does not contain the required columns: {required_columns}')
        
        # Convert energy from input unit to shift unit, and apply the shift
        try:
            data['Energy'] = data['Energy'].apply(lambda x: convert_energy(x, input_unit, shift_unit) + shift)
        except ValueError as e:
            raise ValueError(f'Error in conversion of energy from {filename} from the input unit to the shift unit: {e}')
        
        # Convert energy from shift unit to output unit
        try:
            data['Energy'] = data['Energy'].apply(lambda x: convert_energy(x, shift_unit, output_unit))
        except ValueError as e:
            raise ValueError(f'Error in conversion of energy from {filename} from the shift unit to the output unit: {e}')
        
        # Interpolate onto the final output grid
        f_interp = interp1d(data['Energy'], data['TotalSpectrum'], bounds_error = False, fill_value = 0)
        y_final = f_interp(x_final)

        return y_final

    def process_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting):
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
        summary_data[f'energy / {output_unit}'] = x_final

        # Prepare for plotting and add dynamic resizing to the plotting canvas
        if plotting:
            #need a way to know how many plots to make so that they all will fit on the output figure. This can be from the associated .spectrum files
            total_plots = len([name for name in os.listdir(directory) if name.endswith('.spectrum')])

            plots_per_row = 2
            rows = np.ceil(total_plots / plots_per_row).astype(int)
            cols = plots_per_row
            
            # Adjust figure size dynamically 
            fig_width = (cols * 8) 
            fig_height = (rows * 4) 
            plt.figure(figsize=(fig_width, fig_height))
            
            # Initialize plot counter
            plot_count = 1
                       
            # Add whitespace between panels
            plt.subplots_adjust(wspace=0.30, hspace=0.60, bottom=0.2) 
            
        for filename in os.listdir(directory):
            #get basename of each file from the corresponding .spectrum file
            if filename.endswith('.spectrum'):
                base_name = filename.split('.spectrum')[0]
                sheet_data = {}

                #need regex matching to get a basename. Files named as basename_1 will also accepted basename_12 if using startswith() 
                basename_pattern = re.escape(base_name) + r'(?=\.\s*spectrum\.root)'

                for root_file in os.listdir(directory):
                    if re.match(basename_pattern, root_file):
                        root_num = root_file.split('.spectrum.root')[1]

                        #shift/convert the spectral data by the specified amount/units, then interpolate on the final output range and write to excel file
                        y_final = process_VGFC_spectrum_file(root_file, directory, input_unit, shift_unit, output_unit, shift, x_final)
                        sheet_data[f'root{root_num}'] = y_final

                        #Add each root spectrum to the plot as transparent dots
                        if plotting:
                            ax = plt.subplot(rows, cols, plot_count)
                            plt.plot(x_final, y_final, label=f'{root_num}', linestyle='dotted', alpha=0.75, linewidth=1.5)

                        
                # Include x_final as a column and convert to DataFrame
                final_df = pd.DataFrame(sheet_data, index=x_final)
                final_df.reset_index(inplace=True)
                final_df.rename(columns={'index': f'energy / {output_unit}'}, inplace=True)

                # Calculate the sum of intensities for each energy level by dropping the energy column from the summation
                summed_intensity = final_df.drop(f'energy / {output_unit}', axis=1).sum(axis=1)
                final_df[f'Summed Intensity {base_name}'] = summed_intensity

                #Check the max of the summed intensity - if it is below 1000, then the absorption is very weak and the user has probably specificed an incorrect range to output the data and/or has applied an incorrect shift to their data.
                max_intensity = final_df[f'Summed Intensity {base_name}'].max()

                if max_intensity < 1250:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {filename} is returning a very small absorption between {output_lowervalue} - {output_lowervalue} {output_unit}. Check that you have applied a correct shift, and that you are outputting spectra in a region where your analyte actually absorbs.')

                # Normalize the summed intensities if requested
                if normalize_output:
                    final_df[f'Normalized Intensity {base_name}'] = final_df[f'Summed Intensity {base_name}'] / max_intensity
                    summary_data[f'{base_name}'] = final_df[f'Normalized Intensity {base_name}'].values

                else:
                    summary_data[f'{base_name}'] = final_df[f'Summed Intensity {base_name}'].values

                # Save the root data to the roots Excel file
                final_df.to_excel(writer, sheet_name=base_name, index=False)

                #Add the summed intensity to the plot
                if plotting:
                    plt.plot(final_df[f'energy / {output_unit}'], final_df[f'Summed Intensity {base_name}'], label='Summed Intensity', color='k', linestyle='solid', linewidth = 1.5)
                    plt.title(base_name)
                    plt.xlabel(f'Energy / {output_unit}')
                    plt.ylabel('Intensity')

                    plot_count += 1


        writer.close()

        #Prepare and write the summary file
        final_df_summary = pd.DataFrame(summary_data)
        final_df_summary.to_excel(total_xlsx, index=False)   

        #Prepare and create the plot of the spectra (if requested)
        if plotting:
           
            #Generate a filename for the plot to be saved as that won't interfere w/ previous versions
            save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_wroots.png')
            
            i = 2
            while os.path.isfile(save_file):
                save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_wroots_v{i}.png')
                i += 1
            
            # Save the plot
            plt.tight_layout() 
            plt.savefig(save_file)
            #plt.show()

            #Now for the total spectra plot
            plt.figure(figsize=(20, 10))
            
            #Generate a filename for the plot to be saved as that won't interfere w/ previous versions
            save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum.png')
            i = 2 
            while os.path.isfile(save_file):
                save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_v{i}.png')
                i += 1  
            
            # Plot each total spectrum, skipping first entry in the summary_data dict because that is the energy column 
            for key, y_final in list(summary_data.items())[1:]:
                plt.plot(x_final, y_final, linestyle='solid', alpha=0.75, label=key)
            
            # Add plot details
            plt.title('Total Absorbtion Spectrum')
            plt.xlabel(f'Energy {output_unit}')  # Adjust the label as per your unit
            if normalize_output:
                plt.ylabel('Norm. Intensity')

            else:
                plt.ylabel('Intensity')
            
            plt.legend(loc='best')  # Show legend to identify each spectrum
            plt.grid(True)
            
            #Save the plot
            plt.savefig(save_file)
            #plt.show()

    #Execute the function  
    stime = time.time()
    process_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting) 

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} .spectrum.root files have been processed in {np.round((time.time() - stime), 2)} seconds.')
    return

if __name__ == '__main__':
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\ESD'
    input_unit = 'nm' # Can be 'wavelength', 'wavenumber', or 'eV'
    output_unit = 'nm' # Can be 'wavelength', 'wavenumber', or 'eV'
    shift = -0.45 # Must be a number or float 
    shift_unit = 'eV' # Can be 'wavelength', 'wavenumber', or 'eV'
    output_lowervalue = 200 # Must be a number or float
    output_uppervalue = 400 # Must be a number or float
    output_spacing = 1 # Must be a number or float
    output_file_basename = 'Fentanyl_OProt' #basename only! no file extension!

    normalize_output = False #normalize combined spectra to unity
    plotting = True

    '''Main section'''
    #check for .spectrum.root and .spectrum files, and proceed with processing if they are found
    contains_spectrum_root = any('.spectrum.root' in file for file in os.listdir(directory))
    ends_with_spectrum = any(file.endswith('.spectrum') for file in os.listdir(directory))

    #If directory contains both .spectrum and .spectrum.root files or just .spectrum.root files, the full analysis can be done using only .spectrum.root files:
    if (contains_spectrum_root and ends_with_spectrum):
        
        stime = time.time()
        
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Both .spectrum.root and .spectrum files were found in the provided directory. All spectra will be processed from the root files. Please ensure that all .spectrum.root files are present, otherwise the analysis will be incomplete.')
        extract_ESD_spectrum_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting)
        
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} .spectrum.root files have been processed in {np.round((time.time() - stime), 2)} seconds.')
    
    #If directory does not both .spectrum files and .spectrum.root, inform the user:
    else: 
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .spectrum or .spectrum.root files were found in {directory}. Please ensure that you have specified the correct directory.')
