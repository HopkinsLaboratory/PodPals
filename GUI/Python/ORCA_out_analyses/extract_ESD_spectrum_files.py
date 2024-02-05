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

def extract_ESD_spectrum_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting):
   
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

    #process each .spectrum file that corresponds to each filename, and write each spectrum to a single excel file
    def process_spectrum_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting):
        '''
        Processes each .spectrum file in the specified directory by converting energy units, applying a shift, and interpolating the total spectrum data. The function groups data according to the base name extracted from associated .spectrum files. It requires both .spectrum and .spectrum.root files to be present in the directory. The processed data can be normalized and plotted based on the provided parameters.

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

        The function saves the processed spectra to an Excel file named after `output_file_basename` with '_TotalSpectra.xlsx' appended.
        '''

        summary_data = {}

        #ensure that no file extension is in the output file basename
        output_file_basename = os.path.splitext(output_file_basename)[0]

        #define final interpolation range
        x_final = np.arange(output_lowervalue, output_uppervalue + output_spacing, output_spacing)

        #Write all roots to a unique Excel file name, ensuring that previous versions do not get overwritten
        total_xlsx = os.path.join(directory, f'{output_file_basename}_TotalSpectrum.xlsx')
        i = 2 
        while os.path.isfile(total_xlsx):
            total_xlsx = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_v{i}.xlsx')
            i += 1

        for filename in os.listdir(directory):
            if filename.endswith('.spectrum'):
                base_name = filename.split('.spectrum')[0]
                
                #shift/convert the spectral data by the specified amount/units, then interpolate on the final output range and write to excel file
                y_final = process_VGFC_spectrum_file(filename, directory, input_unit, shift_unit, output_unit, shift, x_final)
                max_intensity = np.max(y_final)

                if max_intensity < 1250:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {filename} is returning a very small absorption between {output_lowervalue} - {output_lowervalue} {output_unit}. Check that you have applied a correct shift, and that you are outputting spectra in a region where your analyte actually absorbs.')
                
                if normalize_output:
                    #normalize y_final to unity
                    y_final = y_final / max_intensity

                summary_data[f'{base_name}'] = y_final

                # Include x_final as a column and convert to DataFrame
                final_df = pd.DataFrame(summary_data, index=x_final)
                final_df.reset_index(inplace=True)
                final_df.rename(columns={'index': f'energy / {output_unit}'}, inplace=True)

                final_df.to_excel(total_xlsx, index=False)                    

        if plotting:
            plt.figure(figsize=(20, 10))
            
            #Generate a filename for the plot to be saved as that won't interfere w/ previous versions
            save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum.png')
            i = 2 
            while os.path.isfile(save_file):
                save_file = os.path.join(directory, f'{output_file_basename}_TotalSpectrum_v{i}.png')
                i += 1  
            
            # Plot each total spectrum 
            for key, y_final in summary_data.items():
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

    #Execute the main function
    stime = time.time()
    process_spectrum_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting)

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} .spectrum files have been processed in {np.round((time.time() - stime), 2)} seconds.')
    return

if __name__ == '__main__':
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\ESD_Test2'
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
   
    #If directory contains .spectrum files, but not .spectrum.root files:
    if ends_with_spectrum and not contains_spectrum_root:
        
        stime = time.time()

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Only .spectrum files were found in the provided directory. Only total spectra will be processed. To see the contributions from each excited state, please add the respective .spectrum.root files to the specified directory.')
        extract_ESD_spectrum_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting)

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} .spectrum files have been processed in {np.round((time.time() - stime), 2)} seconds.')

    #If directory does not contain any .spectrum files or .spectrum.root, inform the user:
    else: 
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .spectrum or .spectrum.root files were found in {directory}. Please ensure that you have specified the correct directory.')