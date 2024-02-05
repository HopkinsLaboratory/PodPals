import os, time, datetime
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

def plot_TDDFT(root_directory, lower_bound, upper_bound, step_size, normalization):
    
    stime = time.time()
    print(f'{datetime.now().strftime('[ %H:%M:%S ]')} Starting extraction of TD-DFT spectra')

    # Get a list of all subdirectories in the root directory
    subdirectories = [f for f in os.listdir(root_directory) if os.path.isdir(os.path.join(root_directory, f))]

    # Iterate through each subdirectory
    for subdirectory in subdirectories:
        # Get a list of all .spectrum files in the current subdirectory
        subdirectory_path = os.path.join(root_directory, subdirectory)
        spectrum_files = [f for f in os.listdir(subdirectory_path) if f.endswith('.spectrum')]

        # Proceed if .spectrum files are found in the subdirectory
        if spectrum_files:
            # Find the maximum value from the first file
            first_file_path = os.path.join(subdirectory_path, spectrum_files[0])
            first_file_data = pd.read_csv(first_file_path, delim_whitespace=True)
            upper_bound = first_file_data['Energy'].max()

            # Create an empty DataFrame to store the interpolated spectra
            common_x_grid = np.arange(lower_bound, upper_bound + step_size, step_size)  # Lower bound hard-coded as 150 nm. Unlikely that users will need to go lower than this
            compiled_data = pd.DataFrame({'Energy': common_x_grid})

            # Process each spectrum file in the subdirectory
            for spectrum_file in spectrum_files:
                # Read data from the current spectrum file
                file_path = os.path.join(subdirectory_path, spectrum_file)
                data = pd.read_csv(file_path, delim_whitespace=True)

                # Interpolate the TotalSpectrum column on the common_x_grid
                interp_func = interp1d(data['Energy'], data['TotalSpectrum'], kind='linear', fill_value=0, bounds_error=False)
                interpolated_spectrum = interp_func(common_x_grid)

                # Add the interpolated spectrum to the DataFrame
                compiled_data[spectrum_file] = interpolated_spectrum

            # Sum all the spectra and add a column for the summed spectrum
            compiled_data['Summed Spectrum'] = compiled_data.iloc[:, 1:].sum(axis=1)

            #define the output file name + location of the respective subdirectory
            output_excel_file = os.path.join(subdirectory_path, f'{subdirectory}.xlsx')
            
            # Normalize the data if requested
            if normalization:
                compiled_data_norm = compiled_data.copy()
                compiled_data_norm.iloc[:, 1:] = compiled_data_norm.iloc[:, 1:] / compiled_data_norm.iloc[:, 1:].max()

                # Add the normalized data as a column
                compiled_data_norm['Normalized Summed Spectrum'] = compiled_data_norm['Summed Spectrum']

                # Write all data + the normalization to an Excel file
                compiled_data_norm.to_excel(output_excel_file, index=False)

            else:
                # Write all data to an Excel file 
                compiled_data.to_excel(output_excel_file, index=False)

    ftime = time.time()

    if normalization:
        print(f'{datetime.now().strftime('[ %H:%M:%S ]')} TD-DFT spectra with normalization have been plotted to excel files in each of the {len(subdirectories)} directories. Process completed in {np.round(ftime - stime, 2)} seconds.')
    else:
        print(f'{datetime.now().strftime('[ %H:%M:%S ]')} TD-DFT spectra have been plotted to excel files in each of the {len(subdirectories)} directories. Process completed in {np.round(ftime - stime, 2)} seconds.')

if __name__ == '__main__':
    # Set the root directory path
    root_directory = r'main_directory'
    lower_bound = 150 #lower bound for the TD-DFT extraction (in nm)
    upper_bound = 400 #upper bound for the TD-DFT extraction (in nm)
    step_size = 1 #step size for the interpolation (in nm)
    normalization = True

    # Call the function 
    plot_TDDFT(root_directory, lower_bound, upper_bound, step_size, normalization)
