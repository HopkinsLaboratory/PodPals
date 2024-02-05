import pandas as pd
import numpy as np
import os
from datetime import datetime

from scipy.constants import R
R = R * 1.E-3  # Convert gas constant to kJ/(mol*K)

def BW_CCS_Analysis(dft_file_path, dlpno_file_path, ccs_file_path, output_file_path, temp):

    def read_csv_and_clean_headers(file_path):
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Clean up the column headers by stripping spaces
        df.columns = df.columns.str.strip()
        
        # Clean up the column entries
        for col in df.columns:
            if df[col].dtype == object:
                df[col] = df[col].str.strip()

        return df

    def check_required_columns(df, required_columns, df_name):
        # Check if all required columns are present in the DataFrame
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"{df_name}: Missing required columns: {', '.join(missing_columns)}")

    def check_filenames(dft_data, dlpno_data, ccs_data):
        # Get the sets of filenames from each DataFrame
        dft_filenames = set(dft_data.index)
        ccs_filenames = set(ccs_data.index)
        filenames_to_check = [dft_filenames, ccs_filenames]
        
        # If DLPNO-CCSD(T) data is provided, add its filenames to the check
        if dlpno_data is not None:
            dlpno_filenames = set(dlpno_data.index)
            filenames_to_check.append(dlpno_filenames)
        
        # Check if all filenames match across datasets
        if not all(filenames == filenames_to_check[0] for filenames in filenames_to_check):
            print(f'DFT filenames: {dft_filenames}')
            print(f'DLPNOfilenames: {dlpno_filenames}')
            print(f'CCS filenames: {ccs_filenames}')
            raise ValueError('Filenames do not match across provided data files. Please ensure all files contain data for the same files.')

    def calculate_BW_CCS(dft_file_path, dlpno_file_path, ccs_file_path, temperature=298.15):
        # Read the DFT CSV file and preprocess filenames, and check that all data needed is present
        dft_data = read_csv_and_clean_headers(dft_file_path)
        check_required_columns(dft_data, ['Filename', 'Electronic energy', 'Gibbs Correction'], 'DFT file')

        dft_data['Filename'] = dft_data['Filename'].str.replace('.out', '', regex=False)
        dft_data.set_index('Filename', inplace=True)
        dft_data.sort_index(inplace=True)  # Order by filename
        dft_data.reset_index(drop=False, inplace=True)  # Reset index and keep 'Filename' as a column

        # Read the CCS CSV file and preprocess filenames
        ccs_data = read_csv_and_clean_headers(ccs_file_path)
        check_required_columns(ccs_data, ['filename', 'CCS [A**2]', 'errCCS [A**2]n'], 'CCS file')

        ccs_data['filename'] = ccs_data['filename'].str.replace('.mout', '', regex=False)
        ccs_data.set_index('filename', inplace=True)
        ccs_data.sort_index(inplace=True)  # Order by filename
        ccs_data.reset_index(drop=False, inplace=True)  # Reset index and keep 'filename' as a column

        # Read the DLPNO-CCSD(T) CSV file and preprocess filenames, if provided
        if not dlpno_file_path is None:
            dlpno_data = read_csv_and_clean_headers(dlpno_file_path)
            check_required_columns(dft_data, ['Filename', 'DLPNO-CCSD(T) energy'], 'CCSDT file')

            dlpno_data['Filename'] = dlpno_data['Filename'].str.replace('.out', '', regex=False)
            dlpno_data.set_index('Filename', inplace=True)
            dlpno_data.sort_index(inplace=True)  # Order by filename
            dlpno_data.reset_index(drop=False, inplace=True)  # Reset index and keep 'Filename' as a column

        # Check if filenames match across all provided data files, and only proceed if they do match
        check_filenames(dft_data, dlpno_data, ccs_data)

        # Initialize a new DataFrame for output data, and initialize with filenames
        output_data = pd.DataFrame(dft_data['Filename'])

        # Gibbs energies
        if not dlpno_file_path is None:
            # Calculate Gibbs energy with DLPNO-CCSD(T)
            output_data['Electronic energy'] = dlpno_data['DLPNO-CCSD(T) energy']  # Update electronic energy to DLPNO-CCSD(T) energy
            output_data['Gibbs energy'] = dlpno_data['DLPNO-CCSD(T) energy'] + dft_data['Gibbs Correction']
        else:
            # Calculate Gibbs energy without DLPNO-CCSD(T)
            output_data['Electronic energy'] = dft_data['Electronic energy'] 
            output_data['Gibbs energy'] = dft_data['Electronic energy'] + dft_data['Gibbs Correction']
        
        # Calculate Relative Gibbs Energy in kJ/mol
        output_data['Relative Gibbs Energy (kJ/mol)'] = (output_data['Gibbs energy'] - output_data['Gibbs energy'].min()) * 2625.5
        
        # Calculate populations, ensuring that the sum of relative populations equations to unity
        output_data['Population'] = np.exp(-output_data['Relative Gibbs Energy (kJ/mol)'] / (R * temperature))
        total_population = output_data['Population'].sum()
        output_data['Relative Population'] = output_data['Population'] / total_population
        
        # Ensure the sum of all relative populations equals 1
        if not np.isclose(output_data['Relative Population'].sum(), 1.0):
            raise ValueError('Sum of all relative populations does not equal 1; please check the data you are providing to the code.')
        
        # Add the CCS data to the dataframe, then calculate the BW CCS
        output_data['CCS (A**2)'] = ccs_data['CCS [A**2]']
        output_data['errCCS (A**2)'] = ccs_data['errCCS [A**2]']

        # Boltzmann weighting of the CCS
        BW_CCS = np.round((output_data['Relative Population'] * output_data['CCS (A**2)']).sum(), 2)
    
        return output_data, BW_CCS

    def calculate_BW_CCS_stdev(output_data):
        
        # The variance of a weighted sum can be calculated as the sum of the squares of the products of the weights and their corresponding variances. 
        # This works b/c the population values act as weights in the calculation of BW_CCS.
        variance_BW_CCS = np.sum(np.square(output_data['Relative Population'] * output_data['errCCS (A**2)']))
        stdev_bw_ccs = np.round(np.sqrt(variance_BW_CCS), 2)

        return stdev_bw_ccs  

    def write_to_excel(output_data, output_file_path, BW_CCS, BW_CCS_stdev):
        # Create a new column for BW_CCS that only contains the value once
        output_data['Boltzmann Weighted CCS (A**2)'] = ''
        output_data.loc[0, 'Boltzmann Weighted CCS (A**2)'] = BW_CCS 

        output_data['Boltzmann Weighted CCS stdev (A**2)'] = ''
        output_data.loc[0, 'Boltzmann Weighted CCS stdev (A**2)'] = BW_CCS_stdev 

        # Write to excel
        with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
            output_data.to_excel(writer, sheet_name='Processed Data', index=False)

    # Execute the function
    output_data, BW_CCS = calculate_BW_CCS(dft_file_path, dlpno_file_path, ccs_file_path, temp)
    BW_CCS_stdev = calculate_BW_CCS_stdev(output_data)
    write_to_excel(output_data, output_file_path, BW_CCS, BW_CCS_stdev)
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} BW CCS calculated successfully: {BW_CCS} +/- {BW_CCS_stdev} A**2')
    return

if __name__ == '__main__':
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T11_BW_CCS_Analyzer\File_sources'
    dft_file_path = os.path.join(directory, 'DFT', 'Thermo_data_298K_100kPa_1-0vibscl.csv')
    dlpno_file_path = os.path.join(directory, 'CCSDT', 'DLPNO_CCSDT_energies.csv') # set to none if not provided
    ccs_file_path = os.path.join(directory, 'CCS', 'Export_mout_CCS_lowfield.csv')
    output_file_path = os.path.join(directory, 'BW_CCS.xlsx')

    temp = 298.15

    BW_CCS_Analysis(dft_file_path, dlpno_file_path, ccs_file_path, output_file_path, temp)

