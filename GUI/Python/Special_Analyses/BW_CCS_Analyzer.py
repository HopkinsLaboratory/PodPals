import pandas as pd
import numpy as np
import os
from datetime import datetime
from PyQt6.QtWidgets import QApplication

from scipy.constants import R
R = R * 1.E-3  #Convert gas constant to kJ/(mol*K)

def BW_CCS_Analysis(dft_file_path, dlpno_file_path, ccs_file_path, output_file_path, temp):

    def read_csv_and_clean_headers(file_path):
        #Read the CSV file
        try:
            df = pd.read_csv(file_path)
            
            #Clean up the column headers by stripping spaces
            df.columns = df.columns.str.strip()
            
            #Clean up the column entries
            for col in df.columns:
                if df[col].dtype == object:
                    df[col] = df[col].str.strip()

            return df

        except FileNotFoundError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The file {os.path.basename(file_path)} could not be found.')

        except PermissionError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The file {os.path.basename(file_path)} could not be accessed due to a permission error. Is it open or in use by another program?.')

        except pd.errors.EmptyDataError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {os.path.basename(file_path)} does not contain any data!')

        except pd.errors.ParserError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {os.path.basename(file_path)} could not be parsed! Is it comma delimited?')
        
        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An unexpected exception occured when reading {os.path.basename(file_path)}: {e}')
            print('Please report this error to the Issues section of the GitHub repo.')

        return None  #On error, return None

    def check_required_columns(df, required_columns, df_name):
        #Check if all required columns are present in the DataFrame
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {df_name} is missing the required columns: {', '.join(missing_columns)}')
            return False
        
        else:
            return True

    def check_filenames(dft_data, dlpno_data, ccs_data):
        
        #Extract the filenames directly from the 'filenames' column of each DataFrame, and remove the file extension from the filename
        dft_filenames = dft_data['Filename'].tolist()
        ccs_filenames = ccs_data['filename'].tolist()
        
        #If DLPNO-CCSD(T) data is provided, add its filenames to the check
        if dlpno_data is not None:
            dlpno_filenames = dlpno_data['Filename'].tolist()
        
        #Initialize a flag to indicate if all filenames match
        all_match = True
        
        #Filename checking procedure for when DLPNO-CCSD(T) data is provided
        if dlpno_data is not None:

            #Check lengths first to ensure they can be compared
            if len(dft_filenames) == len(ccs_filenames) == len(dlpno_filenames):
                
                #Compare each index position across all three lists
                for i in range(len(dft_filenames)):
                    if not (dft_filenames[i] == ccs_filenames[i] == dlpno_filenames[i]):
                        all_match = False
                        break            
            else:
                all_match = False
        
        #If no DLPNO-CCSD(T) data IS GIVEN, compare DFT and CCS only
        else:
            #Check lengths first to ensure they can be compared
            if len(dft_filenames) != len(ccs_filenames):
                all_match = False
            else:
                #Compare each index position across all three lists
                for i in range(len(dft_filenames)):
                    if dft_filenames[i] != ccs_filenames[i]:
                        all_match = False
                        break
            
        #Check if all filenames match across datasets. If they don't, let the user know and halt analysis
        if not all_match:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Filenames do not match across provided data files. Please ensure all files contain data for the same files.')
            print(f'DFT filenames: {", ".join(dft_filenames)}')
            
            if dlpno_data is not None:
                print(f'Coupled cluster filenames: {", ".join(dlpno_filenames)}')
            
            print(f'CCS filenames: {", ".join(ccs_filenames)}')
            
            return False
        
        #If filenames match, return True
        else:
            return True

    def calculate_BW_CCS(dft_file_path, dlpno_file_path, ccs_file_path, temperature=298.15):
        
        #Read the DFT CSV file and preprocess filenames, and check that all data needed is present
        dft_data = read_csv_and_clean_headers(dft_file_path) #returns None on error and the csv data stored in a dataframe on success

        if dft_data is not None:
            #check that all the required columns are present in the DFT thermochem .csv
            dft_check = check_required_columns(dft_data, ['Filename', 'Electronic energy', 'Gibbs Correction'], 'DFT file') #returns false
            
            if dft_check:
                #remove the file extension from the filenames in the Filename column
                dft_data['Filename'] = dft_data['Filename'].str.replace('.out', '', regex=False)
                dft_data.set_index('Filename', inplace=True)
                dft_data.sort_index(inplace=True)  #Order by filename
                dft_data.reset_index(drop=False, inplace=True)  #Reset index and keep 'Filename' as a column

        #Return None for the output data and BW_CCS on any error encountered
            else:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
                return None, None
            
        else:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
            return None, None

        #Read the CCS CSV file and preprocess filenames, and check that all data needed is present
        ccs_data = read_csv_and_clean_headers(ccs_file_path) #returns None on error and the csv data stored in a dataframe on success

        if ccs_data is not None:
            ccs_check = check_required_columns(ccs_data, ['filename', 'CCS [A**2]', 'errCCS [A**2]'], 'CCS file')

            if ccs_check:
                #remove the file extension from the filenames in the Filename column
                ccs_data['filename'] = ccs_data['filename'].str.replace('.mout', '', regex=False)
                ccs_data.set_index('filename', inplace=True)
                ccs_data.sort_index(inplace=True)  #Order by filename
                ccs_data.reset_index(drop=False, inplace=True)  #Reset index and keep 'filename' as a column

            else:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
                return None, None
            
        else:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
            return None, None
        
        #Read the DLPNO-CCSD(T) CSV file and preprocess filenames, if provided
        if not dlpno_file_path is None:
            
            dlpno_data = read_csv_and_clean_headers(dlpno_file_path)
            if dlpno_data is not None:
                dlpno_check = check_required_columns(dlpno_data, ['Filename', 'CCSD(T) energy', 'CCSD energy'], 'CCSDT file') #returns None on error and the csv data stored in a dataframe on success

                if dlpno_check:
                    #remove the file extension from the filenames in the Filename column
                    dlpno_data['Filename'] = dlpno_data['Filename'].str.replace('.out', '', regex=False)
                    dlpno_data.set_index('Filename', inplace=True)
                    dlpno_data.sort_index(inplace=True)  #Order by filename
                    dlpno_data.reset_index(drop=False, inplace=True)  #Reset index and keep 'Filename' as a column

                else:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
                    return None, None
                
            else:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
                return None, None
        
        else:
            dlpno_data = None
            
        #Once .csv files containing data are verified to have the correct information, ensure that the filename basenames are in the same order and have the same basename. Otherwise, the Boltzmann weighting will not work like we want it to!
        filename_check = check_filenames(dft_data, dlpno_data, ccs_data)

        if not filename_check:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
            return None, None

        #Woohoo! Making it to this point means that all hte data is valid, and we can proceed to calcualte the BW CCSs!
        
        #Initialize a new DataFrame for output data, and initialize with filenames (Choosing to initialize w/ DFT data arbitrarily - it could be from any of the .csv files)
        output_data = pd.DataFrame(dft_data['Filename'])

        #Gibbs energies - can be calculated using DFT electronic energies fo DLPNO-CCSD(T) electronic energies; the method depends on the files provided
        #For calculating Gibbs energy using coupled cluster electronic energies
        if not dlpno_file_path is None:
            
            #Calculate Gibbs energy with CCSD(T) if it exists; if not try CCSD energy
            try:
                ccsdt_energy = dlpno_data.get('CCSD(T) energy', 'N/A')  #Set data to 'N/A' or another default if the key doesn't exist, but it should because we checked for it previously - redundant error handling!
                
                #check if ccsdt_energy is numeric by by getting the mean - will throw a type error if there are any non numeric entries present (like N/A)
                ccsdt_energy_check = ccsdt_energy.mean()
                
                #let the user know that CCSD(T) energies are being used
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} CCSD(T) energies and Gibbs corrections were successfully read in from {os.path.basename(dlpno_file_path)} and {os.path.basename(dft_file_path)}, respectively.')
                QApplication.processEvents()

                output_data['Electronic energy'] = dlpno_data['CCSD(T) energy']  #Update electronic energy to CCSD(T) energy
                output_data['Gibbs energy'] = dlpno_data['CCSD(T) energy'] + dft_data['Gibbs Correction']

            except TypeError: 

                try:
                    #Read in ccsd energies and check if they're numeric by getting the mean - will throw a type error if there are any non numeric entries present (like N/A)
                    ccsd_energy = dlpno_data.get('CCSD energy', 'N/A')
                    ccsd_energy_check = ccsd_energy.mean()

                    output_data['Electronic energy'] = dlpno_data['CCSD energy']  #Update electronic energy to CCSD energy if CCSD(T) was not calculated
                    output_data['Gibbs energy'] = dlpno_data['CCSD energy'] + dft_data['Gibbs Correction']

                    #let the user know that CCSD energies are being used instead
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} CCSD(T) energies could not be found in {os.path.basename(dlpno_file_path)}. CCSD energies are being used instead with Gibbs corrections from and {os.path.basename(dft_file_path)}, both of which have been read in correctly.')
                    QApplication.processEvents() 

                except TypeError:
                    #let the user know that both CCSD(T) and CCSD entries are non-numeric and cannot be used - the file needs to be checked for correctness.
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} CCSD(T) or CCSD energies in {os.path.basename(dlpno_file_path)} were both non-numeric. Did the coupled cluster jobs finish correctly?')
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
                    QApplication.processEvents() 

                    return None, None
        
        #Calculate Gibbs energy without coupled cluster energies
        else:            
            try:
                #Read in gibbs corrections and check if they're numeric by getting the mean - will throw a type error if there are any non numeric entries present (like N/A)
                #get gibbs correction instead of electronic energy because the thermochem code always prints out the final electronic energy even if a job didnt finish. Gibbs energies only get written to the .csv if the job finishes.
                gibbs_corr = dft_data.get('Gibbs Correction', 'N/A')
                gibbs_corr_check = gibbs_corr.mean()

                output_data['Electronic energy'] = dft_data['Electronic energy'] 
                output_data['Gibbs energy'] = dft_data['Electronic energy'] + dft_data['Gibbs Correction']

                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} DFT energies and Gibbs corrections were successfully read in from {os.path.basename(dft_file_path)}.')
                QApplication.processEvents() 

            except TypeError:
                #let the user know that both CCSD(T) and CCSD entries are non-numeric and cannot be used - the file needs to be checked for correctness.
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Gibbs corrections in {os.path.basename(dft_file_path)} were non-numeric. Did the DFT jobs finish correctly?')
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
                QApplication.processEvents() 

                return None, None

        #Calculate relative Gibbs Energy in kJ/mol if energies and Gibbs corrections were extracted correctly
        output_data['Relative Gibbs Energy (kJ/mol)'] = (output_data['Gibbs energy'] - output_data['Gibbs energy'].min()) * 2625.5
        
        #Calculate total and relative populations, ensuring that the sum of relative populations equations to unity
        output_data['Population'] = np.exp(-output_data['Relative Gibbs Energy (kJ/mol)'] / (R * temperature))
        total_population = output_data['Population'].sum()
        output_data['Relative Population'] = output_data['Population'] / total_population

        #Ensure the sum of all relative populations equals 1. If not, throw an error.
        if not np.isclose(output_data['Relative Population'].sum(), 1.0, atol=1E-4):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Sum of all relative populations does not equal 1; please check the data you are providing to the code.')
            print(f'Total population: {total_population}')
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Calculation of the BW CCS will now be aborted.')
            return None, None
        
        #Add the CCS data to the dataframe, then calculate the BW CCS
        output_data['CCS (A**2)'] = ccs_data['CCS [A**2]']
        output_data['errCCS (A**2)'] = ccs_data['errCCS [A**2]']

        #Boltzmann weighting of the CCS
        BW_CCS = np.round((output_data['Relative Population'] * output_data['CCS (A**2)']).sum(), 2)
    
        return output_data, BW_CCS

    def calculate_BW_CCS_stdev(output_data):
        
        #The variance of a weighted sum can be calculated as the sum of the squares of the products of the weights and their corresponding variances. 
        #This works b/c the population values act as weights in the calculation of BW_CCS.
        variance_BW_CCS = np.sum(np.square(output_data['Relative Population'] * output_data['errCCS (A**2)']))
        stdev_bw_ccs = np.round(np.sqrt(variance_BW_CCS), 2)

        return stdev_bw_ccs  

    def write_to_excel(output_data, output_file_path, BW_CCS, BW_CCS_stdev):
        #Create a new column for BW_CCS that only contains the value once
        output_data['Boltzmann Weighted CCS (A**2)'] = ''
        output_data.loc[0, 'Boltzmann Weighted CCS (A**2)'] = BW_CCS 

        output_data['Boltzmann Weighted CCS stdev (A**2)'] = ''
        output_data.loc[0, 'Boltzmann Weighted CCS stdev (A**2)'] = BW_CCS_stdev 

        #Write to excel
        with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
            output_data.to_excel(writer, sheet_name='Processed Data', index=False)

    #Execute the function
    output_data, BW_CCS = calculate_BW_CCS(dft_file_path, dlpno_file_path, ccs_file_path, temp) #function returns None for output_data and BW_CCS on any error

    #If a BW CCS was generated, proceed to calculate the st dev and write all data to the output excel file. If a BW CCS was not generated, abort any further operation. 
    if output_data is not None:
        BW_CCS_stdev = calculate_BW_CCS_stdev(output_data)
        write_to_excel(output_data, output_file_path, BW_CCS, BW_CCS_stdev)
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} BW CCS calculated successfully: {BW_CCS} +/- {BW_CCS_stdev} A**2')
        return

    else:
        return
    
if __name__ == '__main__':
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T11_BW_CCS_Analyzer\File_sources'
    dft_file_path = os.path.join(directory, 'DFT', 'Thermo_data_298K_100kPa_1-0vibscl.csv')
    dlpno_file_path = os.path.join(directory, 'CCSDT', 'Coupled_cluster_energies.csv') #set to none if not provided
    #dlpno_file_path = None
    ccs_file_path = os.path.join(directory, 'CCS', 'Export_mout_CCS_lowfield.csv')
    output_file_path = os.path.join(directory, 'BW_CCS_CCSDT.xlsx')

    temp = 298.15

    BW_CCS_Analysis(dft_file_path, dlpno_file_path, ccs_file_path, output_file_path, temp)
