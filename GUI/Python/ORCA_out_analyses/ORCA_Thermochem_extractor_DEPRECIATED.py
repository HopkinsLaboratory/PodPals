import os, re, time
import numpy as np
import pandas as pd
from datetime import datetime

def ORCA_Thermochem_Extractor(directory):

    start = time.time()

    thermo_properties = [
        'Filename',
        'Imaginary Frequencies',
        'Electronic energy',
        'ZPE correction',
        'Thermal correction',
        'Enthalpy correction',
        'Gibbs Correction',
        'Total ZPE',
        'Total Thermal Energy',
        'Total Enthalpy',
        'Total Entropy',
        'Total Gibbs Energy',
        'Relative Gibbs Energy',
    ]

    # Format the header for consistent spacing 
    header = '{},\n'.format(','.join(['{:<25}'] * len(thermo_properties)))

    # Create output file and write header to it
    output_csv = os.path.join(directory, 'Thermo_data.csv')

    try:
        with open(output_csv, 'w') as opf:
            opf.write(header.format(*thermo_properties))

    except (PermissionError, FileExistsError):
        
        i = 1
        if FileExistsError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} A file with the same name already exists. Writing data to a new file name.')

        if PermissionError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} A file with the same name is already open. Writing data to a new file name.')

        while os.path.isfile(output_csv):
            output_csv = os.path.join(directory, f'Thermo_data_{i}.csv')
            i += 1
        
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The new output file is {os.path.basename(output_csv)}')

        with open(output_csv, 'w') as opf:
            opf.write(header.format(*thermo_properties))
            
    # Function to extract a numerical property from the ORCA output file
    def extract_property(data, pattern):
        try:
            result = float(re.findall(pattern, data)[-1].strip())
        except:
            result = 12345.0  # Placeholder for missing values
        return result

    Gibbs_list = []
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x.lower()]
    missing_thermochem = []
    imag_freqs = []

    for filename in filenames:

        with open(os.path.join(directory, filename), 'r') as opf:
            data = opf.read()

        # Extract thermochemical properties from the ORCA output file
        E_el = extract_property(data, 'FINAL SINGLE POINT ENERGY(.*?)\n')
        ZPE_corr = extract_property(data, 'Zero point energy                ...(.*?)Eh')
        Thermal_corr = extract_property(data, 'Total thermal correction(.*?)Eh')
        H_corr = extract_property(data, 'Thermal Enthalpy correction       ...(.*?)Eh')
        S_tot = extract_property(data, 'Final entropy term                ...(.*?)Eh')
        G_corr = extract_property(data, r'G-E\(el\)                           ...(.*?)Eh')
        E_thermal = extract_property(data, 'Total thermal energy(.*?)Eh')
        E_Enthalpy = extract_property(data, 'Total Enthalpy                    ...(.*?)Eh')
        E_Gibbs = extract_property(data, 'Final Gibbs free energy         ...(.*?)Eh')
        E_ZPE = E_el + ZPE_corr
        n_imag = data.count('***imaginary mode***')      
        Erel = 123.0

        Gibbs_list.append(E_Gibbs)

        if 12345.0 in [E_el, ZPE_corr, Thermal_corr, H_corr, G_corr, E_ZPE, E_thermal, E_Enthalpy, S_tot, E_Gibbs]:
            #print(f'{filename} is missing thermochemistry. Writing 12345 as a placeholder for missing value')
            missing_thermochem.append(filename)

        #if a file contians imaginary frequencies, extract their values and write them to a list for printing laster on. This is useful for determining whether a DFT job needs to be resubmitted.
        if n_imag > 0:
            
            try:
                freq_text_pattern = r'([-+]?\d+\.\d+\s+cm\*\*-1)\s+\*\*\*imaginary mode\*\*\*' #regex to pull out <wavenumber>cm**-1 from text 
                imag_freq_matches = re.findall(freq_text_pattern, data)
               
                imag_freqs.append([filename, imag_freq_matches]) #append list of imag freqs to master imag_freqs alongside the associated filename
                
            except Exception as e:
                print(f'Error extracting imaginary frequency data from {filename}: {e}')
                pass

        # Prepare the values to be written
        values = [filename, n_imag, E_el, ZPE_corr, Thermal_corr, H_corr, G_corr, E_ZPE, E_thermal, E_Enthalpy, S_tot, E_Gibbs, Erel]

        # Create a format string for consistent spacing
        format_str = '{}'.format(','.join(['{:<25}'] * len(values)) + ',\n')

        # Append the extracted properties to the CSV file
        with open(output_csv, 'a') as opf:
            opf.write(format_str.format(*values))

    # Calculate the minimum Gibbs energy
    min_Gibbs = np.min(Gibbs_list)

    # Read the CSV into a pandas DataFrame
    df = pd.read_csv(output_csv)

    # Calculate relative energy column and update it in the DataFrame
    df['Relative Gibbs Energy    '] = (df['Total Gibbs Energy       '] - min_Gibbs) * 2625.5

    # Sort the DataFrame based on 'Relative Gibbs Energy'
    df = df.sort_values(by='Relative Gibbs Energy    ')

    # Write the updated DataFrame back to the CSV with consistent spacing
    with open(output_csv, 'w') as opf:
        opf.write(header.format(*thermo_properties))  # Write the header first

        for _, row in df.iterrows():
            values = list(row)
            opf.write(format_str.format(*values))

    format_missing_thermochem = str('\n'.join(missing_thermochem))

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Electronic energies and thermochemistry from {len(filenames)} ORCA .out files have been processed in {np.round(time.time() - start,2)} seconds.')
    
    if len(missing_thermochem) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(missing_thermochem)} files are missing thermochemistry. -12345.0 is being written as a placeholder for the following files:\n{format_missing_thermochem}\n')

    if len(imag_freqs) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(imag_freqs)} file(s) contain imaginary frequencies:')
        for filename, imag_freq_matches in imag_freqs:
            imag_freq_str = ', '.join(map(str, imag_freq_matches)) #imag freqs need to be a string in order to be printed to terminal
            print(f'{filename}: {imag_freq_str}')

if __name__ == "__main__":

    directory = r'G:\Hopkins_Laboratory\Completed Projects\Preferential_Solvation\8NH2_2Me_Quinoline\Amino_prot\MeCN\3MeCN\DFT_2'
    ORCA_Thermochem_Extractor(directory)