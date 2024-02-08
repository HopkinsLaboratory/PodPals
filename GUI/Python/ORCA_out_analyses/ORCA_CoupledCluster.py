import os, re, time
from datetime import datetime
import numpy as np
from PyQt6.QtWidgets import QApplication

def extract_ORCA_coupled_cluster(directory):

    start = time.time()

    #get list of filenames, and check if the directory does not contain any .out files
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x.lower()]

    if len(filenames) == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} There are no .out files in the provided directory.')
        return
        
    properties = [
        'Filename',
        'CCSD(T) energy',
        'CCSD energy',
        'Triples correction',
        'Final correlation energy',
        'T1 diagnostic'
    ]

    #Format the header for consistent spacing 
    header = '{}\n'.format(','.join(['{:<25}'] * len(properties)))

    #Create output file and write header to it, ensuring that previous files of the same name are not overwritten
    output_csv = os.path.join(directory, 'DLPNO_CCSDT_energies.csv')

    i = 2
    while os.path.isfile(output_csv):
        output_csv = os.path.join(directory, f'DLPNO_CCSDT_energies_v{i}.csv')
        i += 1
    
    try:
        with open(output_csv, 'w') as opf:
            opf.write(header.format(*properties))
    
    except IOError as e:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error was encountered when writing to {os.path.basename(output_csv)}: {e}.\nFile processing will not proceed.')
        return
    
    missing_CCSD = [] #list to write files that are missing coupled cluster energies - CCSD will always be printed when doing coupled cluster calcs in ORCA (even when doing CCSD(T))

    #Initialize variables 
    T1_diagnostic = CCSD = triples_corr = final_corr_energy = CCSDT =  None

    for filename in filenames:
        with open(os.path.join(directory, filename), 'r') as opf:
            data = opf.readlines()

        for line in data:
            line = line.strip() #remove any pesky leading and trailing whitespace
            
            #T1 diagnostic from coupled cluster calculation
            if line.startswith('T1 diagnostic'):
                T1_diagnostic = line.split()[3]
                T1_diagnostic_float = float(T1_diagnostic)

                if T1_diagnostic_float > 0.02:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} WARNING: The T1 diagnostic in {filename} is {np.round(T1_diagnostic, 4)}, which is greater than the "safe" threshold of 0.02. The results are not to be trusted and the HF reference might be poor. For more details, consult the ORCA manual.') 

            #CCSD energy
            elif line.startswith('E(TOT)'):
                CCSD = line.split()[2]

            #Triples correction
            elif line.startswith('Triples Correction (T)'):
                triples_corr = line.split()[4]
            
            #Final correlation energy
            elif line.startswith('Final correlation energy'):
                final_corr_energy = line.split()[4]

            #CCSDT energy
            elif line.startswith('E(CCSD(T))'):
                CCSDT = line.split()[2]

        #Append filename to missing-CCSD list if no coupled cluster energy was found
        if CCSD is None:
            missing_CCSD.append(filename)
        
        #Prepare the values to be written
        values = [
            filename,
            CCSDT if CCSDT is not None else 'N/A',
            CCSD if CCSD is not None else 'N/A',
            triples_corr if triples_corr is not None else 'N/A',
            final_corr_energy if final_corr_energy is not None else 'N/A',
            T1_diagnostic if T1_diagnostic is not None else 'N/A'
        ]

        #Create a format string for consistent spacing
        format_str = '{}'.format('{:<25},' * (len(values) - 1) + '{:<25}\n')

        #Append the extracted properties to the CSV file
        with open(output_csv, 'a') as opf:
            opf.write(format_str.format(*values))

    if len(missing_CCSD) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(missing_CCSD)} file(s) did not contain coupled cluster energies:\n{", ".join(missing_CCSD)}\nDid they finish properly?')

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Coupled cluster energies were extracted from {len(filenames)} ORCA .out files in {np.round(time.time() - start,2)} seconds.')

#external testing
if __name__ == '__main__':
    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T8_ORCA_CCSDT_Analyzer'

    extract_ORCA_coupled_cluster(directory)