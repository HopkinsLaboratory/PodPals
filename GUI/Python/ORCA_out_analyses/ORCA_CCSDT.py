import os, re, time
from datetime import datetime
import numpy as np
from PyQt6.QtWidgets import QApplication

def ORCA_CCSDT(directory):

    start = time.time()

    #get list of filenames, and check if the directory does not contain any .out files
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x.lower()]

    if len(filenames) == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} There are no .out files in the provided directory.')
        return
        
    properties = [
        'Filename',
        'DLPNO-CCSD(T) energy',
    ]

    # Format the header for consistent spacing 
    header = '{}\n'.format(','.join(['{:<25}'] * len(properties)))

    # Create output file and write header to it, ensuring that previous files of the same name are not overwritten
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
    
    missing_CCSDT = [] #list to write files that are missing CCSDT energies

    for filename in filenames:
        with open(os.path.join(directory, filename), 'r') as opf:
            data = opf.read()

        # Extract thermochemical properties from the ORCA output file
        try: 
            CCSDT = float(re.findall(r'E\(CCSD\(T\)\)                                 ...(.*?)\n', data)[-1].strip())

        except:
            CCSDT = -12345.0
            missing_CCSDT.append(filename)

        # Prepare the values to be written
        values = [filename, CCSDT]

        # Create a format string for consistent spacing
        format_str = '{}'.format('{:<25},' * (len(values) - 1) + '{:<25}\n')

        # Append the extracted properties to the CSV file
        with open(output_csv, 'a') as opf:
            opf.write(format_str.format(*values))

    if len(missing_CCSDT) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(missing_CCSDT)} file(s) did not contain CCSDT energies:\n{'\n'.join(missing_CCSDT)}')

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} DLPNO-CCSD(T) energies were extracted from {len(filenames)} ORCA .out files in {np.round(time.time() - start,2)} seconds.')

