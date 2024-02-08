import os, time
import numpy as np
import pandas as pd
from datetime import datetime
from PyQt6.QtWidgets import QApplication

def xyz_file_splitter(file,basename):
    
    #Error handling for invalid files and/or files without a .xyz extension
    if not file.endswith('.xyz') or not os.path.isfile(file):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {os.path.basename(file)} is not a .xyz file. Please provide a file with the correct extension.')
        return

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Geometries from {os.path.basename(file)} are now being extracted. Output files will be written to the /Conformers folder with the basename {basename}.')
    QApplication.processEvents()

    start = time.time()
    
    #define prerequisites 
    output_path = os.path.join(os.path.dirname(file), 'Conformers')
    energies_file = os.path.join(output_path, 'Energies.csv')

    #initialize an empty pandas df to store filename and energies
    df_energies = pd.DataFrame(columns=['Filename', 'Energy / hartree'])

    #Read the data from the input file
    with open(file, 'r') as opf:
        data = opf.read()

    #Creating the output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    #Split the data into conformers by looking for the two lines containing the number of atoms and the energy of the conformer
    split_by = str(data.split('\n')[0].strip()) + '\n' + str(data.split('\n')[1].split('-')[0])
    conformers = data.split(split_by)[1:] #first entry will be empty since the file starts with the number of atoms + energy of the conformer

    i = 1  #Set index to one for naming files

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Extracting geometries from {os.path.basename(file)} to /Conformers ...')
    QApplication.processEvents()

    for conformer in conformers:
        conformer_lines = conformer.split('\n')
        filename_prefix = basename if not basename.endswith('_') else basename[:-1] #ensure basename ends with an underscore for readability of resulting filenames

        #Write the .gjf file
        gjf_file_path = os.path.join(output_path, f'{filename_prefix}_{i}.gjf')
        
        with open(gjf_file_path, 'w') as opf:
            opf.writelines(['#opt\n\n', f'{basename}\n\n', '1 1\n'])  #Header with default charge and multiplicity
            opf.writelines('\n'.join(conformer_lines[1:]) + '\n\n')  #Conformer geometry

        #Append filename and energy data to the pandas df
        df_energies = pd.concat([df_energies, pd.DataFrame([{'Filename': f'{filename_prefix}_{i}.gjf', 'Energy / hartree': conformer_lines[0].strip()}])], ignore_index=True)

        i += 1  #Increment file index

    #Ensure the data is sorted by energy, then calculate the relative energies
    df_energies['Energy / hartree'] = pd.to_numeric(df_energies['Energy / hartree'], errors='coerce')
    df_energies = df_energies.sort_values(by='Energy / hartree')
    E_min = df_energies['Energy / hartree'].min()
    df_energies['Relative Energy / kJ mol**-1'] = (df_energies['Energy / hartree'] - E_min) * 2625.5

    #write the Energies.csv
    df_energies.to_csv(energies_file, index=False)

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {i} geometries were extracted from {os.path.basename(file)} and written to .gjf files in {np.round(time.time() - start, 2)} seconds.\n')
    return
    
#external testing
if __name__ == '__main__':
    file = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T1_CREST_xyz_splitter\crest_conformers.xyz'
    basename = 'Fluoxetine'

    xyz_file_splitter(file,basename)