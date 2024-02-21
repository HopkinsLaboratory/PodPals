import os, time
import numpy as np
import pandas as pd
from datetime import datetime
from PyQt6.QtWidgets import QApplication

def xyz_file_splitter(file, basename, export_type):
    
    #Error handling for invalid files and/or files without a .xyz extension
    if not file.endswith('.xyz'):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {os.path.basename(file)} is not a .xyz file. Please provide a file with the correct extension.')
        return
    
    #Error handling for invalid files and/or files without a .xyz extension
    if not os.path.isfile(file):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {os.path.basename(file)} does not exist in the path specified.')
        return    
    
    #Error handling to check that the requesting format to write the extracted geoms to is supported
    if not export_type in ['inp', 'gjf', 'xyz']:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An invalid/unsupported file extension has been requested to write the extracted CREST geometries to.')
        return
    
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Geometries from {os.path.basename(file)} are now being extracted. Output files will be written to the /Conformers folder with the basename {basename}.')
    QApplication.processEvents()

    start = time.time()
    
    #define prerequisites 
    output_path = os.path.join(os.path.dirname(file), f'Conformers_{export_type}')
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

    if len(conformers) < 2:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} WARNING: Only one conformer was found in {os.path.basename(file)} - is this what you expected?')
        QApplication.processEvents()

    i = 1  #Set index to one for naming files

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Extracting geometries from {os.path.basename(file)} to /Conformers ...')
    QApplication.processEvents()

    for conformer in conformers:
        lines = conformer.split('\n')
        
        #first line is always energy
        energy = lines[0].strip()

        #the only entry with 4 splits, where instance 1 is an atom symbol (sometimes an atom number!), and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
        geom_lines = [line for line in lines if len(line.split()) == 4 and (line.split()[0].isalpha() or line.split()[0].isdigit()) and all('.' in x for x in line.split()[1:3])]
        filename_prefix = basename if not basename.endswith('_') else basename[:-1] #ensure basename ends with an underscore for readability of resulting filenames

        #Write the .gjf file
        export_file_path = os.path.join(output_path, f'{filename_prefix}_{i}.{export_type}')
        
        #write to specified file type
        with open(export_file_path, 'w') as opf:
            if export_type == 'gjf':
                opf.writelines(['#opt pm7\n\n', f'{filename_prefix}_{i}\n\n', '1 1\n'])  #Header with default charge and multiplicity
                opf.writelines('\n'.join(geom_lines) + '\n\n')  #Conformer geometry
            
            elif export_type == 'xyz':
                opf.write(f'{len(geom_lines)}\n{filename_prefix}_{i}\n')  #Header with default charge and multiplicity
                opf.writelines('\n'.join(geom_lines) + '\n\n')  #Conformer geometry

            elif export_type == 'inp':
                opf.writelines(['! PM3 Opt Def2-SVP\n', '%maxcore 1000\n\n', '%pal nprocs 1\nend\n\n', '* xyz 1 1\n'])  #Header with default charge and multiplicity
                opf.writelines('\n'.join(geom_lines) + '\n*\n\n\n')  #Conformer geometry

        #Append filename and energy data to the pandas df
        df_energies = pd.concat([df_energies, pd.DataFrame([{'Filename': f'{filename_prefix}_{i}.{export_type}', 'Energy / hartree': energy}])], ignore_index=True)

        i += 1  #Increment file index

    #Ensure the data is sorted by energy, then calculate the relative energies
    df_energies['Energy / hartree'] = pd.to_numeric(df_energies['Energy / hartree'], errors='coerce')
    df_energies = df_energies.sort_values(by='Energy / hartree')
    E_min = df_energies['Energy / hartree'].min()
    df_energies['Relative Energy / kJ mol**-1'] = (df_energies['Energy / hartree'] - E_min) * 2625.5

    #write the Energies.csv
    df_energies.to_csv(energies_file, index=False)

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {i} geometries were extracted from {os.path.basename(file)} and written to .{export_type} files in {np.round(time.time() - start, 2)} seconds.\n')
    return
    
#external testing
if __name__ == '__main__':
    file = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T1_CREST_xyz_splitter\crest_conformers.xyz'
    basename = 'Fluoxetine'
    export_type = 'xyz'

    xyz_file_splitter(file, basename, export_type)