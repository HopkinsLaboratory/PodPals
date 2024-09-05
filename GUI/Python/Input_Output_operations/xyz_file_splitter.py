import os, re, time
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

    def split_xyz_file(data):
        '''Extract geometries from an xyz file containing multiple conformers into individual files.
        Returns a list of lists, where each sub-list is the line of the xyz file'''

        lines = data.strip().splitlines()

        #First line is the number of atoms 
        num_atoms = int(lines[0].strip())

        #then, the size of each block is always num_atoms + 2 because an xyz file format is:
        # num_atoms
        # title car
        # geom data

        chunk_size = num_atoms + 2

        #initialize lists for holding conformers
        conformers = []
        current_conformer = []
        line_count = 0

        # Go through line by line, appending each line to the current conformer list until the chunck size is reached
        for line in lines:
            if line_count == 1:
                current_conformer.append(line)  #Always include the second line (title card), even if blank
            
            elif line.strip() or line_count < chunk_size: #exclude all other blank lines, append all remaining lines to the current conformer list
                current_conformer.append(line.strip())

            line_count += 1

            #Once we reach the size of the xyz file (2 lines plus the number of atoms), save the conformer 
            if line_count == chunk_size:
                
                #There can be empty lines after the xyz data and before the next xyz file starts, so we need to remove these
                filtered_conformer = [current_conformer[0]]  #Always include the first line
                filtered_conformer.append(current_conformer[1])  #Always include the second line (even if blank)
                filtered_conformer.extend([l for l in current_conformer[2:] if l.strip()]) 

                conformers.append('\n'.join(filtered_conformer))

                #Reinitialize variables for the next conformer
                current_conformer = []
                line_count = 0                
        
        #Handle the last conformer if the file doesn't end perfectly on a chunk boundary (natoms + 2)
        if current_conformer:
            filtered_conformer = [current_conformer[0]]  # Always include the first line
            filtered_conformer.append(current_conformer[1])  # Always include the second line (even if blank)
            filtered_conformer.extend([l for l in current_conformer[2:] if l.strip()])
            conformers.append('\n'.join(filtered_conformer))
        
        return conformers
    
    def extract_energy(line):
        '''Regex to match a number, which may include a negative sign and/or decimal point
        Written for the purpose of extracting the energy from the second line of an xyz file. 
        Returns the number if found, or None if not found.
        '''
        
        match = re.match(r'^-?\d+(\.\d+)?', line.strip())
        if match:
            return match.group(0)
        return None

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

    #Split the geometries in the xyz file into conformers 
    conformers = split_xyz_file(data)

    if len(conformers) < 2:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} WARNING: Only one conformer was found in {os.path.basename(file)} - is this what you expected?')
        QApplication.processEvents()

    i = 1  #Set index to one for naming files
    missing_energy = [] #list to write conformer numbers that are missing energy

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Extracting geometries from {os.path.basename(file)} to /Conformers ...')
    QApplication.processEvents()

    for conformer in conformers:
        lines = conformer.split('\n')
        
        #try to extract the energy from the second line of conformer
        energy = extract_energy(lines[1].strip()) if len(lines) > 1 else None

        #if energy is not found, append the conformer number to the missing energy list to let the user know later. Also assign a placeholder value of 12345.0 for writing to file
        if not energy:
            missing_energy.append(i)
            energy = 12345.0

        #the only entry with 4 splits, where instance 1 is an atom symbol (sometimes an atom number!), and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
        geom_lines = [line for line in lines if len(line.split()) == 4 and (line.split()[0].isalpha() or line.split()[0].isdigit()) and all('.' in x for x in line.split()[1:3])]
        
        #Format the geom line for aesthetics
        formatted_geom_lines = []
        for line in geom_lines:
            parts = line.split()
            label = parts[0].ljust(4)
            x = f"{float(parts[1]):>{15}.8f}"  
            y = f"{float(parts[2]):>{15}.8f}"  
            z = f"{float(parts[3]):>{15}.8f}"  
            formatted_geom_lines.append(f"{label}{x}{y}{z}\n")
        
        filename_prefix = basename if not basename.endswith('_') else basename[:-1] #ensure basename ends with an underscore for readability of resulting filenames

        #Write the .gjf file
        export_file_path = os.path.join(output_path, f'{filename_prefix}_{i}.{export_type}')
        
        #write to specified file type
        with open(export_file_path, 'w') as opf:
            if export_type == 'gjf':
                opf.writelines(['#opt pm7\n\n', f'{filename_prefix}_{i}\n\n', '1 1\n'])  #Header with default charge and multiplicity
                opf.writelines(formatted_geom_lines)
                opf.write('\n\n\n')  #Conformer geometry
            
            elif export_type == 'xyz':
                opf.write(f'{len(geom_lines)}\n{filename_prefix}_{i} {energy}\n')  #Header with default charge and multiplicity
                opf.write(formatted_geom_lines)
                opf.write('\n\n')  #Conformer geometry                

            elif export_type == 'inp':
                opf.writelines(['! TPSS D3BJ TightOpt Def2-SVP\n', '%maxcore 3000\n\n', '%pal nprocs 8\nend\n\n', '*xyz 1 1\n'])  #Header with default charge and multiplicity
                opf.writelines(formatted_geom_lines) #Conformer geometry
                opf.write('*\n\n\n')  

        #Append filename and energy data to the pandas df
        df_energies = pd.concat([df_energies, pd.DataFrame([{'Filename': f'{filename_prefix}_{i}.{export_type}', 'Energy / hartree': energy}])], ignore_index=True)

        i += 1  #Increment file index

    if missing_energy:
        
        #convert missing energy list to str
        missing_energy_str = ', '.join(map(str, missing_energy))
        
        if len(missing_energy) == (i - 1):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Energies could not be found for any conformer in the .xyz file. Energies of 12345.0 have been written as placeholders to Energies.csv.')
        
        else:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Energies could not be found in the following conformers - energies of 12345.0 have been written as placeholders to Energies.csv: {missing_energy_str}')
        
    #Ensure the data is sorted by energy, then calculate the relative energies
    df_energies['Energy / hartree'] = pd.to_numeric(df_energies['Energy / hartree'], errors='coerce')
    df_energies = df_energies.sort_values(by='Energy / hartree')
    E_min = df_energies['Energy / hartree'].min()
    df_energies['Relative Energy / kJ mol**-1'] = (df_energies['Energy / hartree'] - E_min) * 2625.5

    #write the Energies.csv
    df_energies.to_csv(energies_file, index=False)

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {i-1} geometries were extracted from {os.path.basename(file)} and written to .{export_type} files in {np.round(time.time() - start, 2)} seconds.\n')
    return
    
#external testing
if __name__ == '__main__':
    file = r'E:\OneDrive - University of Waterloo\Waterloo\Manuscripts\2024\Leipzig_PIC\Calculated_IR_Spectra\Ser\RS-Ser\GOAT\RS_Ser_GOAT_0deg.finalensemble.xyz'
    basename = 'RS_Ser_GOAT_0deg'
    export_type = 'inp'

    xyz_file_splitter(file, basename, export_type)