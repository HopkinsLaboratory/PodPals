# Created by CI 2021-05-17
# Converts Gaussian input (.gjf) files to ORCA input files

import os
import time
import numpy as np
from datetime import datetime

def Gaussian_gjf_to_ORCA_input(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked):

    # Generate a directory to write new files to
    new_dir = os.path.join(directory, 'New_Inputs')
    os.makedirs(new_dir, exist_ok=True)  # Create the directory if it doesn't exist

    # Generate a list of .gjf files from the directory
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.gjf')]
    num_files = len(filenames)

    start = time.time()

    #remove any leadng or trailing whitespace from the calc line
    calc_line = calc_line.strip()

    #Error handling for calc_line
    # check if calc line begns wth a !
    if not calc_line.startswith('!'):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The calc line does not start with a !')
        calc_line = f'! {calc_line}'
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Fixing it now; the new calc line is:\n {calc_line}')
    
    #Check for CHELPG or chelpg in the calc line if ESP charges are requested
    if esp_charges_checked and not 'chelpg' in calc_line.lower():
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} You have requested that ESP charges are calculcated, but CHELPG is not in the calulcation line. Adding it now.')
        calc_line = f'{calc_line} CHELPG'
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The new calc line is:\n {calc_line}')

    #likewise, if the calc line contains the CHELPG keyword and ESP charges are not requested, we need to assign the custom grid for accurate partial charges
    if 'chelpg' in calc_line.lower() and not esp_charges_checked:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} CHELPG was found in the calc_line but you have not requested that ESP charges are calculcated.')
        esp_charges_checked = True
        grid = 0.1
        rmax = 3.0
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The ESP calculation block will be written to the .inp file with a grid size of {grid} angrstroms and a rmax of {rmax} angstroms.')
    
    #check that the memory allocated per core is an integer
    if not isinstance(mpp, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of memory is not an integer. You cannot have fractional CPUs! Fixing it now')
        mpp = round(mpp)
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The mpp being written to each .inp is {mpp}.')

    #check that the number of cores requested for the job is an integer
    if not isinstance(ncores, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores specified is not an integer. You cannot have fractional CPUs! Fixing it now')
        ncores = round(ncores)   
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores being written to each .inp is {ncores}.')

    # Create ORCA files
    for filename in filenames:
        orca_filename = os.path.join(new_dir, f'{filename[:-4]}.inp')
        
        with open(orca_filename, 'w') as opf:
            # Input block
            opf.write(f'{calc_line}\n')

            # Memory
            opf.write(f'%maxcore {mpp}\n\n')

            # Number of cores
            opf.write(f'%pal nprocs {ncores} \nend\n\n')

            # ESP Charges
            if esp_charges_checked:
                opf.write(f'%chelpg\n')
                opf.write(f'grid {grid}\n')
                opf.write(f'rmax {rmax}\n')
                opf.write('end\n\n')

            # Open shell calcs
            if multiplicity != 1:
                opf.write(f'%scf HFTyp UHF\n')
                opf.write('end\n\n')

            # Calc Hessian during first step of the optimization (if requested)
            if calc_hess_checked:
                opf.write(f'%geom\n')
                opf.write('Calc_Hess true\n')
                opf.write('end\n\n')           
            
            # Calculate dipole and quadrupole moments (if requested)
            if polarization_checked:
                opf.write(f'%elprop\n')
                opf.write('Dipole true\n')
                opf.write('Quadrupole true\n')
                opf.write('Polar 1\n')
                opf.write('end\n\n') 

        #read and extract geometry from the .gjf files
        with open(os.path.join(directory, filename), 'r') as opf:
            lines = opf.readlines()

        geometry = []  # Get geometry of atoms from file
        for line in lines:
            split_line = line.split() #split by whitespace
            #the only entry with 4 splits and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
            if len(split_line) == 4 and all('.' in x for x in split_line[1:3]):
                geometry.append(line.split())

        #write geom to a formatted string for printing to a file
        fs = '\n'.join([f'{i[0]:<5s}  {i[1]:>15s}  {i[2]:>15s}  {i[3]:>15s}' for i in geometry])

        #user can call a geometry from an external xyz file within the .inp (if requested)
        if write_xyz_checked:
            xyz_file = os.path.join(new_dir, f'{filename[:-4]}.xyz')
            with open(xyz_file, 'w') as opf:
                opf.write(f'{len(geometry)}\n\n{fs}\n\n')
            with open(orca_filename, 'a') as opf:
                opf.write(f'*xyzfile {charge} {multiplicity} {os.path.basename(xyz_file)}\n')
        
        #otherwise, write the geometry to the .inp
        else:
            with open(orca_filename, 'a') as opf:
                opf.write(f'*xyz {charge} {multiplicity}\n{fs}\n*\n\n')
    
    #write a completion message to the GUI output window
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {num_files} .gjf files were converted to ORCA .inp files in {np.round(time.time() - start,2)} seconds.')
