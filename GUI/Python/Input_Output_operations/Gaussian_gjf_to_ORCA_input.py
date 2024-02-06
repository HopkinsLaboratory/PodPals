import os
import time
import numpy as np
from datetime import datetime
from PyQt6.QtWidgets import QApplication

def Gaussian_gjf_to_ORCA_input(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked):

    #Generate a list of .gjf files from the directory
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.gjf')]
    num_files = len(filenames)

    #Generate a directory to write new files to

    if write_xyz_checked:
        new_dir = os.path.join(directory, 'New_Inputs_xyz')

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {num_files} .gjf files to ORCA .inp files.\nORCA .inp and .xyz files will be written to {os.path.basename(new_dir)} ...')
        QApplication.processEvents()
        
    else: 
        new_dir = os.path.join(directory, 'New_Inputs')

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {num_files} .gjf files to ORCA .inp files.\nORCA .inp files will be written to {os.path.basename(new_dir)} ...')
        QApplication.processEvents()

    os.makedirs(new_dir, exist_ok=True)  #Create the directory if it doesn't exist

    start = time.time()

    #remove any leadng or trailing whitespace from the calc line
    calc_line = calc_line.strip()

    #Error handling for calc_line
    #check if calc line begns wth a !
    if not calc_line.startswith('!'):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The calc line does not start with a !')
        calc_line = f'! {calc_line}'
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Fixing it now; the new calc line is:\n {calc_line}')
        QApplication.processEvents()
    
    #Check for CHELPG or chelpg in the calc line if ESP charges are requested
    if esp_charges_checked and not 'chelpg' in calc_line.lower():
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} You have requested that ESP charges are calculcated, but CHELPG is not in the calulcation line. Adding it now.')
        calc_line = f'{calc_line} CHELPG'
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The new calc line is:\n {calc_line}')
        QApplication.processEvents()

    #likewise, if the calc line contains the CHELPG keyword and ESP charges are not requested, we need to assign the custom grid for accurate partial charges
    if 'chelpg' in calc_line.lower() and not esp_charges_checked:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} CHELPG was found in the calc_line but you have not requested that ESP charges are calculcated.')
        esp_charges_checked = True
        grid = 0.1
        rmax = 3.0
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The ESP calculation block will be written to the .inp file with a grid size of {grid} angrstroms and a rmax of {rmax} angstroms.')
        QApplication.processEvents()
    
    #check that the memory allocated per core is an integer
    if not isinstance(mpp, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of memory is not an integer. You cannot have fractional CPUs! Fixing it now')
        mpp = round(mpp)
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The mpp being written to each .inp is {mpp}.')
        QApplication.processEvents()

    #check that the number of cores requested for the job is an integer
    if not isinstance(ncores, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores specified is not an integer. You cannot have fractional CPUs! Fixing it now')
        ncores = round(ncores)   
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores being written to each .inp is {ncores}.')
        QApplication.processEvents()

    #Create ORCA files
    for filename in filenames:
    
        #read and extract geometry from the .gjf files
        try:
            with open(os.path.join(directory, filename), 'r') as opf:
                lines = opf.readlines()

        #If there is an error opening the .gjf file, infomr that user that it will be skipped.
        except IOError as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error was encountered when trying to open {filename}: {e}.\n Processing of this file will be skipped.')
            QApplication.processEvents()
            continue
    
        geometry = []  #Initialize the geom list

        for line in lines:
            split_line = line.split() #split by whitespace
            
            #the only entry with 4 splits, where instance 1 is an atom symbol (sometimes an atom number!), and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
            if len(split_line) == 4 and (split_line[0].isalpha() or split_line[0].isdigit()) and all('.' in x for x in split_line[1:3]):
                geometry.append(line.split())

        if not geometry:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No atomic coordinate data was found in {filename}. Processing of this file will be skipped.')
            QApplication.processEvents()
            continue
        
        #If geometry is extracted, start to write the ORCA .inp file
        else:
            
            #Define the ORCA .inp filename
            orca_filename = os.path.join(new_dir, f'{filename[:-4]}.inp')

            #write geom to a formatted string for printing to a file
            fs = '\n'.join([f'{i[0]:<5s}  {i[1]:>15s}  {i[2]:>15s}  {i[3]:>15s}' for i in geometry])    

            #Write the ORCA .inp file
            with open(orca_filename, 'w') as opf:
                #Input block
                opf.write(f'{calc_line}\n')

                #Memory
                opf.write(f'%maxcore {mpp}\n\n')

                #Number of cores
                opf.write(f'%pal nprocs {ncores} \nend\n\n')

                #ESP Charges
                if esp_charges_checked:
                    opf.write(f'%chelpg\n')
                    opf.write(f'grid {grid}\n')
                    opf.write(f'rmax {rmax}\n')
                    opf.write('end\n\n')

                #Open shell calcs
                if multiplicity != 1:
                    opf.write(f'%scf HFTyp UHF\n')
                    opf.write('end\n\n')

                #Calc Hessian during first step of the optimization (if requested)
                if calc_hess_checked:
                    opf.write(f'%geom\n')
                    opf.write('Calc_Hess true\n')
                    opf.write('end\n\n')           
                
                #Calculate dipole and quadrupole moments (if requested)
                if polarization_checked:
                    opf.write(f'%elprop\n')
                    opf.write('Dipole true\n')
                    opf.write('Quadrupole true\n')
                    opf.write('Polar 1\n')
                    opf.write('end\n\n') 

                #user can call a geometry from an external xyz file within the .inp (if requested)
                if write_xyz_checked:
                    #Create the .xyz file containing the geometry and add a reference to it in the ORCA .inp file
                    xyz_file = os.path.join(new_dir, f'{filename[:-4]}.xyz')
                    opf.write(f'*xyzfile {charge} {multiplicity} {os.path.basename(xyz_file)}\n')

                    with open(xyz_file, 'w') as opf2xyz:
                        opf2xyz.write(f'{len(geometry)}\n\n{fs}\n\n')
                        
                #otherwise, write the geometry to the .inp
                else:
                    with open(orca_filename, 'a') as opf:
                        opf.write(f'*xyz {charge} {multiplicity}\n{fs}\n*\n\n')

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {num_files} .gjf files were converted to ORCA .inp files in {np.round(time.time() - start,2)} seconds.')
    return
