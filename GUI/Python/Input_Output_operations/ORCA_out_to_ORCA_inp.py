import os, re, time
import numpy as np
from datetime import datetime
from PyQt6.QtWidgets import QApplication

def ORCA_out_to_ORCA_inp(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked, write_gjf_checked):

    #Generate a list of .gjf files from the directory. Note that pseudopotentials write _atom##to the filename, so we filter these out. 
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x]
    num_files = len(filenames)

    if num_files == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .out files were found in {os.path.basename(directory)}.')
        return
    
    #Generate the directory/directories to write new files to
    if write_xyz_checked:
        new_dir = os.path.join(directory, 'New_Inputs_xyz')

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {num_files} ORCA .out files to ORCA .inp files.\nORCA .inp and .xyz files will be written to {os.path.basename(new_dir)} ...')
        QApplication.processEvents()
        
    else: 
        new_dir = os.path.join(directory, 'New_Inputs')
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {num_files} ORCA .out files to ORCA .inp files.\nORCA .inp and will be written to {os.path.basename(new_dir)} ...')
        QApplication.processEvents()

    os.makedirs(new_dir, exist_ok=True)  

    if write_gjf_checked:
        gjf_dir = os.path.join(new_dir, 'gjfs')
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} .gjf files will be written to {os.path.join(os.path.basename(new_dir), "gjfs")}')
        os.makedirs(gjf_dir, exist_ok=True) 

    start = time.time()

    #Error handling 

    #remove any leadng or trailing whitespace from the calc line
    calc_line = calc_line.strip()

    #check is calc line begns wth a !
    if not calc_line.startswith('!'):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The input line in {filename} does not start with a !')
        calc_line = f'!  + {calc_line}'
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
            #read and extract geometry from the .gjf files
            with open(os.path.join(directory, filename), 'r') as opf:
                data = opf.read()

        #If there is an error opening the .gjf file, infomr that user that it will be skipped.
        except IOError as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error was encountered when trying to open {filename}: {e}.\n Processing of this file will be skipped.')
            QApplication.processEvents()
            continue

        #Get the geometry and write to a list called geometry
        geometry = [] #Initialize the geometry list
        XYZ_data = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)', data)
        
        if not XYZ_data:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No atomic coordinate data was found in {filename}. Skipping...')
            QApplication.processEvents()
            continue
        
        #Loop through each line in the FINAL geom block, split each line, then append to the geometry list
        try:
            for line in XYZ_data[-1].split('\n')[2:-3]:
                split_line = line.split()
                #the only entry with 4 splits, where instance 1 is an atom symbol (sometimes an atom number!), and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
                if len(split_line) == 4 and (split_line[0].isalpha() or split_line[0].isdigit()) and all('.' in x for x in split_line[1:3]):
                    geometry.append(split_line)

        #sometimes files with a .out extension from other programs are found in the target directory - lets handle those!
        except IndexError as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Index error encountered in {filename}. Skipping...')
            QApplication.processEvents()
            continue

        #If geometry is extracted, start to write the ORCA .inp file
        else:

            #Define the ORCA .inp filename
            orca_filename = os.path.join(new_dir, f'{filename[:-4]}.inp')

            #Write the geometry to a formatted string
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
                    opf.write(f'*xyz {charge} {multiplicity}\n{fs}\n*\n\n')
        
        #make the .gjf file
        if write_gjf_checked: 
            with open(os.path.join(gjf_dir, f'{filename[:-4]}.gjf'), 'w') as opf_gjf:
                opf_gjf.write('#opt\n\n')
                opf_gjf.write(f'{filename[:-4]}\n\n') #filename excluding its extenstion 
                opf_gjf.write(f'{charge} {multiplicity}\n')
                opf_gjf.write(fs)
                opf_gjf.write('\n\n')
   
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {num_files} ORCA .out files were converted to ORCA .inp files in {np.round(time.time() - start,2)} seconds.')
    return

#external testing
if __name__ == '__main__':

    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T4_ORCA_out_to_ORCA_inp'
    mpp = 3500
    ncores = 8
    charge = 1
    multiplicity = 1
    calc_line = '! wB97X-D3 TightOpt Freq def2-TZVPP def2/J RIJCOSX TightSCF defgrid3'
    esp_charges_checked = False
    grid = 0.1
    rmax = 3.0
    calc_hess_checked = True
    polarization_checked = True
    write_xyz_checked = False
    write_gjf_checked = True

    ORCA_out_to_ORCA_inp(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked, write_gjf_checked)
