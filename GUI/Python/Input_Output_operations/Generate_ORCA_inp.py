import os, re
import time
import numpy as np
from datetime import datetime
from PyQt6.QtWidgets import QApplication

def Generate_ORCA_inp(directory, mpp, ncores, charge, multiplicity, cm_fromPrev, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked, write_gjf_checked):

    #Generate a list of files from the directory - currently accepts Gaussian .gjf, .xyz, and ORCA .inp. 
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.gjf') or x.lower().endswith('.inp') or x.lower().endswith('.xyz')]
    num_files = len(filenames)

    if num_files == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .gjf, .inp, or .xyz files were found in {os.path.basename(directory)}.')
        return

    #Generate a directory to write new files to

    if write_xyz_checked:
        new_dir = os.path.join(directory, 'New_Inputs_xyz')

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {num_files} files to ORCA .inp files.ORCA .inp and .xyz files will be written to {os.path.basename(new_dir)} ...')
        QApplication.processEvents()
        
    else: 
        new_dir = os.path.join(directory, 'New_Inputs')

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {num_files} files to ORCA .inp files.ORCA .inp files will be written to {os.path.basename(new_dir)} ...')
        QApplication.processEvents()

    os.makedirs(new_dir, exist_ok=True)  #Create the directory if it doesn't exist

    #write directory for gjf files (if it doesn't already exist)
    if write_gjf_checked:
        gjf_dir = os.path.join(new_dir, 'gjfs')
        os.makedirs(gjf_dir, exist_ok=True) 

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

    #check is user requests that charge/mult be read from a directory that contains xyz files, then inform that xyz does not contain charge/mult data
    if any(file.lower().endswith('.xyz') for file in filenames) and cm_fromPrev:
        print(
            f'{datetime.now().strftime("[ %H:%M:%S ]")} xyz files were found in the directory, '
            'which do not contain charge/mult data. As such, charge/mult will be read in from '
            'the entries in the GUI window. If alternative numbers are desired, please uncheck '
            '"Read charge/mult from file" and enter the desired charge/mult values.'
        )
        QApplication.processEvents()
        cm_fromPrev = False

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
    
        #Use charge/mult definition from ORCA out file if requested
        if cm_fromPrev:        
            
            #extract from .gjf files (Gaussian)
            if filename.lower().endswith('.gjf'):

                #regex to match a line containing exactly two integers separated by spaces
                gjf_cm_pattern = re.compile(r"^\s*(-?\d+)\s+(-?\d+)\s*$")
                
                for line in lines:
                    gjf_cm_match = gjf_cm_pattern.match(line.strip())
                    if gjf_cm_match:
                        #overwrite charge/mult variables if found in file
                        charge = int(gjf_cm_match.group(1)) 
                        multiplicity = int(gjf_cm_match.group(2))
                        break        

            #extract from .inp files (ORCA)
            if filename.lower().endswith('.inp'):
                #regex to match lines starting with '*' and containing two integers
                orcainp_cm_pattern = re.compile(r"^\*\S*\s+(-?\d+)\s+(-?\d+)", line.strip())

                for line in lines:
                    orcainp_cm_match = orcainp_cm_pattern.match(line.strip())
                    if orcainp_cm_match:
                        #overwrite charge/mult variables if found in file
                        charge = int(orcainp_cm_match.group(1))
                        multiplicity = int(orcainp_cm_match.group(2))
                        break 

        geometry = []  #Initialize the geom list

        for line in lines:
            #Remove the contents within parentheses and any trailing whitespace - this is a thing for .pdb files created with Gaussview
            if '(' or ')' in line:
                clean_line = re.sub(r'\(.*?\)', '', line).strip()
                split_line = clean_line.split()
            else:
                clean_line = line.strip()
                split_line = clean_line.split()
            
            #the only entry with 4 splits, where instance 1 is an atom symbol (sometimes an atom number!), and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
            if len(split_line) == 4 and (split_line[0].isalpha() or split_line[0].isdigit()) and all('.' in x for x in split_line[1:3]):
                geometry.append(split_line)

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

                #GOAT block - new to ORCA 6. hard-written values for now but will be made customizable in the future
                if 'goat' in calc_line.lower():
                    opf.write(f'%goat\n')
                    opf.write('ALIGN true\n')
                    #opf.write('GFNUPHILL GFNFF\n')
                    opf.write('TEMPLIST 3000, 2000, 750, 500\n')
                    opf.write('FREEZECISTRANS TRUE\n')
                    opf.write('MAXEN 5.0\n') 
                    if 'entropy' in calc_line.lower():
                        opf.write('CONFTEMP 298.15\n') 
                        opf.write('CONFDEGEN AUTO\n')                                      
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

                #user can request .gjf files be written - useful for generating .gjf inputs for visualization in Gaussview
                if write_gjf_checked: 
                    with open(os.path.join(gjf_dir, f'{filename[:-4]}.gjf'), 'w') as opf_gjf:
                        opf_gjf.write('#opt pm7\n\n')
                        opf_gjf.write(f'{filename[:-4]}\n\n') #filename excluding its extenstion 
                        opf_gjf.write(f'{charge} {multiplicity}\n')
                        opf_gjf.write(fs)
                        opf_gjf.write('\n\n')

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {num_files} .gjf files were converted to ORCA .inp files in {np.round(time.time() - start,2)} seconds.')
    return

#external testing
if __name__ == '__main__':

    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T3_gjf_to_ORCA_inp\inp_modification'
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
    write_gjf_checked = False

    Generate_ORCA_inp(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked, write_gjf_checked)