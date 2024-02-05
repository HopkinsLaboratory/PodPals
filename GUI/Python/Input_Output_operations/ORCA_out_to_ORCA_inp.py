import os, re, time
import numpy as np
from datetime import datetime

def ORCA_out_to_ORCA_inp(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked, write_gjf_checked):

    # Generate a directory to write new files to
    new_dir = os.path.join(directory, 'New_Inputs')
    os.makedirs(new_dir, exist_ok=True)  # Create the directory if it doesn't exist

    # Generate a list of .gjf files from the directory. Note that pseudopotentials write _atom## to the filename, so we filter these out. 
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x]
    num_files = len(filenames)

    start = time.time()

    #Error handling 

    #remove any leadng or trailing whitespace from the calc line
    calc_line = calc_line.strip()

    #check is calc line begns wth a !
    if not calc_line.startswith('!'):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The input line in {filename} does not start with a !')
        calc_line = '! ' + calc_line
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
            data = opf.read()

        # Get the geometry and write to a list called geometry
        geometry = [] # Get geometry of atoms from file
        GEOM = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)', data)

        # Loop through each line in the geom block, split each line, then append to the geometry list
        for i in GEOM[-1].split('\n')[2:-3]:
            i = i.split()
            if len(i) > 0:
                geometry.append(i)

        # Write the geometry to a formatted string
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

        if write_gjf_checked: 
            
            #directory to write .gjf files to
            gjf_dir = os.path.join(new_dir, 'gjfs')

            #make the .gjf file
            with open(os.path.join(gjf_dir, f'{filename[:-4]}.gjf'), 'w') as opf:
                opf.write('#opt\n\n')
                opf.write(f'{filename[:-4]}\n\n') #filename exclusing its extenstion 
                opf.write(f'{charge} {multiplicity}\n')
                opf.write(fs)
                opf.write('\n\n')
    
    #write a completion message to the GUI output window
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {num_files} ORCA .out files were converted to ORCA .inp files in {np.round(time.time() - start,2)} seconds.')





'''



def ORCA_out_to_ORCA_inp(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, calc_hess_checked, polarization_checked, write_xyz_checked):
    
    directory = r'G:\Hopkins_Laboratory\Protonation_Induced_Chirality_v2\Verapamil_Impurity_A_z1\SS\Final'

    # Important things for the input block
    mpp = 3900          # Memory per processor in MB. Must be an integer
    ncores = 8          # Number of cores to use in the calculation
    charge = 1          # What is the charge?
    multiplicity = 1    # What is the multiplicity?

    #Typical workflow
    calc_line = '! wB97X-D3 TightOpt Freq def2-TZVP def2/J RIJCOSX TightSCF defgrid3 CHELPG'    # DFT calcs, thermochemistry, and partial charges
    calc_line = '! DLPNO-CCSD(T) def2-TZVP def2-TZVP/C VeryTightSCF'                            # DLPNO-CCSD(T) calcs

    # optional flags

    ESP_Charges = False  # Do you want to calculate partial charges via ESP mapping using the ChelpG scheme?

    Polarization = False  # Do you want to calculate dipole and quadrupole moments?

    write_gjf = True    # Would you like .gjf files to be made in addition to the ORCA inputs? This option is useful for visualizing files in Gaussview. 
    calc_hess = False # Do you want to calculate the Hessian for the first opt step?
    
    ###########################################################################
    # No touchy past this point

    import os, re

    # Generate a directory to write new files to
    new_dir = os.path.join(directory, 'New_Inputs')
    os.makedirs(new_dir, exist_ok=True)  # Create the directory if it doesn't exist

    # Generate a list of ORCA output files from the directory, excluding xyz files that are trajectories from an optimization
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x]

    # Create ORCA files
    for filename in filenames:
        orca_filename = os.path.join(new_dir, f'{filename[:-4]}.inp')
        with open(orca_filename, 'w') as opf:
            # Input block
            if 'CHELPG' in calc_line and not ESP_Charges:
                ESP_Charges = True
            if 'CHELPG' not in calc_line and ESP_Charges:
                calc_line = f'{calc_line} CHELPG'

            if not calc_line.startswith('!'):
                calc_line = f'! {calc_line}'
            opf.write(f'{calc_line}\n')

            # Memory
            if not isinstance(mpp, int):
                mpp = round(mpp)
                print('The number of memory is not an integer. You cannot have fractional CPUs! Fixing it now')
            opf.write(f'%maxcore {mpp}\n\n')

            # Number of cores
            if not isinstance(ncores, int):
                ncores = round(ncores)
                print('The number of cores specified is not an integer. You cannot have fractional CPUs! Fixing it now')
            opf.write(f'%pal nprocs {ncores} \nend\n\n')

            # ESP Charges
            if ESP_Charges:
                
                grid = 0.1 # Default is 0.3
                rmax = 3.0 # Default is 2.8

                opf.write(f'%chelpg'+'\n')
                opf.write(f'grid {grid}\n')
                opf.write(f'rmax {rmax}\n')
                opf.write('end\n\n')

            # Open shell calcs
            if multiplicity != 1:
                opf.write(f'%scf HFTyp UHF\n')
                opf.write('end\n\n')

            # Calc Hessian during optimization
            if calc_hess:
                opf.write(f'%geom\n')
                opf.write('Calc_Hess true\n')
                opf.write('end\n\n')

            # Polarization
            if Polarization:
                opf.write(f'%elprop\n')
                opf.write('Dipole true\n')
                opf.write('Quadrupole true\n')
                opf.write('Polar 1\n')
                opf.write('end\n\n')

        with open(os.path.join(directory, filename), 'r') as opf:
            data = opf.read()

        # Get the geometry and write to a list called geometry
        geometry = [] # Get geometry of atoms from file
        GEOM = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)', data)

        # Loop through each line in the geom block, split each line, then append to the geometry list
        for i in GEOM[-1].split('\n')[2:-3]:
            i = i.split()
            if len(i) > 0:
                geometry.append(i)

        # Write the geometry to a formatted string
        fs = '\n'.join([f'{i[0]:<5s}  {i[1]:>15s}  {i[2]:>15s}  {i[3]:>15s}' for i in geometry])

        # At long last, we can write all the info to the .xyz file
        with open(orca_filename, 'a') as opf:
            opf.write(f'*xyz {charge} {multiplicity}\n')
            opf.write(fs)
            opf.write('\n*\n\n') #terminate geom block with an asterisk, followed by two blank lines

        if write_gjf:
            with open(os.path.join(new_dir, f'{filename[:-4]}.gjf'), 'w') as opf:
                opf.write('#opt\n\n')
                opf.write(f'{filename[:-4]}\n\n') #filename exclusing its extenstion 
                opf.write(f'{charge} {multiplicity}\n')
                opf.write(fs)
                opf.write('\n\n')

    print('Done')
'''