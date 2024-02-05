import os, re, time
import numpy as np
from datetime import datetime
from natsort import natsorted
import shutil

def ORCA_out_to_ORCA_TDDFT_VG(directory, mpp, ncores, charge, multiplicity, calc_line, states, write_gjf_checked):

    # Generate a list of .gjf files from the directory. Note that pseudopotentials write _atom## to the filename, so we filter these out. Also sort in natural order for pairwise comparison 
    out_filenames = natsorted([x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x])
    hess_filenames = natsorted([x for x in os.listdir(directory) if x.lower().endswith('.hess') and '_atom' not in x])
    
    if not out_filenames:
        raise ValueError(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .out files could be found in {directory}. Please ensure that you have entered the correct directory and/or .out files are present where you think they are!')

    if not hess_filenames:
        raise ValueError(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .hess files could be found in {directory}. Please ensure that you have entered the correct directory and/or .hess files are present where you think they are!')

    if not len(out_filenames) == len(hess_filenames):
        raise ValueError(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of .out files is different than the number of .hess files in {directory}. Please ensure that all .hess files are present')

    start = time.time() 
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {len(out_filenames)} ORCA .out files to TD-DFT .inp files.')

    # Check if output directory exists, create if not
    new_dir = os.path.join(directory, f'TDDFT')
    i = 1
    while os.path.exists(new_dir):
        try:
            new_dir = f'{os.path.join(directory, f"TDDFT_{i}")}' # Alex
            os.makedirs(new_dir)
            break  # Break the loop if the directory is successfully created
        except FileExistsError:
            i += 1

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    #remove any leadng or trailing whitespace from the calc line
    calc_line = calc_line.strip()

    #check is calc line begns wth an exclamation point
    if not calc_line.startswith('!'):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The calculation line does not start with an exclamation point (!)')
        calc_line = f'! {calc_line}'
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Fixing it now; the new calc line is:\n{calc_line}')
       
    #check that ESD(ABS) is in the calc line, and that other variations of ESD(ABS) are not
    if 'ESD(ABS)' not in calc_line.upper():
        # Define a regular expression pattern to find any variation of ESD(), where the contents of the parentesis can be anything
        
        pattern = re.compile(r'\bESD\([^)]*\)\b')
        bad_ESD = re.findall(pattern, calc_line)
        
        if not bad_ESD: # if ESD is given in the calc_line without any parenthesis
            pattern = re.compile(r'\bESD\b')
            bad_ESD = re.findall(pattern, calc_line)

        #replace the erronous ESD() with ESD(ABS)
        calc_line = re.sub(pattern, 'ESD(ABS)', calc_line)

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {"".join(bad_ESD)} was found in the calc line when the only acceptable input is ESD(ABS) for this implementation of generating inputs via the VG approximation. Please consult the ORCA manual for other variations, and generate these files manually.')
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} In case this was unintentional, the calc line has been modified. The new calc line is:\n{calc_line}')

    #check that the memory allocated per core is an integer
    if not isinstance(mpp, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of memory is not an integer. Fixing it now')
        mpp = round(mpp)
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The mpp being written to each .inp is {mpp}.')

    #check that the number of cores requested for the job is an integer
    if not isinstance(ncores, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores specified is not an integer. You cannot have fractional CPUs! Fixing it now')
        ncores = round(ncores)   
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores being written to each .inp is {ncores}.')

    # Create ORCA .inp files
    for out, hess in zip(out_filenames, hess_filenames):

        out_base = out.split('.out')[0]
        hess_base = hess.split('.hess')[0]

        #if the filenames do not match, raise a value error
        if not out_base == hess_base: 
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The .out file {out} has a different name than the assoicated .hess file from the same list {hess}. Please verify that each .out file has a corresponding .hess file with the same basename.')
            raise ValueError()

        #read and extract geometry from the .out files
        with open(os.path.join(directory, out), 'r') as opf:
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
        geom_fs = '\n'.join([f'{i[0]:<5s}  {i[1]:>15s}  {i[2]:>15s}  {i[3]:>15s}' for i in geometry])
            
        #copy the .hess file to the new directory
        try:
            shutil.copy2(os.path.join(directory, hess), os.path.join(new_dir, hess))
        except FileNotFoundError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {hess} could not be found for some weird reason.')
            raise FileNotFoundError()

        except PermissionError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {hess} could not be copied due to a permission issue.')
            raise PermissionError()

        #write TD-DFT .inp file to the new directory and populate its contents
        orca_filename = os.path.join(new_dir, f'{out[:-4]}.inp')

        with open(orca_filename, 'w') as opf:
            # Input block
            opf.write(f'{calc_line}\n')

            # Memory
            opf.write(f'%maxcore {mpp}\n\n')

            # Number of cores
            opf.write(f'%pal nprocs {ncores} \nend\n\n')

            # Open shell calcs
            if multiplicity != 1:
                opf.write(f'%scf HFTyp UHF\n')
                opf.write('end\n\n')

            #TD-DFT block
            opf.write('%tddft\n')
            opf.write(f'nroots={states}\n')
            opf.write('maxdim 5\n') #Davidson expansion space = MaxDim * nroots. Use MaxDim 5-10 for favorable convergence. Note that the larger MaxDim is, the more disk space is required
            opf.write('tda false\n')
            opf.write('end\n\n')

            #Excited state dynamics block - default spectrum printing to nm, but this doesn't matter. It can be converted to eV, wavenumber, etc. during extraction alongside assigned custom line broadening
            opf.write('%esd\n')
            opf.write(f'GSHESSIAN "{hess}"\n')
            opf.write('hessflag VG\n')
            opf.write(f'STATES {",".join(str(state) for state in range(1, states+1))}\n')
            opf.write('unit nm\n')
            opf.write('LINES gaussian\n')
            opf.write('USEJ FALSE\n')
            opf.write('end\n\n')
            
            #write the geometry to the .inp
            opf.write(f'*xyz {charge} {multiplicity}\n{geom_fs}\n*\n\n')

        #write .gjf files to the directory with the .inp if the user has requested it
        if write_gjf_checked: 
            
            #make the .gjf file and write the geom info to it
            with open(os.path.join(new_dir, f'{out[:-4]}.gjf'), 'w') as opf:
                opf.write('#opt\n\n') #calc line irrelevant for the .gjf file since this is only for visualization
                opf.write(f'{out[:-4]}\n\n') #filename exclusing its extenstion 
                opf.write(f'{charge} {multiplicity}\n')
                opf.write(geom_fs)
                opf.write('\n\n')
        
    #write a completion message to the GUI output window
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(out_filenames)} ORCA .out files were converted to ORCA .inp files in {np.round(time.time() - start,2)} seconds.')

#external testing
if __name__ == '__main__':
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Manual\Appendix_A\DFT'
    mpp = 3400
    ncores = 8
    charge = 1
    multiplicity = 1
    calc_line = r'! wB97X-D3 Def2-TZVPP ESD(ABS)' #dummy test line - this doesn't actualyl work if given to ORCA
    states = 15
    write_gjf_checked = True

    #run the code
    ORCA_out_to_ORCA_TDDFT_VG(directory, mpp, ncores, charge, multiplicity, calc_line, states, write_gjf_checked)