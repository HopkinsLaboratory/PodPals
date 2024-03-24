import os, re, time
import numpy as np
from datetime import datetime
from natsort import natsorted
import shutil
from PyQt6.QtWidgets import QApplication

def ORCA_out_to_ORCA_TDDFT_VG(directory, mpp, ncores, charge, multiplicity, calc_line, states, write_gjf_checked):

    #Generate a list of .gjf files from the directory. Note that pseudopotentials write _atom##to the filename, so we filter these out. Also sort in natural order for pairwise comparison 
    out_filenames = natsorted([x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x])
    hess_filenames = natsorted([x for x in os.listdir(directory) if x.lower().endswith('.hess') and '_atom' not in x])
    
    if not out_filenames:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .out files could be found in {directory}. Please ensure that you have entered the correct directory and/or .out files are present where you think they are!')
        return

    if not hess_filenames:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .hess files could be found in {directory}. Please ensure that you have entered the correct directory and/or .hess files are present where you think they are!')
        return

    if not len(out_filenames) == len(hess_filenames):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of .out files is different than the number of .hess files in {directory}. Please ensure that all .hess files are present')
        return

    start = time.time() 
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting conversion of {len(out_filenames)} ORCA .out files to VG-FC ORCA .inp files.')
    QApplication.processEvents()

    #Check if output directory exists, create if not
    new_dir = os.path.join(directory, f'TDDFT_VGFC')
    os.makedirs(new_dir, exist_ok=True)  
       
    if write_gjf_checked:
        gjf_dir = os.path.join(new_dir, 'gjfs')
        print(f"""{datetime.now().strftime("[ %H:%M:%S ]")} .gjf files will be written to {os.path.join(os.path.basename(new_dir), 'gjfs')}""")
        QApplication.processEvents()

        os.makedirs(gjf_dir, exist_ok=True) 

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    #remove any leadng or trailing whitespace from the calc line
    calc_line = calc_line.strip()

    #check is calc line begns wth an exclamation point
    if not calc_line.startswith('!'):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The calculation line does not start with an exclamation point (!)')
        calc_line = f'! {calc_line}'
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Fixing it now; the new calc line is:\n{calc_line}')
        QApplication.processEvents()
       
    #check that ESD(ABS) is in the calc line, and that other variations of ESD(ABS) are not
    if 'ESD(ABS)' not in calc_line.upper():
        #Define a regular expression pattern to find any variation of ESD(), where the contents of the parentesis can be anything
        
        pattern = re.compile(r'\bESD\([^)]*\)\b')
        bad_ESD = re.findall(pattern, calc_line)
        
        if not bad_ESD: #if ESD is given in the calc_line without any parenthesis
            pattern = re.compile(r'\bESD\b')
            bad_ESD = re.findall(pattern, calc_line)

        #replace the erronous ESD() with ESD(ABS)
        calc_line = re.sub(pattern, 'ESD(ABS)', calc_line)

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {"".join(bad_ESD)} was found in the calc line when the only acceptable input is ESD(ABS) for this implementation of generating inputs via the VG approximation. Please consult the ORCA manual for other variations, and generate these files manually.')
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} In case this was unintentional, the calc line has been modified. The new calc line is:\n{calc_line}')
        QApplication.processEvents()

    #check that the memory allocated per core is an integer
    if not isinstance(mpp, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of memory is not an integer. Fixing it now')
        mpp = round(mpp)
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The mpp being written to each .inp is {mpp}.')
        QApplication.processEvents()

    #check that the number of cores requested for the job is an integer
    if not isinstance(ncores, int):
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores specified is not an integer. You cannot have fractional CPUs! Fixing it now')
        ncores = round(ncores)   
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The number of cores being written to each .inp is {ncores}.')
        QApplication.processEvents()

    #Create ORCA .inp files
    for out, hess in zip(out_filenames, hess_filenames):

        out_base = out.split('.out')[0]
        hess_base = hess.split('.hess')[0]

        #if the filenames do not match, raise a value error
        if not out_base == hess_base: 
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The .out file {out} has a different name than the assoicated .hess file from the same list {hess}.\nAs the .hess file created after an ORCA freq calc finishes has the name basename as the .out file, it is likely that it got deleted somewhere.\n Processing of this file pairing will be skipped.')
            QApplication.processEvents()
            continue

        #read and extract geometry from the .gjf files
        try:
            #read and extract geometry from the .gjf files
            with open(os.path.join(directory, out), 'r') as opf:
                data = opf.read()

        #If there is an error opening the .gjf file, infomr that user that it will be skipped.
        except IOError as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error was encountered when trying to open {filename}: {e}.\n Processing of this file will be skipped.')
            data = None
            continue

        #Get the geometry and write to a list called geometry
        geometry = [] #Initialize the geometry list
        XYZ_data = re.findall(r'CARTESIAN COORDINATES \(ANGSTROEM\)([\s\S]*?)CARTESIAN COORDINATES \(A.U.\)', data)

        #Loop through each line in the FINAL geom block, split each line, then append to the geometry list
        for line in XYZ_data[-1].split('\n')[2:-3]:
            split_line = line.split()
            #the only entry with 4 splits, where instance 1 is an atom symbol (sometimes an atom number!), and a period in instances 1-3 will be the xyz coordiante lines. Write these to the geom list
            if len(split_line) == 4 and (split_line[0].isalpha() or split_line[0].isdigit()) and all('.' in x for x in split_line[1:3]):
                geometry.append(split_line)

        if not geometry:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No atomic coordinate data was found in {filename}. Processing of this file will be skipped.')
            QApplication.processEvents()
            continue
        
        #If geometry is extracted, start to write the ORCA .inp file
        else:

            #copy the .hess file to the new directory
            try:
                shutil.copy2(os.path.join(directory, hess), os.path.join(new_dir, hess))

            except FileNotFoundError:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {hess} could not be found for some weird reason. The {out} and {hess} pairing will be skipped.')
                QApplication.processEvents()
                continue

            except PermissionError:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {hess} could not be copied due to a permission issue. The {out} and {hess} pairing will be skipped.')
                QApplication.processEvents()
                continue

            #write TD-DFT .inp file to the new directory and populate its contents
            orca_filename = os.path.join(new_dir, f'{out[:-4]}.inp')

            #Write the geometry to a formatted string
            fs = '\n'.join([f'{i[0]:<5s}  {i[1]:>15s}  {i[2]:>15s}  {i[3]:>15s}' for i in geometry])
            
            with open(orca_filename, 'w') as opf:
                
                #Input block
                opf.write(f'{calc_line}\n')

                #Memory
                opf.write(f'%maxcore {mpp}\n\n')

                #Number of cores
                opf.write(f'%pal nprocs {ncores} \nend\n\n')

                #Open shell calcs
                if multiplicity != 1:
                    opf.write(f'%scf HFTyp UHF\n')
                    opf.write('end\n\n')

                #TD-DFT block
                opf.write('%tddft\n')
                opf.write(f'nroots={states}\n')
                #opf.write('maxdim 5\n') #Davidson expansion space = MaxDim * nroots. Use MaxDim 5-10 for favorable convergence. Note that the larger MaxDim is, the more disk space is required
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
                opf.write(f'*xyz {charge} {multiplicity}\n{fs}\n*\n\n')

            #write .gjf files to the directory with the .inp if the user has requested it
            if write_gjf_checked: 
                with open(os.path.join(gjf_dir, f'{out[:-4]}.gjf'), 'w') as opf_gjf:
                    opf_gjf.write('#opt\n\n') #calc line irrelevant for the .gjf file since this is only for visualization
                    opf_gjf.write(f'{out[:-4]}\n\n') #filename exclusing its extenstion 
                    opf_gjf.write(f'{charge} {multiplicity}\n')
                    opf_gjf.write(fs)
                    opf_gjf.write('\n\n')
            
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(out_filenames)} ORCA .out files were converted to VG-FC ORCA .inp files in {np.round(time.time() - start,2)} seconds.')
    return

#external testing
if __name__ == '__main__':
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Manual\Appendix_A\DFT'
    mpp = 3400
    ncores = 8
    charge = 1
    multiplicity = 1
    calc_line = r'! wB97X-D3 Def2-TZVPP ESD(ABS)' #dummy test line 
    states = 15
    write_gjf_checked = True

    #run the code
    ORCA_out_to_ORCA_TDDFT_VG(directory, mpp, ncores, charge, multiplicity, calc_line, states, write_gjf_checked)
