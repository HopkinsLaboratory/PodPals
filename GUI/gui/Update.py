import git, os, sys, shutil, subprocess, stat, time
from datetime import datetime
from PyQt6.QtWidgets import QApplication

def Update_GUI_files(repo_url, root, ID_file, repo_SHA, delete_dir_function, ensure_update):
    '''Updates the .py files associated with the MobCal-MPI GUI. Inputs are the repo URL,  the directory where the GUI .py launcher is located, and a .txt file containing the SHA value of the user's local clone of the GUI.'''

    #Define the repository URL
    temp_dir = os.path.join(root, 'temp')
    top_dir = os.path.dirname(root)

    #only clone the repo once when updating the Update.py file
    if ensure_update:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Cloning {repo_url} ...')
        QApplication.processEvents()

        #Remove the temporary directory if it exists
        if os.path.isdir(temp_dir):
            delete_dir_function(temp_dir)
        
        #Clone the GitHub repository to the temporary directory
        try: 
            repo = git.Repo.clone_from(repo_url, temp_dir, branch='master')
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Cloning complete.')
            QApplication.processEvents()
        
        except Exception as e: 
            print(f'An exception occured when attempting to clone the ORCA Analysis GUI repo: {e}')
            print('This block is most commonly entered due to lack of an internet connection, although the error message above may indicate otherwise... Use your best judgement :)')
            print('Please feel free to use the local version of the GUI and update once your internet connection is retored.')
            QApplication.processEvents()
            return

    #A handy dictionary to hold the paths of the files to be updated for subsequent looping. Syntax is as follows- Path to local file : Path to cloned GitHub file
    if ensure_update:
        print(f"""{datetime.now().strftime("[ %H:%M:%S ]")} Ensuring that the local version of {os.path.join(root, 'gui', 'Update.py')} is up to date...""")
        QApplication.processEvents()
        
        #create a location to write a python script to update Update.py to the most recent version
        update_script_path = os.path.join(temp_dir, 'initial_update.py')
        
        #ensure that the Update.py script is replaced with the most recent version before updating any other files, as dependencies can change over time
        update_files = {str(os.path.join(root, 'gui', 'Update.py')): str(os.path.join(temp_dir, 'GUI', 'gui', 'Update.py'))} #Update function
    
    else:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Updating all .py modules to their most recent version...')
        QApplication.processEvents()

        #create a location to write a python script to update Update.py to the most recent version
        update_script_path = os.path.join(temp_dir, 'Allfiles_update.py')

        update_files = {
            str(os.path.join(top_dir, 'Sample_Files')): str(os.path.join(temp_dir, 'Sample_Files')), #Sample files to accompany to GUI
            str(os.path.join(top_dir, 'README.md')): str(os.path.join(temp_dir, 'README.md')), #GUI documentation
            str(os.path.join(root, 'Launcher.py')): str(os.path.join(temp_dir, 'GUI', 'Launcher.py')), #GUI launcher
            #str(os.path.join(root, 'gui', 'Update.py')): str(os.path.join(temp_dir, 'GUI', 'gui', 'Update.py')), #Update function
            str(os.path.join(root, 'gui', 'ORCA_Analysis_GUI.py')): str(os.path.join(temp_dir, 'GUI', 'gui', 'ORCA_Analysis_GUI.py')), #GUI layout file
            str(os.path.join(root, 'gui', 'icon.png')): str(os.path.join(temp_dir, 'GUI', 'gui', 'icon.png')), #The icon for the GUI window

            str(os.path.join(root, 'Python', 'atom_mass.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'atom_mass.py')), #py file with Atom masses 
            str(os.path.join(root, 'Python', 'constants_and_conversions.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'constants_and_conversions.py')), #py file with fundamental contants & functions for unit conversions 

            str(os.path.join(root, 'Python', 'Input_Output_operations', 'xyz_file_splitter.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'xyz_file_splitter.py')), #T1 - xyz file splitter
            str(os.path.join(root, 'Python', 'Input_Output_operations', 'Generate_ORCA_inp.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'Generate_ORCA_inp.py')), #T2 - Generate ORCA .inp
            str(os.path.join(root, 'Python', 'Input_Output_operations', 'cosine_sim.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'cosine_sim.py')), #T3 - cosine similarity sorting
            str(os.path.join(root, 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_inp.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_inp.py')), #T4- ORCA .out to ORCA .inp
            str(os.path.join(root, 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_TDDFT_VG.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_TDDFT_VG.py')), #T5 - ORCA .out to ORCA .inp to VGFC simluations (TD-DFT)
            
            str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'ORCA_opt_plt.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'ORCA_opt_plt.py')), #T6 - plotting the ORCA optimiation routine's progress 
            str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'ORCA_Thermochem_Calculator.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'ORCA_Thermochem_Calculator.py')), #T7 - Calculates thermochem from ORCA .out files 
            str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'ORCA_CoupledCluster.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'ORCA_CoupledCluster.py')), #T8 - Extracts coupled cluster energies from ORCA .out files 
            str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'extract_IR.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'extract_IR.py')), #T9 - Exacts and plots IR spectra from ORCA .out files
            str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_root_files.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_root_files.py')), #T10 - Extract & plot UV spectra from .spectrum.root files
            str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_files.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_files.py')), #T10 - Extract & plot UV spectra from .spectrum files
        
            str(os.path.join(root, 'Python', 'Special_Analyses', 'BW_CCS_Analyzer.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Special_Analyses', 'BW_CCS_Analyzer.py')), #T11 - Boltzmann-weighted CCS calculator
            str(os.path.join(root, 'Python', 'Special_Analyses', 'LED_Analyzer.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Special_Analyses', 'LED_Analyzer.py')), #T12 - ORCA LED analysis tool
    }
    
    #update process for Windows users
    if os.getenv('APPDATA') is not None:
               
        with open(update_script_path, 'w') as opf:
            opf.write('import os, sys\n')
            opf.write('import shutil\n')
            
            #Write the logic to update each file
            for local_path, github_path in update_files.items():
                opf.write(f'if os.path.isfile(r"{local_path}"): os.remove(r"{local_path}")\n')
                opf.write(f'elif os.path.isdir(r"{local_path}"): shutil.rmtree(r"{local_path}")\n')
                opf.write(f'shutil.move(r"{github_path}", r"{os.path.dirname(local_path)}")\n')
            opf.write('sys.exit(0)')
            
            #ensure file updates; without this, users can encounter issues if the .py file is created on cloud services like OneDrive
            opf.flush()
            os.fsync(opf.fileno())
        
        time.sleep(5) #a short delay so that Cloud-based services have enough time to update.

        #Execute the update script and exit
        try:
            subprocess.Popen(['python', update_script_path], shell=True)

        except subprocess.SubprocessError as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} A subprocess error occurred while updating files from the script {update_script_path}: {e}')
            return

        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An unexpected error occurred when trying execute {update_script_path}: {e}')
            print('Please report this error and the sequence of events that caused it to appear to the issues section of the GitHub repo.')
            return

        with open(ID_file,'w') as opf:
            opf.write(repo_SHA)

        return
    
    #Update process for Mac/Linux users
    else:
        for local_path, github_path in update_files.items():
            #Similar to Windows, but using direct Python commands instead of writing to a script
            if os.path.isfile(local_path):
                os.remove(local_path)
            shutil.move(github_path, os.path.dirname(local_path))

        #Update the ID file
        with open(ID_file, 'w') as opf:
            opf.write(repo_SHA)

        return
