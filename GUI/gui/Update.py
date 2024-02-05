import git, os, sys, shutil, subprocess, stat, time

def Update_GUI_files(repo_url, root, ID_file, repo_SHA, delete_dir_function):
    '''Updates the .py files associated with the MobCal-MPI GUI. Inputs are the repo URL,  the directory where the GUI .py launcher is located, and a .txt file containing the SHA value of the user's local clone of the GUI.'''

    # Define the repository URL
    temp_dir = os.path.join(root, 'temp')

    # Remove the temporary directory if it exists
    if os.path.isdir(temp_dir):
        delete_dir_function(temp_dir)
    
    # Clone the GitHub repository to the temporary directory
    try: 
        repo = git.Repo.clone_from(repo_url, temp_dir, branch='master')
    except Exception as e: 
        print(f'Exception: {e}')
        print('Unable to access github to check for updates, likely due to lack of an internet connection...')
        answer = input('Do you wish to proceed using the current version (y/n)?')
        if answer == 'y':
            return
        else:
            raise KeyboardInterrupt('User has opted not to open the GUI without checking for updates. The GUI launcher will now be closed.')
    
    top_dir = os.path.dirname(root)

    # A handy dictionary to hold the paths of the files to be updated for subsequent looping. Syntax is as follows- Path to local file : Path to cloned GitHub file

    update_files = {
        str(os.path.join(top_dir, 'Sample_Files')): str(os.path.join(temp_dir, 'Sample_Files')), #Sample files to accompany to GUI
        str(os.path.join(top_dir, 'Documentation.docx')): str(os.path.join(temp_dir, 'Documentation.docx')), #GUI documentation
        str(os.path.join(root, 'Launcher.py')): str(os.path.join(temp_dir, 'GUI', 'Launcher.py')), #GUI launcher
        str(os.path.join(root, 'gui', 'Update.py')): str(os.path.join(temp_dir, 'GUI', 'gui', 'Update.py')), #Update function
        str(os.path.join(root, 'gui', 'ORCA_Analysis_GUI.py')): str(os.path.join(temp_dir, 'GUI', 'gui', 'ORCA_Analysis_GUI.py')), #GUI layout file

        str(os.path.join(root, 'Python', 'atom_mass.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'atom_mass.py')), #py file with Atom masses 
        str(os.path.join(root, 'Python', 'constants_and_conversions.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'constants_and_conversions.py')), #py file with fundamental contants & functions for unit conversions 

        str(os.path.join(root, 'Python', 'Input_Output_operations', 'xyz_file_splitter.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'xyz_file_splitter.py')), # T1 - xyz file splitter
        str(os.path.join(root, 'Python', 'Input_Output_operations', 'Gaussian_gjf_to_ORCA_input.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'Gaussian_gjf_to_ORCA_input.py')), # T2 - gjf to ORCA .inp
        str(os.path.join(root, 'Python', 'Input_Output_operations', 'cosine_sim.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'cosine_sim.py')), #T3 - cosine similarity sorting
        str(os.path.join(root, 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_inp.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_inp.py')), #T4- ORCA .out to ORCA .inp
        str(os.path.join(root, 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_TDDFT_VG.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Input_Output_operations', 'ORCA_out_to_ORCA_TDDFT_VG.py')), # T5 - ORCA .out to ORCA .inp to VGFC simluations (TD-DFT)
        
        str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'ORCA_opt_plt.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'ORCA_opt_plt.py')), # T6 - plotting the ORCA optimiation routine's progress 
        str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'ORCA_Thermochem_Calculator.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'ORCA_Thermochem_Calculator.py')), # T7 - Calculates thermochem from ORCA .out files 
        str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'ORCA_CCSDT.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'ORCA_CCSDT.py')), # T8 - Extracts CCSDT energies from ORCA .out files 
        str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'extract_IR.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'extract_IR.py')), # T9 - Exacts and plots IR spectra from ORCA .out files
        str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_root_files.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_root_files.py')), # T10 - Extract & plot UV spectra from .spectrum.root files
        str(os.path.join(root, 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_files.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'ORCA_out_analyses', 'extract_ESD_spectrum_files.py')), # T10 - Extract & plot UV spectra from .spectrum files
    
    
        str(os.path.join(root, 'Python', 'Special_Analyses', 'BW_CCS_Analyzer.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Special_Analyses', 'BW_CCS_Analyzer.py')), # T11 - Boltzmann-weighted CCS calculator
        str(os.path.join(root, 'Python', 'Special_Analyses', 'LED_Analyzer.py')): str(os.path.join(temp_dir, 'GUI', 'Python', 'Special_Analyses', 'LED_Analyzer.py')), # T12 - ORCA LED analysis tool
    }
    
    #update process for Windows users
    if os.getenv('APPDATA') is not None:
        
        #create a python script to update the relevant files
        update_script_path = os.path.join(temp_dir, 'update.py')
        
        with open(update_script_path, 'w') as opf:
            opf.write('import os, sys\n')
            opf.write('import shutil\n')
            
            # Write the logic to update each file
            for local_path, github_path in update_files.items():
                opf.write(f'if os.path.isfile(r"{local_path}"): os.remove(r"{local_path}")\n')
                opf.write(f'elif os.path.isdir(r"{local_path}"): shutil.rmtree(r"{local_path}")\n')
                opf.write(f'shutil.move(r"{github_path}", r"{os.path.dirname(local_path)}")\n')
            opf.write('sys.exit(0)')
            
            #ensure file updates; without this, users can encounter issues if the .py file is created on cloud services like OneDrive
            opf.flush()
            os.fsync(opf.fileno())
        
        time.sleep(5) #a short delay so that Cloud-based services have enough time to update.

        # Execute the update script and exit
        try:
            subprocess.Popen(['python', update_script_path], shell=True)

        except subprocess.SubprocessError as e:
            print(f'An subprocess error occurred while executing {update_script_path}: {e}')

        except Exception as e:
            print(f'An unexpected error occurred when trying execute {update_script_path}: {e}')

        with open(ID_file,'w') as opf:
            opf.write(repo_SHA)

        return
    
    #Update process for Linux users
    else:
        for local_path, github_path in update_files.items():
            # Similar to Windows, but using direct Python commands instead of writing to a script
            if os.path.isfile(local_path):
                os.remove(local_path)
            shutil.move(github_path, os.path.dirname(local_path))

        #Update the ID file
        with open(ID_file, 'w') as opf:
            opf.write(repo_SHA)

        return
