#Import modules needed for preparation of GUI launch
import importlib, platform, os, sys, shutil, subprocess
from pathlib import Path

# Before the GUI launches, check that the user has the required packages to run the MobCal-MPI GUI
#The most troublesome package is Git, which also requires GitHub desktop to be on the user's machine. First, we check if it is installed.

def find_github_desktop():
    '''A function to dynamically locate GitHub Desktop'''
    os_type = platform.system()
    paths_to_check = []

    #Check process for Windows users
    if os_type == 'Windows':
        
        # First, try to find the path in the registry
        try:
            #Irrespective of install location of Git desktop, there should always be a registry key here for windows users
            import winreg
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, r'Software\Microsoft\Windows\CurrentVersion\Uninstall\GitHubDesktop') as key:
                path, _ = winreg.QueryValueEx(key, 'InstallLocation')
                git_exe = Path(path) / 'git.exe'
                if git_exe.exists():
                    return str(git_exe.parent)
        except FileNotFoundError:
            pass
        except OSError:
            pass

        #If unsucessful, try some other default locations.  
        paths_to_check.extend([
            Path(os.environ.get('LOCALAPPDATA', '')) / 'GitHubDesktop',
            Path(os.environ.get('PROGRAMFILES', '')) / 'GitHub Desktop',
            Path(os.environ.get('PROGRAMFILES(X86)', '')) / 'GitHub Desktop',
            Path(os.environ.get('USERPROFILE', '')) / 'AppData' / 'Local' / 'GitHubDesktop'
            # Users can add more Windows-specific paths here if they installed GitHub Desktop to a non-default location.
        ])

    # Check process for Mac users (I don't have a MAc so I have not been able to check if this works; I'm going off of stack exchange here)
    elif os_type == 'Darwin':
        paths_to_check.extend([
            Path('/Applications/GitHub Desktop.app'),
            Path.home() / 'Applications' / 'GitHub Desktop.app'
        ])
        # Users can add more Mac-specific paths here if they installed GitHub Desktop to a non-default location.

    # Check process for Linux users (I don't have have a Linux machine so I have not been able to check if this works; I'm going off of stack exchange here)
    elif os_type == 'Linux':
        paths_to_check.extend([
            Path('/usr/bin/github-desktop'),
            Path('/usr/local/bin/github-desktop')
        ])
        # Users can add more Mac-specific paths here if they installed GitHub Desktop to a non-default location.
    
    for path in paths_to_check:
        if (os_type == 'Windows' and path.is_dir()) or (os_type in ['Darwin', 'Linux'] and path.exists()):
            return str(path)
    
    return None

def find_git_executable():
    os_type = platform.system()
    
    if os_type == 'Windows':
        # Try to find the Git path in the registry
        try:
            import winreg
            # Check both HKEY_LOCAL_MACHINE and HKEY_CURRENT_USER
            for hkey in [winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER]:
                with winreg.OpenKey(hkey, r'Software\GitForWindows') as key:
                    path, _ = winreg.QueryValueEx(key, 'InstallPath')
                    git_exe = Path(path) / 'bin' / 'git.exe'
                    if git_exe.exists():
                        return str(git_exe)
        except FileNotFoundError:
            pass
        except OSError:
            pass

    # For macOS and Linux, as well as a failsafe for Windows
    git_path = shutil.which('git')
    if git_path is not None:
        return git_path

    # Common paths to check for Git on Unix-like systems
    unix_paths = [
        '/usr/local/bin/git',  # Common for macOS and some Linux distros
        '/opt/local/bin/git',  # Common for installations via MacPorts
        '/usr/bin/git',        # Common for Linux distros
        '/bin/git',            # Less common, but worth checking
        # more can be added as needed
    ]

    for path in unix_paths:
        if Path(path).is_file():
            return path

    # If Git executable not found, idk you probably didn't install it
    return None
    
def add_to_path(new_path):
    '''Takes a path as input and adds it to the system's PATH environment variable if it isn't already there'''
    os_type = platform.system()

    # Windows uses semicolons (;) as path separator
    if os_type == 'Windows':
        path_separator = ';'
    else:
        # Both macOS (Darwin) and Linux use colons (:) as path separator
        path_separator = ':'

    # Check if new_path is already in PATH
    if new_path not in os.environ['PATH'].split(path_separator):
        os.environ['PATH'] = new_path + path_separator + os.environ['PATH']

def check_git():
    '''Checks if GitHub Desktop and Git are installed on the Users PC and adds it to PATH, and if it isn't, promts them to instal it before continuing'''
    
    #Check for GitHub Desktop
    git_desktop_path = find_github_desktop()

    if not git_desktop_path:
        print('GitHub Desktop is not installed. Please install from the following URL to the default directory:')
        print('https://desktop.github.com/')
        sys.exit(0)
    
    else:
        #If GitHub desktop is found, add it to the system's PATH
        add_to_path(git_desktop_path)

    #Check for Git
    git_path = find_git_executable()

    if not git_path:
        print('Git is not installed. Please install from the following URL to the default directory:')
        print('https://git-scm.com/')
        sys.exit(0)
    
    else:
        #If GitHub is found, add it to the system's PATH
        add_to_path(git_path)

#Now with those functions set up, we can check for the required python modules 
def check_python_packages(required_packages):
    for package in required_packages:
        try:
            importlib.import_module(package)
        except ModuleNotFoundError:
            offer_package_install(package)

def offer_package_install(package):
    print(f'{package} cannot be found. Would you like to install it (y/n)?')
    if input().lower() == 'y':
        print(f'Installing {package}...')
        #git imports as gitpython despite is being called git (github why????????????????????)
        package_name = 'gitpython' if package == 'git' else package
        subprocess.run([sys.executable, '-m', 'pip', 'install', package_name])
    else:
        sys.exit('Exiting: Required package not installed.')

if __name__ == '__main__':
    #Check that github desktop is installed
    check_git()

    #check that all required python modules are installed
    required_packages = ['numpy', 'scipy', 'matplotlib', 'git', 'lxml', 'PyQt6', 'pandas', 'openpyxl', 'csv', 'natsort']
    check_python_packages(required_packages)

#Now that all the required packages are installed, we can import the modules/functions used by the GUI

def check_dependencies():
    #Setup a dictionary showing of the module dependencies of required modules and their functions, organized by functionality
    dependencies = {
        'Input_Output_operations': [
            'Python.Input_Output_operations.xyz_file_splitter.xyz_file_splitter',
            'Python.Input_Output_operations.Generate_ORCA_inp.Generate_ORCA_inp',
            'Python.Input_Output_operations.cosine_sim.cosine_sim',
            'Python.Input_Output_operations.ORCA_out_to_ORCA_inp.ORCA_out_to_ORCA_inp',
            'Python.Input_Output_operations.ORCA_out_to_ORCA_TDDFT_VG.ORCA_out_to_ORCA_TDDFT_VG',
        ],
        'ORCA_out_analyses': [
            'Python.ORCA_out_analyses.ORCA_opt_plt.ORCA_opt_plt',
            'Python.ORCA_out_analyses.ORCA_Thermochem_Calculator.ORCA_Thermochem_Calculator',
            'Python.ORCA_out_analyses.ORCA_CoupledCluster.extract_ORCA_coupled_cluster',
            'Python.ORCA_out_analyses.extract_IR.extract_IR_spectra',
            'Python.ORCA_out_analyses.extract_ESD_spectrum_files.extract_ESD_spectrum_files',
            'Python.ORCA_out_analyses.extract_ESD_spectrum_root_files.extract_ESD_spectrum_root_files',
        ],
        'Special_Analyses': [
            'Python.Special_Analyses.LED_Analyzer.LED_Analysis',
            'Python.Special_Analyses.BW_CCS_Analyzer.BW_CCS_Analysis',
            #'Python.Special_Analyses.NEB_Analyzer.NEB_Analyzer',
        ],
        'General_utilities': [
            'Python.atom_mass.atom_masses',
            'Python.constants_and_conversions.c',
            'Python.constants_and_conversions.convert_energy',
        ]
    }

    missing_dependencies = []

    #Attempt to import each dependency
    for category, funcs in dependencies.items():
        for func in funcs:
            try:
                module_name, function_name = func.rsplit('.', 1) #split at first occurance of the period
                module = __import__(module_name, fromlist=[function_name])
                getattr(module, function_name) #test import
            
            except ImportError as e:
                missing_dependencies.append(func)
            
            except AttributeError as e:
                missing_dependencies.append(func)

    if missing_dependencies:
        print('The following Python modules are missing. Please reclone the repo from Github and try again, ensuring the missing Python modules are present:\n', ",\n".join(missing_dependencies))
        return False
    
    else:
        return True

#import everything from the main GUI .py file, including the check_for_update_and_prompt function called below
from gui.ORCA_Analysis_GUI import *

if __name__ == '__main__':
    #Check for all dependencies
    dependencies_check = check_dependencies()

    #If dependencies are present, proceed to launch the GUI
    if dependencies_check:
        app = QApplication(sys.argv)
        main_win = ORCAAnalysisSuite()
        main_win.show() #show the main window
        #main_win.check_for_update_and_prompt()  #Check for updates after showing the window    
        sys.exit(app.exec())