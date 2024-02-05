from PyQt6.QtWidgets import QApplication, QHBoxLayout, QComboBox, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QMessageBox, QLabel, QGroupBox, QLineEdit, QPushButton, QFileDialog, QTextEdit, QCheckBox, QSpinBox, QSizePolicy, QDoubleSpinBox
from PyQt6.QtCore import QCoreApplication, Qt
from PyQt6.QtGui import QScreen

import os, time, sys, subprocess, traceback
from datetime import datetime
from io import StringIO

#Import the functions for the I/O tab
from Python.Input_Output_operations.xyz_file_splitter import xyz_file_splitter
from Python.Input_Output_operations.Gaussian_gjf_to_ORCA_input import Gaussian_gjf_to_ORCA_input
from Python.Input_Output_operations.cosine_sim import cosine_sim
from Python.Input_Output_operations.ORCA_out_to_ORCA_inp import ORCA_out_to_ORCA_inp
from Python.Input_Output_operations.ORCA_out_to_ORCA_TDDFT_VG import ORCA_out_to_ORCA_TDDFT_VG

#Import the functions for the output analysis Tab
from Python.ORCA_out_analyses.ORCA_opt_plt import ORCA_opt_plt
from Python.ORCA_out_analyses.ORCA_Thermochem_Calculator import ORCA_Thermochem_Calculator
from Python.ORCA_out_analyses.ORCA_CCSDT import ORCA_CCSDT
from Python.ORCA_out_analyses.extract_IR import extract_IR_spectra
from Python.ORCA_out_analyses.extract_ESD_spectrum_files import extract_ESD_spectrum_files
from Python.ORCA_out_analyses.extract_ESD_spectrum_root_files import extract_ESD_spectrum_root_files

#Import the functions for the Additional Analyses tab
from Python.Special_Analyses.LED_Analyzer import LED_Analysis
from Python.Special_Analyses.BW_CCS_Analyzer import BW_CCS_Analysis
#from Python.Special_Analyses.NEB_Analyzer import NEB_Analyzer

#Import the GitHub update function
from gui.Update import Update_GUI_files

class TextRedirect(StringIO):
    #Constructor (__init__ method) for the custom stream class
    def __init__(self, update_output=None, *args, **kwargs):
        
        #Call the constructor of the parent class (StringIO)
        super().__init__(*args, **kwargs)
        
        #Store a callback function for updating external components with the written text
        self.update_output = update_output

    #Override the write method of the parent class (StringIO)
    def write(self, text):
        
        #Call the write method of the parent class to perform the default writing
        super().write(text)
        
        #Invoke the stored callback function to notify external components with the written text
        self.update_output(text)

class ORCAAnalysisSuite(QMainWindow):
    def __init__(self):
        super().__init__()

        self.output_text_edit = QTextEdit()
        self.text_redirector = TextRedirect(update_output=self.update_output_text)

        self.initUI()

    def initUI(self):

        #Set up GUI main window, define tab layout, and generate output text box 
        self.setWindowTitle('ORCA Analysis Suite')

        #Dynamically adjust the size of the GUI window with the size of the user's screen
        #Get the primary screen's geometry
        screen = QApplication.primaryScreen().geometry()
        
        #Set GUI size to 30% of screen width x 40% of screen height, centered on the screen
        width = screen.width() * 0.3
        height = screen.height() * 0.4
        self.setGeometry(int((screen.width() - width) / 2), int((screen.height() - height) / 2), int(width), int(height))

        main_widget = QWidget(self)
        main_layout = QVBoxLayout(main_widget)

        main_tabs = QTabWidget()

        #Main Tabs - Input/Output operations and output analysis
        IO_tab = QTabWidget()
        output_analysis_tab = QTabWidget()
        extras_tab = QTabWidget()

        #Sub-tabs for I/O operations, then add to main IO tab

        self.xyz_file_splitter_tab = XYZFileSplitterTab(self.text_redirector)
        IO_tab.addTab(self.xyz_file_splitter_tab, 'XYZ File Splitter')

        sub_tab2 = GJFtoORCAInputTab(self.output_text_edit)
        IO_tab.addTab(sub_tab2, 'GJF to ORCA Input')

        sub_tab3 = CosineSimTab(self.output_text_edit)
        IO_tab.addTab(sub_tab3, 'Cosine Sim')

        sub_tab4 = ORCAOut_ORCAInputTab(self.output_text_edit)
        IO_tab.addTab(sub_tab4, 'ORCA out to ORCA inp')

        sub_tab5 = ORCAOut_ORCA_TDDFT_Tab(self.output_text_edit)
        IO_tab.addTab(sub_tab5, 'ORCA out to ORCA TD-DFT')

        main_tabs.addTab(IO_tab, 'Input/Output operations')

        #Sub-tabs for output analysis, then add to main output analysis tab
        sub_tab6 = ORCA_Optim_plot_Tab(self.output_text_edit)
        output_analysis_tab.addTab(sub_tab6, 'Plot ORCA Opt Trajectory')   

        sub_tab7 = DFTThermochemTab(self.output_text_edit)
        output_analysis_tab.addTab(sub_tab7, 'Thermochem Analyzer')

        sub_tab8 = CCSDTTab(self.output_text_edit)
        output_analysis_tab.addTab(sub_tab8, 'CCSDT Analyzer')

        sub_tab9 = Extract_IR_SpectraTab(self.output_text_edit)
        output_analysis_tab.addTab(sub_tab9, 'Extract Vib Spectra')

        sub_tab10 = Extract_TDDFT_VG_SpectraTab(self.output_text_edit)
        output_analysis_tab.addTab(sub_tab10, 'Extract TD-DFT Spectra (VG)')       

        main_tabs.addTab(output_analysis_tab, 'Output analysis')

        #Sub tabs for additional analyses
        sub_tab11 = BWCCS_Tab(self.output_text_edit)
        extras_tab.addTab(sub_tab11, 'Boltzmann-weight CCSs')

        sub_tab12 = LED_Analysis_Tab(self.output_text_edit)
        extras_tab.addTab(sub_tab12, 'Analyze LED')

        sub_tab13 = NEB_Analysis_Tab(self.output_text_edit)
        extras_tab.addTab(sub_tab13, 'Analyze NEB')             
        
        main_tabs.addTab(extras_tab, 'Additional ORCA functions')

        #Add all the tabs to the GUI
        main_layout.addWidget(main_tabs)

        #Create and add a window for print statements to show up within the GUI interface
        output_title_label = QLabel("Code printout window")
        output_title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        main_layout.addWidget(output_title_label)

        self.output_text_edit.setMinimumHeight(150)  

        main_layout.addWidget(self.output_text_edit)
        main_widget.setLayout(main_layout)
        
        self.setCentralWidget(main_widget)

        sys.stdout = self.text_redirector

    def update_output_text(self, text):
        self.output_text_edit.insertPlainText(text)

    def check_for_update(self):

        #The GUI is likely to the updated throughout the years, so its best practice to implement some update functionality - users may not check GitHub frequently. 
        def get_latest_commit_sha(repo_url, branch='HEAD'):
            '''Grabs the SHA value associated with the latest commit to a GitHub repo. Function takes a URL as input.''' 
            try:
                # Run the git ls-remote command
                result = subprocess.run(['git', 'ls-remote', repo_url, branch], capture_output=True, text=True)

                # Check if the command was successful
                if result.returncode != 0:
                    raise Exception(f'Error: {result.stderr}')

                # Parse the output to get the SHA
                output = result.stdout.split()
                return output[0] if output else None
            
            except Exception as e:
                raise Exception(f'An error occurred: {e}')

        def delete_dir(directory):
            '''Windows has special permissions on read-only folders - here's a function that deals with it using shutil.rmtree'''

            os_type = platform.system()

            def handleRemoveReadonly(func, path, exc_info):
                '''A function to deal with the error encourntered when trying to delete a file using shutil.rmtree that is read-only in Windows.'''
                
                #from: https://stackoverflow.com/questions/4829043/how-to-remove-read-only-attrib-directory-with-python-in-windows
                os.chmod(path, stat.S_IWRITE)
                os.unlink(path)

            try:
                if os_type == 'Windows':
                    shutil.rmtree(directory, onexc=handleRemoveReadonly)
                
                #no special read/write perms needed to delete read only folders on Mac/Linux platforms
                else:
                    shutil.rmtree(directory)
            
            #If using Python versions <3.9, syntax for encountering errors when deleting directories is different.
            except TypeError:

                try:
                    shutil.rmtree(directory, onerror=handleRemoveReadonly) 

                except:
                    print(f'Encountered error when attempting to delete {directory}. This is probably because your python version is <3.9. Please update to Python 3.12+ to rectify this - you can also manually delete the temp folder after the GUI loads. ')
                    pass

            except PermissionError:
                print(f'Encountered permission error when attempting to delete {directory}. This is probably caused by a OneDrive sync issue - you can manually delete the GUI/temp folder after the GUI loads. ')
                pass

            except Exception as e: 
                print(f'An unexpected error encountered trying to remove {directory}: {e}.\n Please report this error and your workflow / system specs to the Issues section on Github. You can delete still delete the /temp manually.')         
                # Get the current working directory and define the temporary directory path
                pass

        '''Main Update function'''

        root = os.getcwd()

        #URL of the MobCal-MPI repo
        repo_url = 'https://github.com/HopkinsLaboratory/ORCA_Analysis_GUI'

        #get SHA value of repo-ID
        repo_SHA = get_latest_commit_sha(repo_url)

        # File to store the local version's SHA value
        ID_file = os.path.join(root, 'ID.txt') 
        try:
            with open(ID_file,'r') as opf:
                file_content = opf.read()

                # Check if the file is in the correct format
                if '\n' in file_content or '\r' in file_content or ' ' in file_content:
                    print(f'{ID_file} is not in the correct format. Please re-download this file from the ORCA Analysis GUI GitHub repo and re-run the GUI launcher.')
                    return

                local_SHA = file_content.strip()

        except FileNotFoundError:
            print(f'{ID_file} could not be found. Please re-download from the ORCA Analysis GUI GitHub repo and re-run the GUI launcher.')
            return
            
        # Remove the temporary directory if it exists and only if the local and repo SHAs match
        if repo_SHA == local_SHA: 
            temp_dir = os.path.join(os.getcwd(), 'temp')
            if os.path.isdir(temp_dir):
                delete_dir(temp_dir)

        else:
            choice_title = 'Update available!'
            choice_prompt = 'An update to the ORCA Analysis GUI is available. Would you like to update now?'
            choice = QMessageBox.question(self, choice_title, choice_prompt, QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)

            if choice == QMessageBox.StandardButton.Yes:
                print('The ORCA Analysis GUI is being updated. Any errors encountered during the update process will be printed below.')
                
                #Run the update function
                Update_GUI_files(repo_url, root, ID_file, repo_SHA, delete_dir)

                print('The ORCA Analysis GUI files have been succesfully updated to their current version. The GUI will now close. Please reload the GUI')
                sys.exit(0)

            else:
                print('The user has opted to use their local version of ORCA Analysis GUI.')

class XYZFileSplitterTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
      
        layout = QVBoxLayout(self)

        #XYZ file input section - title, input box, and explorer option for browsing xyz files
        file_layout = QHBoxLayout()

        file_label = QLabel('.xyz file:')
        self.file_path_input = QLineEdit()
        file_button = QPushButton('Browse')
        
        file_layout.addWidget(file_label)
        file_layout.addSpacing(55) #spacing for alignment of input boxes
        file_layout.addWidget(self.file_path_input)
        file_layout.addWidget(file_button)

        file_button.clicked.connect(self.browse_xyz_file)  #Connect the button click to the browse_xyz_file function
       
        layout.addLayout(file_layout)
        layout.addSpacing(5)  

        #Basenanme input section
        basename_layout = QHBoxLayout()
        basename_label = QLabel('Output basename:')
        self.basename_input = QLineEdit()
        self.basename_input.setPlaceholderText('e.g., "Fluoxetine" yields Fluoxetine_1.gjf, Fluoxetine_2.gjf, etc.')

        basename_layout.addWidget(basename_label)
        basename_layout.addWidget(self.basename_input)

        layout.addLayout(basename_layout)
        layout.addSpacing(10)
        
        #Run Button
        run_button = QPushButton('Run')
        layout.addWidget(run_button)
        run_button.clicked.connect(self.run_xyz_file_split)

        layout.addStretch(1)

        layout.setContentsMargins(30, 30, 30, 30) 
        self.setLayout(layout)

    def browse_xyz_file(self):

        file_path, _ = QFileDialog.getOpenFileName(self, 'Select .xyz file', '', 'XYZ Files (*.xyz)')

        #Check if a file path was selected
        if file_path:
            self.file_path_input.setText(file_path)

    def run_xyz_file_split(self):
        
        file_path = self.file_path_input.text()
        basename = self.basename_input.text()

        try:
            xyz_file_splitter(file_path, basename)
        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

class GJFtoORCAInputTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        #Directory selection
        directory_layout = QHBoxLayout()
        
        directory_label = QLabel('Directory:')
        self.directory_input = QLineEdit()
        directory_button = QPushButton('Browse')
        directory_button.clicked.connect(self.browse_directory)

        directory_layout.addWidget(directory_label)
        directory_layout.addWidget(self.directory_input)
        directory_layout.addWidget(directory_button)
        
        layout.addLayout(directory_layout)
        layout.addSpacing(5)  

        #Memory, Cores, Charge, Multiplicity
        performance_layout = QHBoxLayout()
        mpp_label = QLabel('Mem per core (MB):')
        self.mpp_input = QSpinBox()
        self.mpp_input.setMaximum(20000)
        self.mpp_input.setValue(3400)
        self.mpp_input.setMinimumWidth(60)

        ncores_label = QLabel('#cores:')
        self.ncores_input = QSpinBox()
        self.ncores_input.setMaximum(32)
        self.ncores_input.setValue(8)
        self.ncores_input.setMinimumWidth(60)

        charge_label = QLabel('Charge:')
        self.charge_input = QSpinBox()
        self.charge_input.setValue(1)
        self.charge_input.setMinimumWidth(60)

        multiplicity_label = QLabel('Multiplicity:')
        self.multiplicity_input = QSpinBox()
        self.multiplicity_input.setValue(1)
        self.multiplicity_input.setMinimumWidth(60)

        performance_layout.addWidget(mpp_label)
        performance_layout.addWidget(self.mpp_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(ncores_label)
        performance_layout.addWidget(self.ncores_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(charge_label)
        performance_layout.addWidget(self.charge_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(multiplicity_label)
        performance_layout.addWidget(self.multiplicity_input)
        performance_layout.addStretch(1)
        
        layout.addLayout(performance_layout)
        layout.addSpacing(5)

        #Calculation Method         
        calc_type_label = QLabel('Method:')
        layout.addWidget(calc_type_label)

        #Selection box for DFT/CCSDT/Custom methods
        calc_layout = QHBoxLayout()
        self.calc_type_combobox = QComboBox()
        self.calc_type_combobox.addItems(['DFT', 'CCSDT', 'Custom'])
        self.calc_type_combobox.setCurrentIndex(0)
        self.calc_type_combobox.currentIndexChanged.connect(self.update_calc_line)

        calc_layout.addWidget(self.calc_type_combobox)
        calc_layout.addSpacing(10)

        #calc_line_label = QLabel('Calculation Line (start with !):')
        self.calc_line_input = QLineEdit()
        self.calc_line_input.setText('! wB97X-D3 TightOpt Freq def2-TZVPP def2/J RIJCOSX TightSCF defgrid3') #default DFT method 
        calc_layout.addWidget(self.calc_line_input)

        layout.addLayout(calc_layout)
        layout.addSpacing(5)

        #ESP Charges Checkbox
        grid_layout = QHBoxLayout()
        self.esp_charges_checkbox = QCheckBox('Compute ESP Charges?')
        self.esp_charges_checkbox.stateChanged.connect(self.toggle_esp_charge_inputs)
        
        #Grid and Rmax
        self.grid_label = QLabel('Grid:')
        self.grid_input = QDoubleSpinBox()
        self.grid_input.setValue(0.1)
        self.grid_input.setEnabled(False)

        self.rmax_label = QLabel('Rmax:')
        self.rmax_input = QDoubleSpinBox()
        self.rmax_input.setValue(3.0)
        self.rmax_input.setEnabled(False)

        grid_layout.addWidget(self.esp_charges_checkbox)
        grid_layout.addSpacing(30)
        grid_layout.addWidget(self.grid_label)
        grid_layout.addWidget(self.grid_input)
        grid_layout.addSpacing(10)
        grid_layout.addWidget(self.rmax_label)
        grid_layout.addWidget(self.rmax_input)
        grid_layout.addStretch(1)
        
        layout.addLayout(grid_layout)

        #Checkboxes for Hessian, Polarization, XYZ Coordinates
        self.calc_hess_checkbox = QCheckBox('Calculate Hessian on first optimization step?')
        self.polarization_checkbox = QCheckBox('Calculate dipole/quadrupole moments?')
        self.write_xyz_checkbox = QCheckBox('Call XYZ coordinates from external .xyz file?')

        layout.addWidget(self.calc_hess_checkbox)
        layout.addWidget(self.polarization_checkbox)
        layout.addWidget(self.write_xyz_checkbox)
        layout.addSpacing(10)

        #Run button
        run_button = QPushButton('Run')
        run_button.clicked.connect(self.run_gjf_to_orca_input)

        layout.addWidget(run_button)
        layout.addStretch(1)

        layout.setContentsMargins(30, 30, 30, 30)
        self.setLayout(layout)

    def update_calc_line(self, index):
        if index == 0:  #DFT
            self.calc_line_input.setText("! wB97X-D3 TightOpt Freq def2-TZVPP def2/J RIJCOSX TightSCF defgrid3 CHELPG")
            self.calc_line_input.setReadOnly(True)
        elif index == 1:  #CCSDT
            self.calc_line_input.setText("! DLPNO-CCSD(T) def2-TZVPP def2-TZVP/C VeryTightSCF")
            self.calc_line_input.setReadOnly(True)
        else:  #Custom
            self.calc_line_input.setReadOnly(False)
    
    def toggle_esp_charge_inputs(self, state):
        self.grid_input.setEnabled(state == 2)
        self.rmax_input.setEnabled(state == 2)

        #Check if the ESP charges checkbox is checked
        if state == 2:  #Qt.Checked state is 2
            current_text = self.calc_line_input.text()
            
            #Check if 'CHELPG' or 'chelpg' is not in the current method line if the box is checked. If it's not, add it. 
            if 'CHELPG' not in current_text.upper():  
                updated_text = f'{current_text.strip()} CHELPG'
                self.calc_line_input.setText(updated_text)
    
    def browse_directory(self):
        directory_path = QFileDialog.getExistingDirectory(self, "Select Directory")

        if directory_path:
            self.directory_input.setText(directory_path)

    def run_gjf_to_orca_input(self):
        try:
            #Gather input values from the GUI
            directory = self.directory_input.text()
            mpp = self.mpp_input.value()
            ncores = self.ncores_input.value()
            charge = self.charge_input.value()
            multiplicity = self.multiplicity_input.value()
            calc_line = self.calc_line_input.text()

            #Ensure calc_line starts with "!"
            if not calc_line.startswith("!"):
                calc_line = "! " + calc_line

            #Check the state of ESP_Charges checkbox
            esp_charges_checked = self.esp_charges_checkbox.isChecked()
            
            #If ESP_Charges is checked, gather values from additional fields
            if esp_charges_checked:
                grid = self.grid_input.value()
                rmax = self.rmax_input.value()
            else:
                grid = 0.1
                rmax = 3.0

            #checkbox logic
            if self.calc_hess_checkbox.isChecked(): self.calc_hess_checked = True
            else: self.calc_hess_checked = False

            if self.polarization_checkbox.isChecked(): self.polarization_checked = True
            else: self.polarization_checked = False

            if self.write_xyz_checkbox.isChecked(): self.write_xyz_checked = True
            else: self.write_xyz_checked = False

            try:
                #Call the Gaussian_gjf_to_ORCA_input function
                Gaussian_gjf_to_ORCA_input(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, self.calc_hess_checked, self.polarization_checked, self.write_xyz_checked)
            except Exception as e:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

class CosineSimTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)

        #Directory input
        directory_layout = QHBoxLayout()

        directory_label = QLabel('Directory or CSV File:')
        self.directory_input = QLineEdit()

        browse_button = QPushButton('Browse')
        browse_button.clicked.connect(self.browse_directory_or_csv)
        
        directory_layout.addWidget(directory_label)
        directory_layout.addWidget(self.directory_input)
        directory_layout.addWidget(browse_button)

        layout.addLayout(directory_layout)  
        layout.addSpacing(5)

        #Similiarity input
        sim_layout = QHBoxLayout()

        #add the input for defining the similiarity threshold
        self.sim_label = QLabel('Similiarity threshold (%):')
        self.sim_input = QDoubleSpinBox()
        self.sim_input.setValue(97.0)
        self.sim_input.setMaximum(100.0)  
        self.sim_input.setMinimumWidth(70)
        
        sim_layout.addWidget(self.sim_label)
        sim_layout.addWidget(self.sim_input)
        sim_layout.addStretch(1)

        layout.addLayout(sim_layout)
        layout.addSpacing(5)

        #add a checkbox for writing all pairwise similarities to a file
        self.write_sim_checkbox = QCheckBox('Write pairwise similiarities to a .csv file?')
        layout.addWidget(self.write_sim_checkbox)

        layout.addSpacing(10)

        #Run Button
        run_button = QPushButton('Run Cosine Similarity')
        run_button.clicked.connect(self.run_cosine_similarity)
        
        layout.addWidget(run_button)
        layout.addStretch(1)

        layout.setContentsMargins(30, 30, 30, 30) 
        self.setLayout(layout)


    def browse_directory_or_csv(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter('CSV Files (*.csv)')
        file_dialog.setDefaultSuffix('csv')
        file_dialog.setFileMode(QFileDialog.FileMode.ExistingFile)

        file_path, _ = file_dialog.getOpenFileName()  #Use getOpenFileName to get the file path
        
        if file_path:
            self.directory_input.setText(file_path)

    def run_cosine_similarity(self):
        #Get the selected directory or CSV file, as well as the similarity threshold
        input_path = self.directory_input.text()
        similarity_threshold = self.sim_input.value()

        #checkbox logic
        if self.write_sim_checkbox.isChecked(): self.write_sim = True
        else: self.write_sim = False

        try:
            if not os.path.isdir(input_path) and not os.path.isfile(input_path):
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {input_path} is not a valid directory or file. Please check that  you have specified is file/location that exists, then try again.')
            
            else:
                #if proper files are given, run the code
                cosine_sim(input_path, similarity_threshold, self.write_sim)

        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

class ORCAOut_ORCAInputTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        #Directory selection
        directory_layout = QHBoxLayout()
        directory_label = QLabel('Directory:')
        self.directory_input = QLineEdit()
        directory_button = QPushButton('Browse')
        directory_button.clicked.connect(self.browse_directory)

        directory_layout.addWidget(directory_label)
        directory_layout.addWidget(self.directory_input)
        directory_layout.addWidget(directory_button)
        layout.addLayout(directory_layout)
        layout.addSpacing(5) 

        #Memory, Cores, Charge, Multiplicity
        performance_layout = QHBoxLayout()
        mpp_label = QLabel('Mem per core (MB):')
        self.mpp_input = QSpinBox()
        self.mpp_input.setMaximum(20000)
        self.mpp_input.setValue(3400)
        self.mpp_input.setMinimumWidth(60)

        ncores_label = QLabel('#cores:')
        self.ncores_input = QSpinBox()
        self.ncores_input.setMaximum(32)
        self.ncores_input.setValue(8)
        self.ncores_input.setMinimumWidth(60)

        charge_label = QLabel('Charge:')
        self.charge_input = QSpinBox()
        self.charge_input.setValue(1)
        self.charge_input.setMinimumWidth(60)

        multiplicity_label = QLabel('Multiplicity:')
        self.multiplicity_input = QSpinBox()
        self.multiplicity_input.setValue(1)
        self.multiplicity_input.setMinimumWidth(60)

        performance_layout.addWidget(mpp_label)
        performance_layout.addWidget(self.mpp_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(ncores_label)
        performance_layout.addWidget(self.ncores_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(charge_label)
        performance_layout.addWidget(self.charge_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(multiplicity_label)
        performance_layout.addWidget(self.multiplicity_input)
        performance_layout.addStretch(1)
        
        layout.addLayout(performance_layout)
        layout.addSpacing(5)

        #Calculation Method         
        calc_type_label = QLabel('Method:')
        layout.addWidget(calc_type_label)

        #Selection box for DFT/CCSDT/Custom methods
        calc_layout = QHBoxLayout()
        self.calc_type_combobox = QComboBox()
        self.calc_type_combobox.addItems(['DFT', 'CCSDT', 'Custom'])
        self.calc_type_combobox.setCurrentIndex(0)
        self.calc_type_combobox.currentIndexChanged.connect(self.update_calc_line)

        calc_layout.addWidget(self.calc_type_combobox)
        calc_layout.addSpacing(10)

        #calc_line_label = QLabel('Calculation Line (start with !):')
        self.calc_line_input = QLineEdit()
        self.calc_line_input.setText('! wB97X-D3 TightOpt Freq def2-TZVPP def2/J RIJCOSX TightSCF defgrid3 CHELPG') #default DFT method 
        calc_layout.addWidget(self.calc_line_input)
        #calc_layout.addStretch(1)

        layout.addLayout(calc_layout)
        layout.addSpacing(5)

        #ESP Charges Checkbox
        grid_layout = QHBoxLayout()
        self.esp_charges_checkbox = QCheckBox('Compute ESP Charges?')
        self.esp_charges_checkbox.stateChanged.connect(self.toggle_esp_charge_inputs)
        
        #Grid and Rmax
        self.grid_label = QLabel('Grid:')
        self.grid_input = QDoubleSpinBox()
        self.grid_input.setValue(0.1)
        self.grid_input.setEnabled(False)

        self.rmax_label = QLabel('Rmax:')
        self.rmax_input = QDoubleSpinBox()
        self.rmax_input.setValue(3.0)
        self.rmax_input.setEnabled(False)

        grid_layout.addWidget(self.esp_charges_checkbox)
        grid_layout.addSpacing(30)
        grid_layout.addWidget(self.grid_label)
        grid_layout.addWidget(self.grid_input)
        grid_layout.addSpacing(10)
        grid_layout.addWidget(self.rmax_label)
        grid_layout.addWidget(self.rmax_input)
        grid_layout.addStretch(1)
        
        layout.addLayout(grid_layout)

        #Checkboxes for Hessian, Polarization, XYZ Coordinates
        self.calc_hess_checkbox = QCheckBox('Calculate Hessian on first optimization step?')
        self.polarization_checkbox = QCheckBox('Calculate dipole/quadrupole Moments?')
        self.write_xyz_checkbox = QCheckBox('Call XYZ coordinates from external .xyz file?')
        self.write_gjf_checkbox = QCheckBox('Write final coordinates to a .gjf file for visualization in GaussView?')

        layout.addWidget(self.calc_hess_checkbox)
        layout.addWidget(self.polarization_checkbox)
        layout.addWidget(self.write_xyz_checkbox)
        layout.addWidget(self.write_gjf_checkbox)
        layout.addSpacing(10)

        #Run button
        run_button = QPushButton('Run')
        run_button.clicked.connect(self.run_ORCA_out_to_ORCA_inp)
        layout.addWidget(run_button)

        layout.addStretch(1)

        layout.setContentsMargins(30, 30, 30, 30)
        self.setLayout(layout)

    def update_calc_line(self, index):
        #Options for calculation line
        if index == 0:  #DFT
            self.calc_line_input.setText('! wB97X-D3 TightOpt Freq def2-TZVPP def2/J RIJCOSX TightSCF defgrid3')
            self.calc_line_input.setReadOnly(True)
        elif index == 1:  #CCSDT
            self.calc_line_input.setText('! DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C VeryTightSCF')
            self.calc_line_input.setReadOnly(True)
        else:  #Custom
            self.calc_line_input.setReadOnly(False)
    
    def toggle_esp_charge_inputs(self, state):
        self.grid_input.setEnabled(state == 2)
        self.rmax_input.setEnabled(state == 2)

        #Check if the ESP charges checkbox is checked
        if state == 2:  #Qt.Checked state is 2
            current_text = self.calc_line_input.text()
            
            #Check if 'CHELPG' or 'chelpg' is not in the current text if the box is checked. If it's not, add it. 
            if 'CHELPG' not in current_text.upper():  
                updated_text = f'{current_text.strip()} CHELPG'
                self.calc_line_input.setText(updated_text)
    
    def browse_directory(self):
        directory_path = QFileDialog.getExistingDirectory(self, "Select Directory")

        if directory_path:
            self.directory_input.setText(directory_path)

    def run_ORCA_out_to_ORCA_inp(self):
        try:
            #Gather input values from the GUI
            directory = self.directory_input.text()
            mpp = self.mpp_input.value()
            ncores = self.ncores_input.value()
            charge = self.charge_input.value()
            multiplicity = self.multiplicity_input.value()
            calc_line = self.calc_line_input.text()

            #Ensure calc_line starts with "!"
            if not calc_line.startswith("!"):
                calc_line = "! " + calc_line

            #Check the state of ESP_Charges checkbox
            esp_charges_checked = self.esp_charges_checkbox.isChecked()
            
            #If ESP_Charges is checked, gather values from additional fields
            if esp_charges_checked:
                grid = self.grid_input.value()
                rmax = self.rmax_input.value()
            else:
                #default values for the ESP grid
                grid = 0.1
                rmax = 3.0

            #checkbox logic
            if self.calc_hess_checkbox.isChecked(): self.calc_hess_checked = True
            else: self.calc_hess_checked = False

            if self.polarization_checkbox.isChecked(): self.polarization_checked = True
            else: self.polarization_checked = False

            if self.write_xyz_checkbox.isChecked(): self.write_xyz_checked = True
            else: self.write_xyz_checked = False

            if self.write_gjf_checkbox.isChecked(): self.write_gjf_checked = True
            else: self.write_gjf_checked = False

            try:
                #Call the Gaussian_gjf_to_ORCA_input function
                ORCA_out_to_ORCA_inp(directory, mpp, ncores, charge, multiplicity, calc_line, esp_charges_checked, grid, rmax, self.calc_hess_checked, self.polarization_checked, self.write_xyz_checked, self.write_gjf_checked)
            except Exception as e:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

class ORCAOut_ORCA_TDDFT_Tab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        #Directory selection title
        directory_label = QLabel('Directory containing .out and .hess files:')
        layout.addWidget(directory_label)

        #Directory input box
        directory_layout = QHBoxLayout()
        self.directory_input = QLineEdit()
        directory_button = QPushButton('Browse')
        directory_button.clicked.connect(self.browse_directory)

        directory_layout.addWidget(self.directory_input)
        directory_layout.addWidget(directory_button)
        layout.addLayout(directory_layout)
        layout.addSpacing(5) 

        #Memory, Cores, Charge, Multiplicity
        performance_layout = QHBoxLayout()
        mpp_label = QLabel('Mem per core (MB):')
        self.mpp_input = QSpinBox()
        self.mpp_input.setMaximum(20000)
        self.mpp_input.setValue(3400)
        self.mpp_input.setMinimumWidth(60)

        ncores_label = QLabel('#cores:')
        self.ncores_input = QSpinBox()
        self.ncores_input.setMaximum(32)
        self.ncores_input.setValue(8)
        self.ncores_input.setMinimumWidth(60)

        charge_label = QLabel('Charge:')
        self.charge_input = QSpinBox()
        self.charge_input.setValue(1)
        self.charge_input.setMinimumWidth(60)

        multiplicity_label = QLabel('Multiplicity:')
        self.multiplicity_input = QSpinBox()
        self.multiplicity_input.setValue(1)
        self.multiplicity_input.setMinimumWidth(60)

        performance_layout.addWidget(mpp_label)
        performance_layout.addWidget(self.mpp_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(ncores_label)
        performance_layout.addWidget(self.ncores_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(charge_label)
        performance_layout.addWidget(self.charge_input)
        performance_layout.addSpacing(10)

        performance_layout.addWidget(multiplicity_label)
        performance_layout.addWidget(self.multiplicity_input)
        performance_layout.addStretch(1)
        
        layout.addLayout(performance_layout)
        layout.addSpacing(5)

        #calc line
        calc_layout = QHBoxLayout()
        calc_line_label = QLabel('Method:')
        self.calc_line_input = QLineEdit()
        self.calc_line_input.setText('! wB97X-D3 def2-TZVPP defgrid3 ESD(ABS)') #default DFT method 

        calc_layout.addWidget(calc_line_label)
        calc_layout.addWidget(self.calc_line_input)

        layout.addLayout(calc_layout)

        #number of states
        states_layout = QHBoxLayout()
        n_states_label = QLabel('#states:')
        self.nstates_input = QSpinBox()
        self.nstates_input.setValue(10)
        self.nstates_input.setMinimumWidth(50)

        states_layout.addWidget(n_states_label)
        states_layout.addWidget(self.nstates_input)
        states_layout.addStretch(1)

        layout.addLayout(states_layout)
        layout.addSpacing(5)

        #write .gjf checkbox
        self.write_gjf_checkbox = QCheckBox('Write final coordinates to a .gjf file for visualization in GaussView?')
        layout.addWidget(self.write_gjf_checkbox)

        #Add some vertical spacing
        layout.addSpacing(10)

        run_button = QPushButton('Run')
        run_button.clicked.connect(self.run_ORCA_out_to_ORCA_TDDFT_VG)

        layout.addWidget(run_button)

        #Set vertical and horizontal spacing for margins
        layout.setContentsMargins(30, 30, 30, 30)  
        
        self.setLayout(layout)

        #Add some vertical spacing
        layout.addStretch(1)
   
    def browse_directory(self):
        directory_path = QFileDialog.getExistingDirectory(self, "Select Directory")

        if directory_path:
            self.directory_input.setText(directory_path)

    def run_ORCA_out_to_ORCA_TDDFT_VG(self):
        try:
            #Gather input values from the GUI
            directory = self.directory_input.text()
            mpp = self.mpp_input.value()
            ncores = self.ncores_input.value()
            charge = self.charge_input.value()
            multiplicity = self.multiplicity_input.value()
            calc_line = self.calc_line_input.text()
            states = self.nstates_input.value()

            #checkbox logic for .gjf writing
            if self.write_gjf_checkbox.isChecked(): self.write_gjf_checked = True
            else: self.write_gjf_checked = False

            try:
                #Call the ORCA .out to ORCA .inp for Vertical Gradient (VG) TD-DFT calcs
                ORCA_out_to_ORCA_TDDFT_VG(directory, mpp, ncores, charge, multiplicity, calc_line, states, self.write_gjf_checked)
            except Exception as e:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

class ORCA_Optim_plot_Tab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)

        #Layout for specifying the output file to plot
        directory_layout = QHBoxLayout()

        directory_label = QLabel('ORCA.out file:')
        directory_layout.addWidget(directory_label)

        self.file_input = QLineEdit()
        directory_layout.addWidget(self.file_input)

        browse_button = QPushButton('Browse')
        browse_button.clicked.connect(self.browse_OUT_file)
        directory_layout.addWidget(browse_button)

        layout.addLayout(directory_layout)
        layout.addSpacing(5)

        #Save optim plot checkbox
        self.save_plot_checkbox = QCheckBox('Save the plot to an external file?')
        layout.addWidget(self.save_plot_checkbox)
        layout.addSpacing(10)

        #Run Button
        run_button = QPushButton('Plot')
        run_button.clicked.connect(self.run_plt)
        layout.addWidget(run_button)

        layout.addStretch(1)
        layout.setContentsMargins(30, 30, 30, 30) 
        self.setLayout(layout)

    def browse_OUT_file(self):
        #Open a file dialog to select a .out file
        file_path, _ = QFileDialog.getOpenFileName(self, "Select .out file", "", "OUT Files (*.out)")

        #Check if a file was selected
        if file_path:
            self.file_input.setText(file_path)

    def run_plt(self):
        #Get the selected directory or CSV file
        input_path = self.file_input.text()

        #Set save plot option depending on the checkbox
        if self.save_plot_checkbox.isChecked(): save_plot = True
        else: save_plot = False

        if not os.path.isfile(input_path):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {input_path} is not a valid .out file. Please check that the directory you have specified is correct, then try again.')
            return
            
        else:
            ORCA_opt_plt(input_path, save_plot)

class DFTThermochemTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)

        #Directory layout
        directory_layout = QHBoxLayout()

        directory_label = QLabel('Directory:')
        self.directory_input = QLineEdit()
        browse_button = QPushButton('Browse')
        browse_button.clicked.connect(self.browse_directory)
        
        directory_layout.addWidget(directory_label)
        directory_layout.addWidget(self.directory_input)
        directory_layout.addWidget(browse_button)
        
        layout.addLayout(directory_layout)
        layout.addSpacing(5)

        #Temperature, pressure, and vibrational scaling factor input
        temp_pressure_vib_layout = QHBoxLayout()
        
        temp_label = QLabel('Temp (K):')
        self.temp_input = QDoubleSpinBox()
        self.temp_input.setRange(0.0, 100000.0)
        self.temp_input.setSingleStep(10)
        self.temp_input.setValue(298.15)
        self.temp_input.setMinimumWidth(60)  
        
        pressure_label = QLabel('Pressure (Pa):')
        self.pressure_input = QDoubleSpinBox()
        self.pressure_input.setRange(0.0, 100000000.0) 
        self.pressure_input.setSingleStep(100)    
        self.pressure_input.setValue(100000.0)
        self.pressure_input.setMinimumWidth(60)      

        vib_scaling_label = QLabel('Vib. Scaling Factor:')
        self.vib_scaling_input = QDoubleSpinBox()
        self.vib_scaling_input.setRange(0, 10)  
        self.vib_scaling_input.setSingleStep(0.01)    
        self.vib_scaling_input.setValue(1.00)
        self.vib_scaling_input.setMinimumWidth(60)        

        temp_pressure_vib_layout.addWidget(temp_label)
        temp_pressure_vib_layout.addWidget(self.temp_input)
        temp_pressure_vib_layout.addSpacing(10)

        temp_pressure_vib_layout.addWidget(pressure_label)
        temp_pressure_vib_layout.addWidget(self.pressure_input)
        temp_pressure_vib_layout.addSpacing(10)

        temp_pressure_vib_layout.addWidget(vib_scaling_label)
        temp_pressure_vib_layout.addWidget(self.vib_scaling_input)
        temp_pressure_vib_layout.addStretch(1)
        
        layout.addLayout(temp_pressure_vib_layout)
        layout.addSpacing(5)

        #Sorting Options
        sort_by_layout = QHBoxLayout()
        sort_by_layout.addWidget(QLabel('Sort files by:'))
        self.sort_by_combo = QComboBox()
        self.sort_by_combo.addItems(['Filename', 'Gibbs (G)', 'Enthalpy (H)', 'Entropy (S)', 'Thermal Energy (E)', 'Zero-point energy (ZPE)'])
        sort_by_layout.addWidget(self.sort_by_combo)
        sort_by_layout.addStretch()  #Needed to align text with box on the left
        layout.addLayout(sort_by_layout)
        layout.addSpacing(10)

        #Run Button
        run_button = QPushButton('Extract Thermochemistry')
        run_button.clicked.connect(self.run_dft_thermochemistry)
        
        layout.addWidget(run_button)
        layout.addStretch(1)

        layout.setContentsMargins(30, 30, 30, 30)
        self.setLayout(layout)

    def browse_directory(self):
        directory_path = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory_path:
            self.directory_input.setText(directory_path)

    def run_dft_thermochemistry(self):
        input_path = self.directory_input.text()
        
        #check that T, p, and vib scaling factor inputs are numeric inputs and are all non-negative numbers
        try:
            T = float(self.temp_input.text())
            if T < 0:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: Temperature cannot be negative.')
                return   
        except ValueError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: Temperature input must be numeric.')
            return

        try:
            p = float(self.pressure_input.text())
            if p < 0:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: Pressure cannot be negative.')
                return    
        except ValueError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: Pressure input must be numeric.')
            return

        try:
            vib_scl = float(self.vib_scaling_input.text())

            if vib_scl < 0:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: Vibrational scaling factor cannot be negative.')
                return                
        except ValueError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: Vibrational scaling factor input must be numeric.')
            return
        
        #get sorting option
        sort_option = self.sort_by_combo.currentText()
        sort_by_map = {'Filename': 'F', 'Gibbs': 'G', 'Enthalpy (H)': 'H', 'Entropy (S)': 'S', 'Thermal Energy (E)': 'E', 'Zero-point energy': 'ZPE'}
        sort_by = sort_by_map.get(sort_option, 'F')

        #check that the input path is a real directory, and if so, run the thermochem calculator code
        try:
            if not os.path.isdir(input_path):
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {input_path} is not a valid directory. Please check that the directory you have specified is correct, then try again.')
            else:
                if sort_option == 'Filename':
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Thermochem is being analyzed at {T} K, {p} Pa. Output files will be sorted by energy according to the order of the {sort_option}s.')
                else:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Thermochem is being analyzed at {T} K, {p} Pa. Output files will be sorted by energy according to the order of the {sort_option} energies.')
                
                #run the code
                ORCA_Thermochem_Calculator(input_path, T, p, vib_scl, sort_by)

        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\n \nTraceback: {traceback.format_exc()}')

class CCSDTTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)

        #Directory Input
        directory_layout = QHBoxLayout()  

        directory_label = QLabel('Directory:')
        self.directory_input = QLineEdit()
        
        directory_layout.addWidget(directory_label)
        directory_layout.addWidget(self.directory_input)

        browse_button = QPushButton('Browse')
        browse_button.clicked.connect(self.browse_directory)
        directory_layout.addWidget(browse_button)

        layout.addLayout(directory_layout)
        layout.addSpacing(10)

        #Run Button
        run_button = QPushButton('Extract CCSDT')
        run_button.clicked.connect(self.run_CCSDT)
        layout.addWidget(run_button)

        layout.addStretch(1)
        layout.setContentsMargins(30, 30, 30, 30)
        self.setLayout(layout)

    def browse_directory(self):
        directory_path = QFileDialog.getExistingDirectory(self, "Select Directory")

        if directory_path:
            self.directory_input.setText(directory_path)

    def run_CCSDT(self):
        #Get the selected directory or CSV file
        input_path = self.directory_input.text()

        try:
            if not os.path.isdir(input_path):
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {input_path} is not a valid directory. Please check that the directory you have specified is correct, then try again.')
            else:
                ORCA_CCSDT(input_path)

        except Exception as e:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: {e}\nTraceback: {traceback.format_exc()}')

class Extract_IR_SpectraTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)

        #Directory input
        directory_layout = QHBoxLayout()
        
        directory_label = QLabel('Directory:')
        self.directory_input = QLineEdit()
        
        browse_button = QPushButton('Browse')
        browse_button.clicked.connect(self.browse_directory)

        directory_layout.addWidget(directory_label)
        directory_layout.addWidget(self.directory_input)
        directory_layout.addWidget(browse_button)
        
        layout.addLayout(directory_layout)
        layout.addSpacing(5)

        #Lower and Upper Bound Inputs
        bounds_layout = QHBoxLayout()
        
        lower_bound_label = QLabel('Export from:')
        self.lower_bound_input = QDoubleSpinBox()
        self.lower_bound_input.setRange(0.0, 5500.0)
        self.lower_bound_input.setSingleStep(10)
        self.lower_bound_input.setValue(400.0)
        self.lower_bound_input.setSuffix(' cm-1') 
        self.lower_bound_input.setDecimals(1)
        self.lower_bound_input.setMinimumWidth(100)  
        
        upper_bound_label = QLabel('to:')
        self.upper_bound_input = QDoubleSpinBox()
        self.upper_bound_input.setRange(1.0, 5500.0) 
        self.upper_bound_input.setSingleStep(10)    
        self.upper_bound_input.setValue(4000.0)
        self.upper_bound_input.setSuffix(' cm-1')
        self.upper_bound_input.setDecimals(1) 
        self.upper_bound_input.setMinimumWidth(100)      

        step_size_label = QLabel('Step size:')
        self.step_size_input = QDoubleSpinBox()
        self.step_size_input.setRange(0.001, 4500.0)  
        self.step_size_input.setSingleStep(1)    
        self.step_size_input.setValue(1.0)
        self.step_size_input.setSuffix(' cm-1') 
        self.step_size_input.setDecimals(1) 
        self.step_size_input.setMinimumWidth(100)        

        bounds_layout.addWidget(lower_bound_label)
        bounds_layout.addWidget(self.lower_bound_input)
        bounds_layout.addSpacing(10)

        bounds_layout.addWidget(upper_bound_label)
        bounds_layout.addWidget(self.upper_bound_input)
        bounds_layout.addSpacing(10)

        bounds_layout.addWidget(step_size_label)
        bounds_layout.addWidget(self.step_size_input)
        bounds_layout.addStretch(1)
        
        layout.addLayout(bounds_layout)
        layout.addSpacing(5)

        #FWHM of the Gaussian convolution and scaling factor
        peak_options = QHBoxLayout()

        #FWHM Input
        fwhm_label = QLabel('FWHM:')
        self.fwhm_input = QDoubleSpinBox()
        self.fwhm_input.setRange(0.1, 50.0)  #Set an appropriate range
        self.fwhm_input.setValue(8.00)
        self.fwhm_input.setSingleStep(0.10)
        self.fwhm_input.setMinimumWidth(80)        

        peak_options.addWidget(fwhm_label)
        peak_options.addWidget(self.fwhm_input)
        peak_options.addSpacing(10)

        #Harmonic Scaling Factor Input
        vib_scl_label = QLabel('Scaling Factor:')
        self.vib_scl_input = QDoubleSpinBox()
        self.vib_scl_input.setRange(0.01, 2.0) 
        self.vib_scl_input.setSingleStep(0.01)
        self.vib_scl_input.setValue(1.0)
        self.vib_scl_input.setMinimumWidth(80)        
        
        peak_options.addWidget(vib_scl_label)
        peak_options.addWidget(self.vib_scl_input)
        peak_options.addStretch(1)

        layout.addLayout(peak_options)
        layout.addSpacing(5)

        #Normalization Checkbox
        self.normalization_checkbox = QCheckBox("Normalize Spectra?")
        layout.addWidget(self.normalization_checkbox)

        #Plotting Checkbox
        self.plotting_checkbox = QCheckBox("Plot Spectra externally?")
        layout.addWidget(self.plotting_checkbox)
        layout.addSpacing(10)

        #Run Button
        run_button = QPushButton('Extract/Plot IR Spectra')
        run_button.clicked.connect(self.run_IR_extract)
        layout.addWidget(run_button)

        #Add some vertical spacing
        layout.addStretch(1)

        #Set the layout
        layout.setContentsMargins(30, 30, 30, 30)
        self.setLayout(layout)

    def browse_directory(self):
        directory_path = QFileDialog.getOpenFileName(self, "Select .out file", "", "OUT Files (*.out)")
        if directory_path[0]:
            self.directory_input.setText(directory_path[0])

    def run_IR_extract(self):
        #Get the selected directory or CSV file
        directory = self.directory_input.text()
        fwhm = self.fwhm_input.value()
        lower_bound = self.lower_bound_input.value()
        upper_bound = self.upper_bound_input.value()
        step_size = self.step_size_input.value()
        vib_scl = self.vib_scl_input.value()
        normalization = self.normalization_checkbox.isChecked() #true if checked, false if not
        plotting = self.plotting_checkbox.isChecked() #true if checked, false if not

        #Check if the values meet your specific conditions
        if fwhm <= 0:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} FWHM must be greater than 0.\n')
            return
        if lower_bound < 0:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The lower-bound of the frequency range to the exported cannot be negative.')
            return
        if upper_bound <= lower_bound:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The upper bound must be greater than the lower-bound of the frequency range to the exported.')
            return
        
        if step_size <= 0:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The step size of for the frequency range interpolation must be greater than 0.')
            return
        
        if vib_scl <= 0:  #Assuming vib_scl should be between 0 and 2
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The harmonic scaling factor must be greater than zero.')
            return
        
        if not os.path.isdir(directory):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {directory} is not a valid directory. Please check that the file you have specified is correct, then try again.')
            return
        
        #if all checks pass, proceed to extract the IR spectra from each .out file
        extract_IR_spectra(directory, fwhm, lower_bound, upper_bound, step_size, normalization, plotting, vib_scl)

class Extract_TDDFT_VG_SpectraTab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)

        #Directory Label
        directory_layout = QHBoxLayout()  

        directory_label = QLabel('Directory:')
        directory_layout.addWidget(directory_label)
        
        self.directory_input = QLineEdit()
        directory_layout.addWidget(self.directory_input)

        browse_button = QPushButton('Browse')
        browse_button.clicked.connect(self.browse_directory)
        directory_layout.addWidget(browse_button)

        layout.addLayout(directory_layout)
        layout.addSpacing(5)

        #Unit Specifications
        unit_layout = QHBoxLayout()  

        #Input unit combobox
        input_unit_label = QLabel('Input unit:')
        self.input_unit_combo = QComboBox()
        self.input_unit_combo.addItems(['nm', 'cm**-1', 'eV'])
        self.input_unit_combo.setCurrentText('wavenumber')  #Default to wavenumber
        self.input_unit_combo.setMinimumWidth(30)  
        
        unit_layout.addWidget(input_unit_label)
        unit_layout.addWidget(self.input_unit_combo)
        unit_layout.addSpacing(25)

        #Output Unit ComboBox
        output_unit_label = QLabel('Output unit:')
        self.output_unit_combo = QComboBox()
        self.output_unit_combo.addItems(['nm', 'cm**-1', 'eV'])
        self.output_unit_combo.setCurrentText('wavenumber')  #Default to wavenumber
        self.output_unit_combo.currentTextChanged.connect(self.update_output_values) #Connect output unit to the auto-update method
        self.output_unit_combo.setMinimumWidth(30)  

        unit_layout.addWidget(output_unit_label)
        unit_layout.addWidget(self.output_unit_combo)
        unit_layout.addStretch(1)

        layout.addLayout(unit_layout)
        layout.addSpacing(5)

        #Output Value Range
        output_value_layout = QHBoxLayout()

        #Lower bound
        output_lower_label = QLabel('Export from:')
        self.output_lower_input = QDoubleSpinBox()
        self.output_lower_input.setRange(0,1000000)
        self.output_lower_input.setValue(150.00)
        self.output_lower_input.setMinimumWidth(90)  
        self.output_lower_input.setSuffix(' nm') #default to wavelength

        #Upper bound
        output_upper_label = QLabel('to:')
        self.output_upper_input = QDoubleSpinBox()
        self.output_upper_input.setRange(0,1000000)
        self.output_upper_input.setValue(400.00)
        self.output_upper_input.setMinimumWidth(90)  
        self.output_upper_input.setSuffix(' nm') #default to wavelength

        #Step size
        output_spacing_label = QLabel('Step size:')
        self.output_spacing_input = QDoubleSpinBox()
        self.output_spacing_input.setRange(0,100000)
        self.output_spacing_input.setValue(1.00)
        self.output_spacing_input.setMinimumWidth(90)  
        self.output_spacing_input.setSuffix(' nm') #default to wavelength

        #Bound QHBox layout        
        output_value_layout.addWidget(output_lower_label)
        output_value_layout.addWidget(self.output_lower_input)
        output_value_layout.addSpacing(10)

        output_value_layout.addWidget(output_upper_label)
        output_value_layout.addWidget(self.output_upper_input)
        output_value_layout.addSpacing(10)

        output_value_layout.addWidget(output_spacing_label)
        output_value_layout.addWidget(self.output_spacing_input)
        output_value_layout.addStretch(1)

        layout.addLayout(output_value_layout)
        layout.addSpacing(5)

        #Shift spectra layout
        shift_layout = QHBoxLayout()

        #Shift Input
        shift_label = QLabel('Shift spectra by:')
        self.shift_input = QDoubleSpinBox()
        self.shift_input.setRange(-200000,200000)
        self.shift_input.setValue(-0.45)
        self.shift_input.setMinimumWidth(50)  

        #Shift Unit ComboBox
        self.shift_unit_combo = QComboBox()
        self.shift_unit_combo.addItems(['nm', 'cm**-1', 'eV'])
        self.shift_unit_combo.setCurrentText('eV')  #Default to eV
        self.shift_unit_combo.currentTextChanged.connect(self.update_shift_value) #Connect shift unit to the auto-update method
        self.shift_unit_combo.setMinimumWidth(30) 

        shift_layout.addWidget(shift_label)
        shift_layout.addWidget(self.shift_input)
        shift_layout.addWidget(self.shift_unit_combo)
        shift_layout.addStretch(1)

        layout.addLayout(shift_layout)
        layout.addSpacing(5)

        #Basename input
        basename_layout = QHBoxLayout()
        output_file_label = QLabel('Basename:')
        self.output_file_input = QLineEdit()
        self.output_file_input.setPlaceholderText('Enter the basename of your .spectrum and .spectrum.rootN files.')
        
        basename_layout.addWidget(output_file_label)
        basename_layout.addWidget(self.output_file_input)

        layout.addLayout(basename_layout)
        layout.addSpacing(5)

        #Normalization checkbox
        self.normalization_checkbox = QCheckBox('Normalize Output?')
        self.normalization_checkbox.setChecked(True)  #Default checked
        layout.addWidget(self.normalization_checkbox)

        #Plotting checkbox
        self.plot_checkbox = QCheckBox('Plot spectra externally?')
        layout.addWidget(self.plot_checkbox)
        layout.addSpacing(10)

        #Run Button
        run_button = QPushButton('Extract UV-Vis Spectra (ESD)')
        run_button.clicked.connect(self.run_extract_ESD)
        layout.addWidget(run_button)

        layout.addStretch(1)
        layout.setContentsMargins(30, 30, 30, 30)
        self.setLayout(layout)

    def browse_directory(self):
        directory_path = QFileDialog.getExistingDirectory(self, "Select Directory")

        if directory_path:
            self.directory_input.setText(directory_path)
    
    #update output interpolation range depending on selection of the output unit
    def update_output_values(self, unit):
        if unit == 'nm':
            self.output_lower_input.setValue(150)
            self.output_upper_input.setValue(400)
            self.output_spacing_input.setValue(1)
        elif unit == 'cm**-1':
            self.output_lower_input.setValue(25000)
            self.output_upper_input.setValue(67000)
            self.output_spacing_input.setValue(50)
        elif unit == 'eV':
            self.output_lower_input.setValue(3.1)
            self.output_upper_input.setValue(8.3)
            self.output_spacing_input.setValue(0.05)
        
        #Update the suffix of the spin boxes based on the selected unit
        self.output_lower_input.setSuffix(f' {unit}')
        self.output_upper_input.setSuffix(f' {unit}')
        self.output_spacing_input.setSuffix(f' {unit}')

    #update correction factor depending on selection of the shift unit
    def update_shift_value(self, unit):
        if unit == 'nm':
            self.shift_input.setValue(20.0)
        elif unit == 'cm**-1':
            self.shift_input.setValue(-450.0)
        elif unit == 'eV':
            self.shift_input.setValue(-0.25)

    def run_extract_ESD(self):
        
        #get parameters from GUI interface and assign them to variables
        directory = self.directory_input.text()
        input_unit = self.input_unit_combo.currentText()
        output_unit = self.output_unit_combo.currentText()
        shift = self.shift_input.value()
        shift_unit = self.shift_unit_combo.currentText()
        output_lowervalue = self.output_lower_input.value()
        output_uppervalue = self.output_upper_input.value()
        output_spacing = self.output_spacing_input.value()
        normalize_output = self.normalization_checkbox.isChecked()
        output_file_basename = self.output_file_input.text()
        plotting = self.plot_checkbox.isChecked()

        #check that data provided in the GUI is correctly set up
        if not os.path.isdir(directory):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {directory} is not a valid directory. Please check that the directory you have specified is correct, then try again.')
            return

        if output_lowervalue < 0:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The lower bound of the output window cannot be negative!')
            return

        if output_uppervalue <= output_lowervalue:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The lower bound of the output window must be greater than the upper bound!')
            return
        
        if output_spacing < 0:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The step size for the output interpolation cannot be negative!')
            return
        
        #check for .spectrum.root and .spectrum files, and choose the appropriate analysis depending on what is present
        contains_spectrum_root = any('.spectrum.root' in file for file in os.listdir(directory))
        ends_with_spectrum = any(file.endswith('.spectrum') for file in os.listdir(directory))

        #If directory contains both .spectrum and .spectrum.root files or just .spectrum.root files, the full analysis can be done using only .spectrum.root files:
        if (contains_spectrum_root and ends_with_spectrum):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Both .spectrum.root and .spectrum files were found in the provided directory. All spectra will be processed from the root files. Please ensure that all .spectrum.root files are present, otherwise the analysis will be incomplete.')
            QApplication.processEvents()
            extract_ESD_spectrum_root_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting)

        #If directory contains .spectrum files, but not .spectrum.root files:
        elif ends_with_spectrum and not contains_spectrum_root:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Only .spectrum files were found in the provided directory. Only total spectra will be processed. To see the contributions from each excited state, please add the respective .spectrum.root files to the specified directory.')
            QApplication.processEvents()
            extract_ESD_spectrum_files(directory, input_unit, shift_unit, output_unit, shift, output_lowervalue, output_uppervalue, output_spacing, output_file_basename, normalize_output, plotting)

        #If directory does not contain any .spectrum files or .spectrum.root, inform the user:
        else: 
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No .spectrum or .spectrum.root files were found in {directory}. Please ensure that you have specified the correct directory.')
            return        

class BWCCS_Tab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)
        
        #DFT Thermochem input
        DFT_thermochem_layout = QHBoxLayout()  

        thermochem_label = QLabel('DFT Thermochem .csv:')
        self.thermochem_input = QLineEdit()
        self.thermochem_input.setPlaceholderText('Directory + name of .csv file containing DFT electronic energies and thermochemistry.') 
        
        self.browse_button_thermochem = QPushButton('Browse')
        self.browse_button_thermochem.clicked.connect(self.browse_thermochem_file)
        
        DFT_thermochem_layout.addWidget(thermochem_label)
        DFT_thermochem_layout.addWidget(self.thermochem_input)
        DFT_thermochem_layout.addWidget(self.browse_button_thermochem)

        layout.addLayout(DFT_thermochem_layout)
        layout.addSpacing(5)

        #DLPNO-CCSD(T) input - optional
        CCSDT_layout = QHBoxLayout()  

        CCSDT_label = QLabel('DLPNO-CCSD(T) .csv:')
        self.CCSDT_input = QLineEdit()
        self.CCSDT_input.setPlaceholderText('(Optional) Directory + name of .csv file containing DLPNO-CCSD(T) energies.') 

        self.browse_button_CCSDT = QPushButton('Browse')
        self.browse_button_CCSDT.clicked.connect(self.browse_CCSDT_file)

        CCSDT_layout.addWidget(CCSDT_label)
        CCSDT_layout.addWidget(self.CCSDT_input)
        CCSDT_layout.addWidget(self.browse_button_CCSDT)

        layout.addLayout(CCSDT_layout)
        layout.addSpacing(5)        

        #CCS input from MobCal-MPI 2.0
        CCS_layout = QHBoxLayout()  

        CCS_label = QLabel('CCS .csv:')
        self.CCS_input = QLineEdit()
        self.CCS_input.setPlaceholderText('Directory + name of .csv file containing CCSs computed from MobCal-MPI 2.0') 

        self.browse_button_CCS = QPushButton('Browse')
        self.browse_button_CCS.clicked.connect(self.browse_CCS_file)

        CCS_layout.addWidget(CCS_label)
        CCS_layout.addWidget(self.CCS_input)
        CCS_layout.addWidget(self.browse_button_CCS)

        layout.addLayout(CCS_layout)
        layout.addSpacing(5)  

        #Output filename 
        output_layout = QHBoxLayout()  

        output_label = QLabel('Output .xlsx:')
        self.output_input = QLineEdit()
        self.output_input.setPlaceholderText('Directory + name of a .xlsx to write the data to') 

        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_input)

        layout.addLayout(output_layout)
        layout.addSpacing(5)  

        #Temperature input for calculating populations
        temp_layout = QHBoxLayout()

        #Temp input
        temp_label = QLabel('Temp (K):')
        self.temp_input = QDoubleSpinBox()
        self.temp_input.setMinimum(0.01)
        self.temp_input.setMaximum(10000.00)
        self.temp_input.setValue(298.15)
        self.temp_input.setSingleStep(5)
        self.temp_input.setMinimumWidth(80)        

        temp_layout.addWidget(temp_label)
        temp_layout.addWidget(self.temp_input)
        temp_layout.addStretch(1)

        layout.addLayout(temp_layout)
        layout.addSpacing(10)

        #Run Button
        run_button = QPushButton('Run BW-CCS Analysis')
        run_button.clicked.connect(self.run_BW_CCS)
        layout.addWidget(run_button)

        layout.addStretch(1)

        layout.setContentsMargins(30, 30, 30, 30)  
        self.setLayout(layout)

    def browse_thermochem_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select .csv file', '', 'CSV Files (*.csv)')

        #Check if a file path was selected
        if file_path:
            self.thermochem_input.setText(file_path)
    
    def browse_CCSDT_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select .csv file', '', 'CSV Files (*.csv)')

        #Check if a file path was selected
        if file_path:
            self.CCSDT_input.setText(file_path)

    def browse_CCS_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select .csv file', '', 'CSV Files (*.csv)')

        #Check if a file path was selected
        if file_path:
            self.CCS_input.setText(file_path)

    def run_BW_CCS(self):
        
        #get parameters from GUI interface and assign them to variables
        thermochem_file = self.thermochem_input.text()
        CCSDT_file = self.CCSDT_input.text()
        CCS_file = self.CCS_input.text()
        output_file = self.output_input.text()
        temp = self.temp_input.value()

        #check that data provided in the GUI is correctly set up
        if not os.path.isfile(thermochem_file):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {thermochem_file} is not a valid file. Please check that the directory you have specified is correct, then try again.')
            return

        elif not thermochem_file.lower().endswith('.csv'):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The thermochem file must be a directory path + filename that ends with .csv. This should also be the output from the Thermochem Analyzer function in this GUI.')

        elif not os.path.isfile(CCS_file):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {CCS_file} is not a valid file. Please check that the directory you have specified is correct, then try again.')
            return

        elif not CCS_file.lower().endswith('.csv'):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The CCS file must be a directory path + filename that ends with .csv. This should also be the output from the MobCal-MPI 2.0 GUI moutAnalyzer function.')

        elif not os.path.isdir(os.path.dirname(output_file)):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {os.path.dirname(output_file)} is not a directory to write the output .xlsx file to. Please create it or specifcy the file in a directory that exists, then try again.')
            return

        elif not output_file.lower().endswith('.xlsx'):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The output file must be a directory path + filename that ends with .xlsx.')
            return

        elif not temp > 0:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The temperature specified must be greater than 0.')
            return           
        
        elif not os.path.isfile(CCSDT_file):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {CCS_file} is not a valid file. Please check that the directory you have specified is correct, then try again.')
            return            

        elif not CCSDT_file.lower().endswith('.csv'):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The DLPNO-CCSD(T) file must be a directory path + filename that ends with .csv. This should also be the output from the CCSDT Analyzer function in this GUI.')

        #Check if a DLPNO-CCSD(T) file was given. If not, assign it a value of None
        elif CCSDT_file.strip() == '':
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No DLPNO-CCSD(T) file was provided. The BW-CCS will be calculated using DFT electronic eneriges.')
            CCSDT_file = None

        #if all checks pass, proceed to calcualte the BW CCS
        else:
            BW_CCS_Analysis(thermochem_file, CCSDT_file, CCS_file, output_file, temp)

class LED_Analysis_Tab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)
        
        #DLPNO-CCSD(T) of Fragment 1 .out file inputs
        frag1_layout = QHBoxLayout()  

        frag1_label = QLabel('Fragment-1 .out:')
        self.frag1_input = QLineEdit()
        self.frag1_input.setPlaceholderText('dir + filename of the DLPNO-CCSD(T) .out for Fragment 1 from the LED job.') 
        
        self.browse_button_frag1 = QPushButton('Browse')
        self.browse_button_frag1.clicked.connect(self.browse_frag1_file)
        
        frag1_layout.addWidget(frag1_label)
        frag1_layout.addWidget(self.frag1_input)
        frag1_layout.addWidget(self.browse_button_frag1)

        layout.addLayout(frag1_layout)
        layout.addSpacing(5)

        #DLPNO-CCSD(T) of Fragment 2 .out file inputs
        frag2_layout = QHBoxLayout()  

        frag2_label = QLabel('Fragment-2 .out:')
        self.frag2_input = QLineEdit()
        self.frag2_input.setPlaceholderText('dir + filename of the DLPNO-CCSD(T) .out for Fragment 2 from the LED job.') 
        
        self.browse_button_frag2 = QPushButton('Browse')
        self.browse_button_frag2.clicked.connect(self.browse_frag2_file)
        
        frag2_layout.addWidget(frag2_label)
        frag2_layout.addWidget(self.frag2_input)
        frag2_layout.addWidget(self.browse_button_frag2)

        layout.addLayout(frag2_layout)
        layout.addSpacing(5)    

        #Parent-LED .out file input
        Parent_LED_layout = QHBoxLayout()  

        Parent_LED_label = QLabel('Parent LED .out:')
        self.Parent_LED_input = QLineEdit()
        self.Parent_LED_input.setPlaceholderText('dir + filename of the .out file for the LED job of the parent molecule.') 

        self.browse_button_Parent_LED = QPushButton('Browse')
        self.browse_button_Parent_LED.clicked.connect(self.browse_Parent_LED_file)

        Parent_LED_layout.addWidget(Parent_LED_label)
        Parent_LED_layout.addSpacing(5)
        Parent_LED_layout.addWidget(self.Parent_LED_input)
        Parent_LED_layout.addWidget(self.browse_button_Parent_LED)

        layout.addLayout(Parent_LED_layout)
        layout.addSpacing(10)    

        #Run Button
        run_button = QPushButton('Run LED Analysis')
        run_button.clicked.connect(self.run_LED)
        layout.addWidget(run_button)

        layout.addStretch(1)

        layout.setContentsMargins(30, 30, 30, 30)  
        self.setLayout(layout)

    def browse_frag1_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select .csv file', '', 'CSV Files (*.csv)')

        #Check if a file path was selected
        if file_path:
            self.frag1_input.setText(file_path)
    
    def browse_frag2_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select .csv file', '', 'CSV Files (*.csv)')

        #Check if a file path was selected
        if file_path:
            self.frag2_input.setText(file_path)

    def browse_Parent_LED_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, 'Select .csv file', '', 'CSV Files (*.csv)')

        #Check if a file path was selected
        if file_path:
            self.Parent_LED_input.setText(file_path)

    def run_LED(self):
        
        #get parameters from GUI interface and assign them to variables
        frag1_file = self.frag1_input()
        frag2_file = self.frag2_input()
        Parent_LED_file = self.Parent_LED_input()

        #check that data provided in the GUI is correctly set up
        if not os.path.isfile(frag1_file):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {frag1_file} is not a valid file. Please check that the directory you have specified is correct, then try again.')

        elif not os.path.isfile(frag2_file):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {frag2_file} is not a valid file. Please check that the directory you have specified is correct, then try again.')

        elif not os.path.isfile(Parent_LED_file):
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {Parent_LED_file} is not a valid file. Please check that the directory you have specified is correct, then try again.')

        #if all checks pass, proceed to extract spectra
        else:
            LED_Analysis(thermochem_file, CCSDT_file, CCS_file, output_file = 'LED_Analysis.csv')

class NEB_Analysis_Tab(QWidget):
    def __init__(self, text_redirector, parent=None):
        super().__init__(parent)
        self.text_redirector = text_redirector
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout(self)
        
        message = QLabel('This feature will be added in the next update of the ORCA Analysis GUI. This tab is just a placeholder.')
        
        layout.addWidget(message)
        layout.addStretch(1)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_win = ORCAAnalysisSuite()
    main_win.show()
    sys.exit(app.exec())
