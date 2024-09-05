import os, re, time
import numpy as np
import pandas as pd
from datetime import datetime
from Python.constants_and_conversions import c

#supress a specific warning about a pandas update that isn't going to impact the code
import warnings

# Suppress all FutureWarnings
warnings.filterwarnings("ignore", category=FutureWarning)

from PyQt6.QtWidgets import QApplication
import numpy as np

def ORCA_Thermochem_Calculator(directory, T = 298.15, p = 101325., vib_scl = 1., sort_by = 'F'):

    def read_molecule_data(file, vib_scl):
        '''Reads data from the specified file in the given directory.'''
        
        with open(file, 'r') as opf:
            data = opf.readlines()

        #Initialize variables and arrays
        charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI = None, None, None, None, None, None, None
        total_dipole, dipole_x, dipole_y, dipole_z, polariz, dipole_ax = None, None, None, None, None, None
        vibs_array, ZPE, Total_ZPE = None, None, None
        n_imag = 0
        imag_freqs = None
        spin_contam = None

        #flag to check for normal termination
        normal_term = False

        #Initialize flag to check if thermochem was requested in the method line
        thermochem_flag = False  

        #First, we find the section of the out file that contains the data imported from the .inp - this will contain the freq flag that determined whether thermochem is calcualted
        for line in data:
            #Starting from the top of the file, look for freq in each line. If it is found, break the loop
            if re.search(r'\bfreq\b', line, re.IGNORECASE):
                thermochem_flag = True
                break
            
            #There is no point in searching lines after the termination of the input block, so we're defining a flag that will break this loop when set to true
            elif '*END OF INPUT*' in line:
                break

        #get info from files depending on if thermochem was calculated or not

        for line in data:
            
            #charge
            if line.startswith(' Total Charge'):
                charge = float(line.split()[-1])
            
            #multiplicity
            elif line.startswith(' Multiplicity'):
                multi = float(line.split()[-1])

            #spin contamination
            elif line.startswith('Deviation                       :'):
                spin_contam = float(line.split()[-1])
            
            #final electronic energy
            elif line.startswith('FINAL SINGLE POINT ENERGY'):
                Eelec = float(line.split()[-1])

            #rotational constants
            elif line.startswith('Rotational constants in cm-1'):
                RotA = float(line.split()[4]) * 100.  #in m**-1
                RotB = float(line.split()[5]) * 100.
                RotC = float(line.split()[6]) * 100.
                RotABC = np.array([RotA, RotB, RotC])  #in m**-1

            #point group
            elif line.startswith('Point Group:'):
                sigma_OR = int(line.split()[-1])  #integer

            #mass
            elif line.startswith('Total Mass          ...'):
                mass = float(line.split()[-2])
                m_SI = mass * c['amu2kg']

            #dipole moment
            elif line.startswith('Magnitude (Debye)'):
                total_dipole = float(line.split()[-1])  #in Debye

            #dipole along x, y, and z axes
            elif line.startswith('Total Dipole Moment    :'):
                
                #given in a.u, then converted to Debye
                dipoles = line.split()[4:] #XYZ components are given after the 4th entry because the dipoles are prefixed by "Total Dipole Moment    :", which is split by whitespaces

                #conversion factor - 1 Debye = 3.33564E-30 C*m, and 1 a.u of the electricdipole moment is 8.47835326E-30 Cm
                au_to_Debye = 8.47835326E-30 / 3.33564E-30

                try:
                    dipole_x, dipole_y, dipole_z = [float(dipole) * au_to_Debye for dipole in dipoles]

                except (IndexError, TypeError, ValueError) as e:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} XYZ components of the dipole from {os.path.basename(file)} could not be read properly. Error: {e}')
                    pass

            #polarizability
            elif line.startswith('Isotropic polarizability :'):
                polariz = float(line.split()[-1]) * 1.4818e-31  #in m^3

            #dipole moment along the rotational axes
            elif line.startswith('x,y,z [Debye]'):
                mu_abc = np.array([float(s) for s in line.split()[-3:]])
                
                #max dipole along each of the rotational axes
                dipole_ax = ['A', 'B', 'C'][np.argmax(np.abs(mu_abc))]
            
            elif "****ORCA TERMINATED NORMALLY****" in line:
                normal_term = True
        
        #only look for vibrations if 'freq' was specificed in the method line
        if thermochem_flag:
            
            #Search for imaginary freqs in the last printing of vib frequencies (in case multiple calculations of Hessian is requested by the user during OptTS )
            #get the index of the line where the last vib freq block is printed
            last_vib_header_index = None
            last_vib_footer_index = None
            for i, line in enumerate(data):
                
                #find start of the LAST vib header section
                if line.strip() == '-----------------------' and i + 1 < len(data) and data[i + 1].strip() == 'VIBRATIONAL FREQUENCIES':
                    last_vib_header_index = i + 4  #The line after the header and scaling factor line

                if line.strip() == '------------' and i + 1 < len(data) and data[i + 1].strip() == 'NORMAL MODES':
                    last_vib_footer_index = i  #The lines before the 'NORMAL MODES' block
            
            if last_vib_header_index is not None and last_vib_footer_index is not None:
                vibs_array = []
                for line in data[last_vib_header_index:last_vib_footer_index]:             
                   
                    #get vib freqs
                    if 'cm**-1' in line:                       
                        if '***imaginary mode***' not in line:
                            vib_str = line.split()[1]  # Extract the vib freq from strings of the format:    1:         0.00 cm**-1
                            try:
                                vib = float(line.split()[1]) * vib_scl
                            except (TypeError, ValueError):
                                print(f'Error: Unable to convert {vib_str} to a float from {os.path.basename(file)}.')
                                continue
                            
                            if vib > 0:  # Exclude zero and imaginary modes
                                vibs_array.append(vib)  # in cm**-1
      
                        #SEARCH FOR IMAGINARY MODES
                        if '***imaginary mode***' in line:
                            if imag_freqs is None:
                                imag_freqs = []
                            
                            vib_str = line.split()[1]  # Extract the vib freq from strings of the format:       6:      -502.32 cm**-1 ***imaginary mode***
                            try:
                                vib = float(line.split()[1]) #no scaling applied to imag mobe
                                if vib < 0:
                                    imag_freqs.append(vib)
                                    n_imag += 1

                            except ValueError:
                                print(f'Error: Unable to convert {vib_str} to a float from {os.path.basename(file)}.')
                                continue

        #After all data is extracted from the .out file and thermochem was requested (and was completed!), calculate he ZPE from the now non-None vibs array
        if vibs_array is not None:
            vibs_array = np.array(vibs_array)
            ZPE = 0.5 * c['h_SI'] * c['c_SI'] * 100. * np.sum(vibs_array) * c['J2Eh']  #Hartree
            Total_ZPE = Eelec + ZPE

        return charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, total_dipole, dipole_x, dipole_y, dipole_z, polariz, dipole_ax, vibs_array, ZPE, Total_ZPE, n_imag, imag_freqs, normal_term, thermochem_flag, spin_contam
        
    def calc_partition_function(filename, m_SI, RotABC, sigma_OR, vibs_array, multi, T, p):
        '''Calculates the partition functions (trans, rot, vib, elec) for the molecule at a given p, T. Methodology follows that of “Molecular Thermodynamics” by McQuarrie and Simon (1999)'''

        #Set np to raise exceptions for divide by zero
        np.seterr(divide='raise')

        #translational partition function
        q_trans = (2. * np.pi * m_SI * c['kB'] * T / c['h_SI'] ** 2) ** 1.5 * c['kB'] * T / p

        #rotational partitional function - three cases: non-linear, linear, and atomic
        try:
            
            #non-linear polyatomic
            q_rot = 1. / sigma_OR * (c['kB'] * T / (c['h_SI'] * c['c_SI'])) ** 1.5 * (np.pi / (np.prod(RotABC))) ** 0.5
            rot_type = 'polyatomic'

        #if the molecule is a diatomic, linear polyatomic, or a single atom, a FloatingPointError error will be encountered. In this case, we can evaualte the rotational parititon function using that for a diatomic
        except FloatingPointError:
            
            #any linear molecule
            try: 
                q_rot = 1. / sigma_OR * (c['kB'] * T / (c['h_SI'] * c['c_SI'] * np.max(RotABC))) 
                rot_type = 'linear'

                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {filename} is a linear molecule. Rotational contributions to energy and entropy will be adjusted accordingly.')
                QApplication.processEvents()

            #if still a FloatingPointError, then the calculation is on a single atom, in which case, the q_rot is 1.
            except FloatingPointError:

                q_rot = 1.
                rot_type = 'atom'
                
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {filename} is a single atom. Rotational contributions to energy and entropy will be adjusted accordingly.')
                QApplication.processEvents()

        #Set np back to default warning
        np.seterr(divide='warn')

        #vibrational partition functions - 2 cases: atomic and non-atomic

        #non-atomic entities
        if rot_type != 'atom':
            u = c['h_SI'] * c['c_SI'] * vibs_array * 100. / (c['kB'] * T) #factor of 100 is for conversion of cm-1 to m-1
            q_vib = np.prod(1. / (1. - np.exp(-u)))

        #atoms dont have vibrations, so their parition function is unity. 
        else:
            q_vib = 1.
        
        #electronic partition function - equates to the multiplicity of the ion for species in their ground state
        q_elec = multi
        
        return [q_trans, q_rot, q_vib, q_elec, rot_type]

    def calc_thermochemistry(vibs_array, q_trans, q_rot, q_vib, q_elec, ZPE, mol_Eelec, T, p, rot_type):
        '''Computes standard thermochemical functions at given p, T.
        Methodology follows that of “Molecular Thermodynamics” by McQuarrie and Simon (1999), which was summarized by Ochterski in 2000 (https://gaussian.com/wp-content/uploads/dl/thermo.pdf), and uses slightly different method for vibrational contibutions compared to the methods used in ORCA.
        Consequently, the thermochemical corrections computed here will be slightly different than the printouts in the ORCA .out file. Both are valid, although these values generated by this code more generally applicable.'''

        #Translational contributions
        S_trans = c['kB_Eh'] * (np.log(q_trans) + 1. + (3. / 2.)) #translational contribution to entropy
        E_trans = (3. / 2.) * c['kB_Eh'] * T #translational contribution to energy
        
        #Rotational contributions; 3 cases for rotation: single atom, linear, or non-linear polyatomic
        if rot_type == 'atom':
            #atoms have no rotational constants, so contributions from rotation to entropy and energy are both zero
            S_rot = 0
            E_rot = 0
        
        elif rot_type == 'linear':
            S_rot = c['kB_Eh'] * (np.log(q_rot) + 1.) #rotational contribution to entropy for linear molecule
            E_rot = c['kB_Eh'] * T #rotational contribution to energy for linear molecule

        else:
            S_rot = c['kB_Eh'] * (np.log(q_rot) + (3. / 2.)) #rotational contribution to entropy for a non-linear polyatomic
            E_rot = 3. / 2. * c['kB_Eh'] * T #rotational contribution to energy for a non-linear polyatomic
        
        #Vibrational contributions: atomic vs. non-atomic
        if rot_type != 'atom':
            
            theta_v = c['h_SI'] * c['c_SI'] * vibs_array * 100. / c['kB']
            theta_v_oT = theta_v / T
            
            S_vib = c['kB_Eh'] * np.sum((theta_v_oT / (np.exp(theta_v_oT) - 1.)) - np.log(1. - np.exp(-theta_v_oT))) #vibrational contribution to entropy - slightly different than the default method used by ORCA
            E_vib = c['kB_Eh'] * np.sum(theta_v / (np.exp(theta_v_oT) - 1.)) #vibrational contribution to energy
        
        else:
            S_vib = 0
            E_vib = 0
            ZPE = 0 #needed for calculation of total energy later
        
        #Electronic contributions
        S_elec = c['kB_Eh'] * np.log(q_elec) #electronic contribution to entropy
        E_elec = 0. #always zero for systems in their electronic ground state

        #Contributions to total energy
        Total_S = S_trans + S_rot + S_vib + S_elec
        
        Ecorr = E_trans + E_rot + E_vib + E_elec #total energy
        Total_E = Ecorr + ZPE + mol_Eelec #note mol_Eelec is the electronic energy of the analyte

        Hcorr = (Ecorr + c['kB_Eh'] * T) + ZPE
        Total_H = Hcorr + mol_Eelec

        Gcorr = ((Hcorr - ZPE) - T * Total_S) + ZPE
        Total_G = Gcorr + mol_Eelec

        return Ecorr, Total_E, Hcorr, Total_H, Gcorr, Total_G, Total_S

    '''Main code operation'''
    #get list of filenames, and check if the directory does not contain any .out files. If using pseudopotentials, an _atom will be added to the filename - we don't want those files!
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x.lower()]

    if len(filenames) == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} There are no .out files in the provided directory.')
        return

    start = time.time()

    properties = [
        'Filename',
        'Imaginary Frequencies',
        'Electronic energy',
        'ZPE correction',
        'Thermal correction',
        'Enthalpy correction',
        'Gibbs Correction',
        'Total ZPE',
        'Total Thermal Energy',
        'Total Enthalpy',
        'Total Entropy (T*S)',
        'Total Gibbs Energy',
        'Relative Gibbs Energy',
        'Rot. const. A (cm**-1)',
        'Rot. const. B (cm**-1)',
        'Rot. const. C (cm**-1)',
        'Rot. Symm. Number',
        'Dipole X',
        'Dipole Y',
        'Dipole Z',
        'Total Dipole',
        'Isotropic Polarizability',
        'Trans. Part. func.',
        'Rot. Part. func.',
        'Vib. Part. function',
        'Elec. Part. func.',
        'Spin contam: Deviation from S*(S+1)'
    ]

    #initialize pd FataFrame w/ columns assoicated w/ thermo properties
    df = pd.DataFrame(columns=properties)

    #initialize lists for storing data
    Gibbs_list = []
    thermochem_not_requested = []
    master_imag_freq_list = []
    
    abnormal_term_list = []
    abnormal_term_wVibs_list = []

    for filename in filenames:

        #get data from input file
        charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, total_dipole, dipole_x, dipole_y, dipole_z, polariz, dipole_ax, vibs_array, ZPE, Total_ZPE, n_imag, imag_freqs, normal_term, thermochem_flag, spin_contam = read_molecule_data(os.path.join(directory, filename), vib_scl)

        #evaluate rotational contants if numbers were extracted. otherwise, return N/A
        try:
            RotA = RotABC[0]/100
            RotB = RotABC[1]/100
            RotC = RotABC[2]/100
        
        except (IndexError, TypeError):
            RotA = 'N/A'
            RotB = 'N/A'
            RotC = 'N/A'

        #if the job terminated normally
        if normal_term:
            
            #calculate parition functions from input data if thermochem was requested
            if thermochem_flag:
                
                #try to calculate thermochemistry
                try:
                    q_trans, q_rot, q_vib, q_elec, rot_type = calc_partition_function(filename, m_SI, RotABC, sigma_OR, vibs_array, multi, T, p)

                    #calculate thermochemical corrections
                    Ecorr, Total_E, Hcorr, Total_H, Gcorr, Total_G, Total_S = calc_thermochemistry(vibs_array, q_trans, q_rot, q_vib, q_elec, ZPE, Eelec, T, p, rot_type)
                    Gibbs_list.append(Total_G)
                
                #if thermochem corrections cannot be calculated (such as because of a incomplete job, write placeholders)
                except Exception as e:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error during the calculation of thermochemistry from {filename}: {e}. Please report this issue along with your .out file to the issues section of the Github repo.')

                #if a file contians imaginary frequencies, extract their values and write them to a list for printing laster on. This is useful for determining whether a DFT job needs to be resubmitted.
                if n_imag > 0:
                
                    try:
                        master_imag_freq_list.append([filename, imag_freqs]) #append list of imag freqs to master imag_freqs alongside the associated filename
                        
                    except Exception as e:
                        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error extracting imaginary frequency data from {filename}: {e}. Please report this issue along with your .out file to the issues section of the Github repo.')
                        QApplication.processEvents()
                        pass

                #placeholder for relative energy
                Erel = 123.0

                #Prepare the values to be written
                values = [filename, n_imag, Eelec, ZPE, Ecorr, Hcorr, Gcorr, Total_ZPE, Total_E, Total_H, Total_S * T, Total_G, Erel,  RotA, RotB, RotC, sigma_OR, dipole_x, dipole_y, dipole_z, total_dipole, polariz, q_trans, q_rot, q_vib, q_elec, spin_contam]

            #if thermochem was not requested (like in a single point energy calculation or geometry optimization): 
            else:

                thermochem_not_requested.append(filename)
                
                #if an electronic energy is found
                if Eelec:
                    values = [filename, 'N/A', Eelec, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', RotA, RotB, RotC, 'N/A', dipole_x, dipole_y, dipole_z, total_dipole, polariz, 'N/A', 'N/A', 'N/A', 'N/A', spin_contam]

                #I don't know how a job could terminate normally but also enter this block, but YOLO. 
                else:
                    values = [filename, 'N/A', 12345.0, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', RotA, RotB, RotC, 'N/A', dipole_x, dipole_y, dipole_z, total_dipole, polariz, 'N/A', 'N/A', 'N/A', 'N/A', spin_contam]
                    print(f'Congratulations. {filename} entered a block of code that I did not think was possible. Please report this to the issues section of the Github repo.')
        
        #if the ORCA job did not finish
        else:

            #did they get to the point where they calcualted vib freqs? Sometimes the job terminates due to walltime during the CHELPG charge calculation step
            #We can check for this to see if a non-None ZPE was returned.
            #ZPE requires the extraction of all vib frequencies, and will only be updated to a non-None value if vib freqs are found. 
            if thermochem_flag and ZPE:
                abnormal_term_wVibs_list.append(filename)

            else:
                abnormal_term_list.append(filename)

            #if the job is a single point energy calculation did not finish
            if not thermochem_flag and not Eelec:
                values = [filename, 'N/A', 12345.0, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', RotA, RotB, RotC, sigma_OR, dipole_x, dipole_y, dipole_z, total_dipole, polariz, 'N/A', 'N/A', 'N/A', 'N/A', spin_contam]
            
            #if the job is an optimization that didn't finish
            else:
                values = [filename, 'N/A', Eelec, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', RotA, RotB, RotC, sigma_OR, dipole_x, dipole_y, dipole_z, total_dipole, polariz, 'N/A', 'N/A', 'N/A', 'N/A', spin_contam]
        
        values = [val if val is not None else 'N/A' for val in values]

        new_row = pd.DataFrame([values], columns=properties)
        new_row = new_row.dropna(axis=1, how='all')

        #append values to Dataframe
        df = pd.concat([df, new_row], ignore_index=True)

        #debugging for when problems arise
        '''
        variable_names = ['filename', 'n_imag', 'Eelec', 'ZPE', 'Ecorr', 'Hcorr', 'Gcorr', 'Total_ZPE', 'Total_E', 'Total_H', 'Total_S*T', 'Total_G', 'Erel', 'RotA/100', 'RotB/100', 'RotC/100', 'sigma_OR', 'dipole_x', 'dipole_y', 'dipole_z', 'total_dipole', 'polariz_value', 'q_trans', 'q_rot', 'q_vib', 'q_elec', 'spin_contam_value']
        variable_values = [v if v is not None else 'I am none' for v in values]
        with pd.ExcelWriter(os.path.join(directory, 'debugging.xlsx'), engine='openpyxl') as writer:
            
            df = pd.DataFrame({
                'Variable Name': variable_names,
                'Value': variable_values
            })

            df.to_excel(writer, sheet_name=f'{filename}', index=False)
        '''

    #Calculate the minimum Gibbs energy after all the energies have been extracted and all non entries are removed (liek the N/A strings)
    numeric_gibbs_list = np.array([float(g) for g in Gibbs_list if isinstance(g, (int, float))])

    #ensure that the gibbs list isn't empty - otherwise it won't have a minimum!
    if numeric_gibbs_list.size > 0:
        min_Gibbs = np.min(numeric_gibbs_list)
    else:
        min_Gibbs = 0

    #replace all N/A entries in the Total Gibbs, Enthalpy, entropy, ZPE, etc with infinity for sorting purposes
    columns_to_sort_energies = ['Total Gibbs Energy', 'Total Enthalpy', 'Total Thermal Energy', 'Total Entropy (T*S)', 'Total ZPE']

    for column in columns_to_sort_energies:
        df[column] = pd.to_numeric(df[column], errors='coerce')
        df[column] = df[column].fillna(np.inf)
    
    #Replace the placeholder relative energies with the actual relative energies in the DataFrame, then temporarily fill NaNs with infinity (needed for calculation of rel gibbs)
    df['Relative Gibbs Energy'] = (df['Total Gibbs Energy'] - min_Gibbs) * 2625.4996395 #Conversion for Hartree to kJ/mol - can be changed to user preference
    
    #Sort files based on user input
    if sort_by == 'G':
        #Sort the DataFrame based on 'Relative Gibbs Energy'
        df = df.sort_values(by='Relative Gibbs Energy')

    if sort_by == 'H':
        #Sort the DataFrame based on 'Total Enthalpy'
        df = df.sort_values(by='Total Enthalpy')

    if sort_by == 'E':
        #Sort the DataFrame based on 'Total Thermal Energy'
        df = df.sort_values(by='Total Thermal Energy')

    if sort_by == 'S':
        #Sort the DataFrame based on 'Total Thermal Energy'
        df = df.sort_values(by='Total Entropy (T*S)')

    if sort_by == 'Z':
        df = df.sort_values(by='Total ZPE')
    
    #replace all infinite values with N/A
    for column in columns_to_sort_energies:
        df[column] = df[column].replace(np.inf, 'N/A')

    #Create a unique output file name
    output_csv = os.path.join(directory, f'Thermo_data_{int(T)}K_{int(p)}Pa_{str(np.round(vib_scl, 4)).replace(".","-")}vibscl.csv')

    i = 2
    while os.path.isfile(output_csv):
        output_csv = os.path.join(directory, f'Thermo_data_{int(T)}K_{int(p)}Pa_{str(np.round(vib_scl, 4)).replace(".","-")}vibscl_{i}.csv')
        i += 1
    
    df.to_csv(output_csv, index=False)

    #Inform users via the print window of any abnormal terminations, missing thermochemistry, and imaginary frequencies
    if len(thermochem_not_requested) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Vibrational frequencies were not requested from the following files; N/A is being written as a placeholder for relevant quantities:\n{", ".join(thermochem_not_requested)}.\n')

    if len(abnormal_term_list) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The following files did not terminate normally, and consequently, will have "N/A" being written in place of missing values:\n{", ".join(abnormal_term_list)}.\n')
    
    if len(abnormal_term_wVibs_list) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Of the files that did not terminate normally, the following did contain vibrational frequencies, and thus, likely reached walltime during a subsequent calculation step (e.g., CHELPG charge calcualtion):\n{", ".join(abnormal_term_wVibs_list)}.\n')
    
    if len(master_imag_freq_list) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(master_imag_freq_list)} file(s) contain imaginary frequencies:')
        for filename, imag_freq in master_imag_freq_list:
            print(f'{filename}: {imag_freq}')

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Electronic energies, thermochemistry, and molecular properties have been extracted from {len(filenames)} ORCA .out files in {np.round(time.time() - start,2)} seconds.')

#external testing
if __name__ == '__main__':

    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\Manuscripts\2024\Heterobimetallic_Coordination_Derek\Calcs_ORCA6\Fragments\DFT_OptFreq'
    sort_by = 'Z'
    T = 298.15
    p = 100000.
    vib_scl = 1.00

    '''fundamental constants'''

    pi4eps = 4.*np.pi*8.8541878128e-12 #F/m, 4pi*electric constant
    kB = 1.380649e-23 #J/K boltzmann constant
    N_Av = 6.02214076e+23 #Avogadro's number
    h_SI = 6.62607015e-34 #Js, Planck constant
    h_bar = h_SI/(2.*np.pi) #reduced Planck constant
    c_SI = 299792458. #m/s, speed of light
    elC = 1.602176634e-19 #C, elementary charge
    m_el = 9.1093837015e-31 #kg, electron mass
    r_Bohr = pi4eps*h_bar**2 / (elC**2 * m_el) #m, Bohr radius

    ##unit conversions
    bohr2m = r_Bohr #m/bohr
    bohr2A = bohr2m * 1.0e+10 #Angstrom/bohr
    cm2in = 0.3937 #cm to inch
    J2Eh = (m_el * r_Bohr**2) / (h_bar**2) #1 Joule in hartree
    amu2kg = 1./N_Av * 1e-3 #amu/kg
    D2Cm = 3.336e-30 #Debye to Coulomb*m
    eigval2s2 = 1e20/amu2kg/J2Eh #Eh/Ang**2.amu = 2.62549964E+29 s**-2

    #spectroscopic unit conversions
    eV2J = elC #Conversion factor from eV to Joules
    nm2m = 1E-9 #Conversion factor from nm to meters
    cm2m = 1E-2 #Conversion factor from cm to meters

    ##fundamental constants in other units
    kB_Eh = kB*J2Eh #hartree/kelvin
    h_cm = h_SI/(h_SI*c_SI*100.) #cm**-1 * s, Planck constant in wavenumber*s

    '''dictionary to write all the constants to'''

    c = {
        #constants
        'pi4eps': pi4eps,
        'kB': kB,
        'N_Av': N_Av,
        'h_SI': h_SI,
        'h_bar': h_bar,
        'c_SI': c_SI,
        'elC': elC,
        'm_el': m_el,

        #unit conversions
        'r_Bohr': r_Bohr,
        'bohr2m': bohr2m,
        'bohr2A': bohr2A,
        'cm2in': cm2in,
        'J2Eh': J2Eh,
        'amu2kg': amu2kg,
        'D2Cm': D2Cm,
        'eigval2s2': eigval2s2,

        #spectroscopic unit conversions
        'eV2J': eV2J,
        'nm2m': nm2m,
        'cm2m': cm2m,

        #consants in other units
        'kB_Eh': kB_Eh,
        'h_cm': h_cm
    }

    #run the code
    ORCA_Thermochem_Calculator(directory,T, p, vib_scl, sort_by)