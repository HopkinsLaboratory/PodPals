import os, re, time
import numpy as np
import pandas as pd
from datetime import datetime
from Python.constants_and_conversions import c

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
        vibs_array, ZPE = None, None
        n_imag = 0
        imag_freqs = None

        #get info from files
        for line in data:
            
            #charge
            if line.startswith(' Total Charge'):
                charge = float(line.split()[-1])
            
            #multiplicity
            if line.startswith(' Multiplicity'):
                multi = float(line.split()[-1])
            
            # final electronic energy
            if line.startswith('FINAL SINGLE POINT ENERGY'):
                Eelec = float(line.split()[-1])

            #rotational constants
            if line.startswith('Rotational constants in cm-1'):
                RotA = float(line.split()[4]) * 100.  # in m**-1
                RotB = float(line.split()[5]) * 100.
                RotC = float(line.split()[6]) * 100.
                RotABC = np.array([RotA, RotB, RotC])  # in m**-1

            #point group
            if line.startswith('Point Group:'):
                sigma_OR = int(line.split()[-1])  # integer

            #mass
            if line.startswith('Total Mass          ...'):
                mass = float(line.split()[-2])
                m_SI = mass * c['amu2kg']

            #dipole moment
            if line.startswith('Magnitude (Debye)'):
                total_dipole = float(line.split()[-1])  # in Debye

            #dipole along x, y, and z axes
            if line.startswith('Total Dipole Moment    :'):
                
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
            if line.startswith('Isotropic polarizability :'):
                polariz = float(line.split()[-1]) * 1.4818e-31  # in m^3

            #dipole moment along the rotational axes
            if line.startswith('x,y,z [Debye]'):
                mu_abc = np.array([float(s) for s in line.split()[-3:]])
                
                #max dipole along each of the rotational axes
                dipole_ax = ['A', 'B', 'C'][np.argmax(np.abs(mu_abc))]

            #vibrational freqs
            if line.startswith('freq.   '):
                if vibs_array == None:
                    vibs_array = []
                
                vib_str = line.split()[1]  # Extract the string
                try:
                    vib = float(line.split()[1]) * vib_scl
                except TypeError:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error: Unable to convert {vib_str} to a float from {os.path.basename(file)}.')
                    continue
                    
                if vib > 0:
                    vibs_array.append(vib)  # in cm**-1
                else: #increase counter for imaginary freqs if any freqs are < 0
                    if imag_freqs == None:
                        imag_freqs = []
                
                    imag_freqs.append(vib)
                    n_imag += 1

        if vibs_array is not None:
            vibs_array = np.array(vibs_array)
            ZPE = 0.5 * c['h_SI'] * c['c_SI'] * 100. * np.sum(vibs_array) * c['J2Eh']  # Hartree
            Total_ZPE = Eelec + ZPE

        return charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, total_dipole, dipole_x, dipole_y, dipole_z, polariz, dipole_ax, vibs_array, ZPE, Total_ZPE, n_imag, imag_freqs
        
    def calc_partition_function(RotABC, sigma_OR, vibs_array, multi, T, p):
        '''Calculates the partition functions (trans, rot, vib, elec) for the molecule at a given p, T.'''
        q_trans = (2. * np.pi * m_SI * c['kB'] * T / c['h_SI'] ** 2) ** 1.5 * c['kB'] * T / p

        q_rot = 1. / sigma_OR * (c['kB'] * T / (c['h_SI'] * c['c_SI'])) ** 1.5 * (np.pi / (np.prod(RotABC))) ** 0.5

        u = c['h_SI'] * c['c_SI'] * vibs_array * 100. / (c['kB'] * T)
        q_vib = np.prod(1. / (1. - np.exp(-u)))

        q_elec = multi

        return [q_trans, q_rot, q_vib, q_elec]

    def calc_thermochemistry(vibs_array, q_trans, q_rot, q_vib, q_elec, ZPE, Eelec, T, p):
        '''Computes standard thermochemical functions at given p, T.
        Methodology follows that of “Molecular Thermodynamics” by McQuarrie and Simon (1999), which was summarized by Ochterski in 2000 (https://gaussian.com/wp-content/uploads/dl/thermo.pdf), and uses slightly different method for vibrational contibutions compared to the methods used in ORCA.
        Consequently, the thermochemical corrections computed here will be slightly different than the printouts in the ORCA .out file. Both are valid, although these values generated by this code more generally applicable.'''

        theta_v = c['h_SI'] * c['c_SI'] * vibs_array * 100. / c['kB']
        theta_v_oT = theta_v / T

        S_T = c['kB_Eh'] * (np.log(q_trans) + 1. + 3. / 2.) #translational contribution to entropy
        S_R = c['kB_Eh'] * (np.log(q_rot) + 3. / 2.) #rotational contribution to entropy
        S_V = c['kB_Eh'] * np.sum(theta_v_oT / (np.exp(theta_v_oT) - 1.) - np.log(1. - np.exp(-theta_v_oT))) #vibrational contribution to entropy - slightly different than the default method used by ORCA
        S_E = c['kB_Eh'] * np.log(q_elec) #electronic contribution to entropy
        Total_S = S_T + S_R + S_V + S_E

        E_T = 3. / 2. * c['kB_Eh'] * T #translational contribution to energy
        E_R = 3. / 2. * c['kB_Eh'] * T #rotational contribution to energy
        E_V = c['kB_Eh'] * np.sum(theta_v / (np.exp(theta_v_oT) - 1.)) #vibrational contribution to energy
        E_E = 0. #always zero for systems in their electronic ground state
        
        Ecorr = E_T + E_R + E_V + E_E #total energy
        Total_E = Ecorr + ZPE + Eelec

        Hcorr = (Ecorr + c['kB_Eh'] * T) + ZPE
        Total_H = Hcorr + Eelec

        Gcorr = ((Hcorr - ZPE) - T * Total_S) + ZPE
        Total_G = Gcorr + Eelec

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
        'Rot. const. A',
        'Rot. const. B',
        'Rot. const. C',
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
    ]

    # Format the header for consistent spacing 
    header = '{},\n'.format(','.join(['{:<25}'] * len(properties)))

    # Create output file and write header to it, ensuring that previous files of the same name are not overwritten
    output_csv = os.path.join(directory, f'Thermo_data_{int(T)}K_{int(p/1000)}kPa_{str(vib_scl).replace(".","-")}vibscl.csv')

    i = 2
    while os.path.isfile(output_csv):
        output_csv = os.path.join(directory, f'Thermo_data_{int(T)}K_{int(p)}Pa_{str(vib_scl).replace(".","-")}vibscl_{i}.csv')
        i += 1
    
    try:
        with open(output_csv, 'w') as opf:
            opf.write(header.format(*properties))
    
    except IOError as e:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error was encountered when writing to {os.path.basename(output_csv)}: {e}.\nFile processing will not proceed.')
        return

    Gibbs_list = []
    missing_thermochem = []
    master_imag_freq_list = []

    for filename in filenames:

        #get data from input file
        charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, total_dipole, dipole_x, dipole_y, dipole_z, polariz, dipole_ax, vibs_array, ZPE, Total_ZPE, n_imag, imag_freqs = read_molecule_data(os.path.join(directory, filename), vib_scl)

        #Check 1: if any item above is None, then thermochem was not extracted correctly
        #Check 2: if any item above in not None but is a numpy array [isinstance(x, np.ndarray)]-  ensure that none of its elements are non-zero and not None [np.all()]. If they are, vib freqs were not extracted properly
        if all(x is not None and (not isinstance(x, np.ndarray) or np.all(x)) for x in [charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, total_dipole, dipole_x, dipole_y, dipole_z, dipole_ax, vibs_array, ZPE]):

            #calculate parition functions from input data
            q_trans, q_rot, q_vib, q_elec = calc_partition_function(RotABC, sigma_OR, vibs_array, multi, T, p)

            #calculate thermochemical corrections
            Ecorr, Total_E, Hcorr, Total_H, Gcorr, Total_G, Total_S = calc_thermochemistry(vibs_array, q_trans, q_rot, q_vib, q_elec, ZPE, Eelec, T, p)
            Gibbs_list.append(Total_G)

            #if a file contians imaginary frequencies, extract their values and write them to a list for printing laster on. This is useful for determining whether a DFT job needs to be resubmitted.
            if n_imag > 0:
            
                try:
                    master_imag_freq_list.append([filename, imag_freqs]) #append list of imag freqs to master imag_freqs alongside the associated filename
                    
                except Exception as e:
                    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Error extracting imaginary frequency data from {filename}: {e}')
                    QApplication.processEvents()
                    pass

            #placeholder for relative energy
            Erel = 123.0

            # Check if polariz is None and replace it with a placeholder if so
            polariz_value = 'N/A' if polariz is None else polariz

            # Prepare the values to be written
            values = [filename, n_imag, Eelec, ZPE, Ecorr, Hcorr, Gcorr, Total_ZPE, Total_E, Total_H, Total_S * T, Total_G, Erel, RotABC[0], RotABC[1], RotABC[2], sigma_OR, dipole_x, dipole_y, dipole_z, total_dipole, polariz_value, q_trans, q_rot, q_vib, q_elec]

        #if the file is missing thermochemistry, write -12345.0 as a placeholder'
        else:
            # Check if polariz is None and replace it with a placeholder if so
            polariz_value = 'N/A' if polariz is None else polariz
            
            missing_thermochem.append(filename)
            values = [filename, -12345.0, Eelec, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, RotABC[0], RotABC[1], RotABC[2], sigma_OR, dipole_x, dipole_y, dipole_z, total_dipole, polariz_value, q_trans, q_rot, q_vib, q_elec]
            
        # Create a format string for consistent spacing, then write the data to the output.csv
        format_str = '{}'.format(','.join(['{:<25}'] * len(values)) + ',\n')

        with open(output_csv, 'a') as opf:
            opf.write(format_str.format(*values))

    # Calculate the minimum Gibbs energy after all the energies have been extracted and are located in a convienent spot
    min_Gibbs = np.min(Gibbs_list)

    # Read the CSV into a pandas DataFrame
    df = pd.read_csv(output_csv)

    # Replace the placeholder relative energies with the actual relative energies in the DataFrame
    df['Relative Gibbs Energy    '] = (df['Total Gibbs Energy       '] - min_Gibbs) * 2625.4996395 #Conversion for Hartree to kJ/mol - can be changed to user preference

    #Sort files based on user input
    if sort_by == 'G':
        # Sort the DataFrame based on 'Relative Gibbs Energy'
        df = df.sort_values(by='Relative Gibbs Energy    ')

    if sort_by == 'H':
        # Sort the DataFrame based on 'Total Enthalpy'
        df = df.sort_values(by='Total Enthalpy           ')

    if sort_by == 'E':
        # Sort the DataFrame based on 'Total Thermal Energy'
        df = df.sort_values(by='Total Thermal Energy     ')

    if sort_by == 'S':
        # Sort the DataFrame based on 'Total Thermal Energy'
        df = df.sort_values(by='Total Entropy            ')

    if sort_by == 'Z':
        df = df.sort_values(by='Total ZPE                ')

    # Write the updated DataFrame back to the CSV with consistent spacing
    with open(output_csv, 'w') as opf:
        opf.write(header.format(*properties))  # Write the header first

        for _, row in df.iterrows():
            values = list(row)
            #print(values)
            opf.write(format_str.format(*values))

    format_missing_thermochem = str('\n'.join(missing_thermochem))
    
    if len(missing_thermochem) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(missing_thermochem)} files are missing thermochemistry. -12345.0 is being written as a placeholder for the following files:\n{format_missing_thermochem}\n')

    if len(master_imag_freq_list) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(imag_freqs)} file(s) contain imaginary frequencies:')
        for filename, imag_freq in master_imag_freq_list:
            print(f'{filename}: {imag_freq}')

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Electronic energies and thermochemistry from {len(filenames)} ORCA .out files have been processed in {np.round(time.time() - start,2)} seconds.')

#external testing
if __name__ == "__main__":
    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T7_ORCAThermochemAnalyzer'
    sort_by = 'F'
    T = 298.15
    p = 100000.
    vib_scl = 1.00

    '''fundamental constants'''

    pi4eps = 4.*np.pi*8.8541878128e-12 # F/m, 4pi*electric constant
    kB = 1.380649e-23 # J/K boltzmann constant
    N_Av = 6.02214076e+23 # Avogadro's number
    h_SI = 6.62607015e-34 # Js, Planck constant
    h_bar = h_SI/(2.*np.pi) # reduced Planck constant
    c_SI = 299792458. # m/s, speed of light
    elC = 1.602176634e-19 # C, elementary charge
    m_el = 9.1093837015e-31 # kg, electron mass
    r_Bohr = pi4eps*h_bar**2 / (elC**2 * m_el) # m, Bohr radius

    ## unit conversions
    bohr2m = r_Bohr # m/bohr
    bohr2A = bohr2m * 1.0e+10 # Angstrom/bohr
    cm2in = 0.3937 # cm to inch
    J2Eh = (m_el * r_Bohr**2) / (h_bar**2) # 1 Joule in hartree
    amu2kg = 1./N_Av * 1e-3 # amu/kg
    D2Cm = 3.336e-30 # Debye to Coulomb*m
    eigval2s2 = 1e20/amu2kg/J2Eh # Eh/Ang**2.amu = 2.62549964E+29 s**-2

    #spectroscopic unit conversions
    eV2J = elC # Conversion factor from eV to Joules
    nm2m = 1E-9 # Conversion factor from nm to meters
    cm2m = 1E-2 # Conversion factor from cm to meters

    ## fundamental constants in other units
    kB_Eh = kB*J2Eh # hartree/kelvin
    h_cm = h_SI/(h_SI*c_SI*100.) # cm**-1 * s, Planck constant in wavenumber*s

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


