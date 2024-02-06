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
        dipole, polariz, dipole_ax, vibs_array, ZPE = None, None, None, None, None
        n_imag = 0
        imag_freqs = None

        #get info from files
        for line in data:
            
            if line.startswith(' Total Charge'):
                charge = float(line.split()[-1])

            if line.startswith(' Multiplicity'):
                multi = float(line.split()[-1])

            if line.startswith('FINAL SINGLE POINT ENERGY'):
                Eelec = float(line.split()[-1])

            if line.startswith('Rotational constants in cm-1'):
                RotA = float(line.split()[4]) * 100.  # in m**-1
                RotB = float(line.split()[5]) * 100.
                RotC = float(line.split()[6]) * 100.
                RotABC = np.array([RotA, RotB, RotC])  # in m**-1

            if line.startswith('Point Group:'):
                sigma_OR = int(line.split()[-1])  # integer

            if line.startswith('Total Mass          ...'):
                mass = float(line.split()[-2])
                m_SI = mass * c['amu2kg']

            if line.startswith('Magnitude (Debye)'):
                dipole = float(line.split()[-1])  # in Debye

            if line.startswith('Isotropic polarizability :'):
                polariz = float(line.split()[-1]) * 1.4818e-31  # in m^3

            if line.startswith('x,y,z [Debye]'):
                mu_abc = np.array([float(s) for s in line.split()[-3:]])
                dipole_ax = ['A', 'B', 'C'][np.argmax(np.abs(mu_abc))]

            if line.startswith('freq.   '):
                if vibs_array == None:
                    vibs_array = []
                
                vib_str = line.split()[1]  # Extract the string
                try:
                    vib = float(line.split()[1]) * vib_scl
                except TypeError:
                    raise ValueError(f"Error: Unable to convert '{vib_str}' to a float from {os.path.basename(file)}.")
                    
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

        return charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, dipole, polariz, dipole_ax, vibs_array, ZPE, Total_ZPE, n_imag, imag_freqs
    
    def calc_partition_function(RotABC, sigma_OR, vibs_array, multi, T, p):
        '''Calculates the partition functions (trans, rot, vib, elec) for the molecule at a given p, T.'''
        q_trans = (2. * np.pi * m_SI * c['kB'] * T / c['h_SI'] ** 2) ** 1.5 * c['kB'] * T / p

        q_rot = 1. / sigma_OR * (c['kB'] * T / (c['h_SI'] * c['c_SI'])) ** 1.5 * (np.pi / (np.prod(RotABC))) ** 0.5

        u = c['h_SI'] * c['c_SI'] * vibs_array * 100. / (c['kB'] * T) #vib partition functional can be scaled if needed
        q_vib = np.prod(1. / (1. - np.exp(-u)))

        q_elec = multi

        return [q_trans, q_rot, q_vib, q_elec]


    def calc_thermochemistry(vibs_array, q_trans, q_rot, q_vib, q_elec, ZPE, Eelec, T, p):
        '''Computes standard thermochemical functions at given p, T.'''

        theta_v = c['h_SI'] * c['c_SI'] * vibs_array * 100. / c['kB']
        theta_v_oT = theta_v / T

        S_T = c['kB_Eh'] * (np.log(q_trans) + 1. + 3. / 2.) #translational contribution to entropy
        S_R = c['kB_Eh'] * (np.log(q_rot) + 3. / 2.) #rotational contribution to entropy
        S_V = c['kB_Eh'] * np.sum(theta_v_oT / (np.exp(theta_v_oT) - 1.) - np.log(1. - np.exp(-theta_v_oT))) #vibrational contribution to entropy
        S_E = c['kB_Eh'] * np.log(q_elec) #electronic contribution to entropy
        Total_S = S_T + S_R + S_V + S_E

        #print(S_T*T, S_R*T, S_V*T, S_E*T)

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

    start = time.time()

    thermo_properties = [
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
        'Trans. Part. func.',
        'Rot. Part. func.',
        'Vib. Part. function',
        'Elec. Part. func.'
    ]

    #get list of filenames, and check if the directory does not contain any .out files
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x.lower()]

    if len(filenames) == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} There are no .out files in the provided directory.')
        return

    # Format the header for consistent spacing 
    header = '{},\n'.format(','.join(['{:<25}'] * len(thermo_properties)))

    # Create output file and write header to it
    output_csv = os.path.join(directory, f'Thermo_data_{int(T)}K_{int(p/1000)}kPa_{str(vib_scl).replace(".","-")}vibscl.csv')

    i = 1
    while os.path.isfile(output_csv):
        output_csv = os.path.join(directory, f'Thermo_data_{int(T)}K_{int(p)}Pa_{str(vib_scl).replace(".","-")}vibscl_{i}.csv')
        i += 1
    
    with open(output_csv, 'w') as opf:
        opf.write(header.format(*thermo_properties))

    Gibbs_list = []
    missing_thermochem = []
    master_imag_freq_list = []

    for filename in filenames:

        #get data from input file
        charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, dipole, polariz, dipole_ax, vibs_array, ZPE, Total_ZPE, n_imag, imag_freqs = read_molecule_data(os.path.join(directory, filename), vib_scl)
        
        #check for missing thermochemistry. If all values are there, proceed with thermochem analysis 
        if all(x is not None and (not isinstance(x, np.ndarray) or np.all(x)) for x in [charge, multi, Eelec, RotABC, sigma_OR, mass, m_SI, dipole, dipole_ax, vibs_array, ZPE]):
        #will return false if any element in list is None
        #do not check for polarization because that calculation is optional and not done every time

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

            # Prepare the values to be written
            values = [filename, n_imag, Eelec, ZPE, Ecorr, Hcorr, Gcorr, Total_ZPE, Total_E, Total_H, Total_S * T, Total_G, Erel, q_trans, q_rot, q_vib, q_elec]

            # Create a format string for consistent spacing
            format_str = '{}'.format(','.join(['{:<25}'] * len(values)) + ',\n')

            # Append the extracted properties to the CSV file
            with open(output_csv, 'a') as opf:
                opf.write(format_str.format(*values))

        else:
            #print(f'{filename} is missing thermochemistry. Writing 12345 as a placeholder for missing value')
            missing_thermochem.append(filename)
            values = [filename, -12345.0, Eelec, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0, -12345.0]
            
            # Create a format string for consistent spacing
            format_str = '{}'.format(','.join(['{:<25}'] * len(values)) + ',\n')

            with open(output_csv, 'a') as opf:
                opf.write(format_str.format(*values))

    # Calculate the minimum Gibbs energy
    min_Gibbs = np.min(Gibbs_list)

    # Read the CSV into a pandas DataFrame
    df = pd.read_csv(output_csv)

    # Calculate relative energy column and update it in the DataFrame
    df['Relative Gibbs Energy    '] = (df['Total Gibbs Energy       '] - min_Gibbs) * 2625.5

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
        opf.write(header.format(*thermo_properties))  # Write the header first

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
    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Manual\Appendix_A\GUI\Python\Alex_stuff\Thermochemistry'
    sort_by = 'G'
    T = 298.15
    p = 100000.
    vib_scl = 0.95
    ORCA_Thermochem_Calculator(directory, sort_by, T, p, vib_scl)




