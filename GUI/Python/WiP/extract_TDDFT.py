import os, time, re
import numpy as np
import pandas as pd
from datetime import datetime

def extract_TDDFT_VG(directory, fwhm, lower_bound, upper_bound, step_size, normalization, unit, scale_freq, fixed_shift, output_excel_file):
    
    def extract_TDDFT_spectrum(file, unit, scale_freq, fixed_shift):
        '''Extracts out the spectrum block from an ORCA .out file. Function returns the freq in cm**-1 and the intensity in km/mol as arrays.'''
        
        #constants 
        h = 6.62607015E-34 #J s
        c = 299792458. # m/s
        e = 1.60217663E-19 # coulombs 

        # Open the .out file and read its contents
        with open(file, 'r') as f:
            file_content = f.read()
        
        # Define the start and end patterns using regular expressions
        start_pattern = re.compile(r'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS')
        end_pattern = re.compile(r'ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS')

        # Find the start and end indices of the TDDFT spectrum section
        start_match = start_pattern.search(file_content)
        end_match = end_pattern.search(file_content)

        if start_match and end_match:
            start_idx = start_match.end()
            end_idx = end_match.start()

            # Extract the IR spectrum block
            ir_block = file_content[start_idx:end_idx].strip()

            # Split the block into lines and extract frequency and intensity
            lines = ir_block.split('\n')
            data_lines = [line.split() for line in lines[4:-2]]
            
            # get number of excited state
            state = [line[0] for line in data_lines if len(line) == 8]
            
            # Create a local variable for unit and set it to 'nm' if it's not provided
            local_unit = unit if unit in ['eV', 'cm-1', 'nm'] else 'nm'

            # Get energy of the transition and scale it according to the unit provided
            
            # for eV
            if local_unit == 'eV':
                energy = [((h * c) / (float(line[2]) * 1E-9)) / e - float(scale_freq) for line in data_lines if len(line) == 8]  # s0 to s## transition energy, where ## is the corresponding state. Value is in eV

            #for cm**-1
            elif local_unit == 'cm-1' and not fixed_shift:
                energy = [float(line[1]) - float(scale_freq) for line in data_lines if len(line) == 8]  # s0 to s## transition energy, where ## is the corresponding state. Value is in cm-1

            elif local_unit == 'cm-1' and fixed_shift:
                orig_energy = [float(line[1]) for line in data_lines if len(line) == 8]  # s0 to s## transition energy, where ## is the corresponding state. Value is in nm
                scaled_energies = [((h * c * float(line[1]) * 100) / e) - float(scale_freq) for line in data_lines if len(line) == 8]  # s0 to s## transition energy, where ## is the corresponding state. Value is in nm
                print(orig_energy)
                energy = [(ener * e) / (h * c * 100) for ener in scaled_energies]
                print(energy)

            #for nm
            elif local_unit == 'nm' and not fixed_shift:
                energy = [float(line[2]) for line in data_lines if len(line) == 8]  # s0 to s## transition energy, where ## is the corresponding state. Value is in nm
                

            elif local_unit == 'nm' and fixed_shift:
                orig_energy = [float(line[2]) for line in data_lines if len(line) == 8]  # s0 to s## transition energy, where ## is the corresponding state. Value is in nm
                scaled_energies = [((h * c) / (float(line[2]) * 1E-9)) / e - float(scale_freq) for line in data_lines if len(line) == 8]  # s0 to s## transition energy, where ## is the corresponding state. Value is in nm
                print(orig_energy)
                energy = [(h * c * 1E9) / (ener * e) for ener in scaled_energies]
                print(energy)

            #get the intensity
            intensity = [float(line[3]) for line in data_lines] #intensity treated as oscillator strength

            return state, energy, intensity
        
        else:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The spectrum block could not be found in {os.path.basename(file)}.')
            raise ValueError

    def broaden_spectrum(freq, intensity, fwhm, lower_bound, upper_bound, step_size, normalization):
        '''Broadens each absoprtion in an IR stick strectrum with a user specified Gaussian FWHM.'''
        # Create a grid from the lower_bound to the upper_bound with 1 cm**-1 spacing
        grid = np.arange(lower_bound, upper_bound + step_size, step_size)

        #define sigma for gaussian function
        w = fwhm / (2 * np.sqrt(np.log(4)))

        # Create an empty array to store the broadened spectrum
        broadened_spectrum = np.zeros_like(grid, dtype=float)

        # Broaden the spectrum for each frequency/intensity pair
        for f, i in zip(freq, intensity):

            broadened_intensity = i * np.exp(-0.5 * np.square((grid - f) / w))

            # Add the broadened intensity to the overall spectrum
            broadened_spectrum += broadened_intensity        

        if normalization:
            # Normalize the spectrum to have a maximum value of 1
            max_value = np.max(broadened_spectrum)
            
            if max_value > 0: #avoid dividing by zero
                broadened_spectrum /= max_value

            elif max_value == 0:
                raise ValueError('The extracted spectrum has a maximum value of zero.')

            else:
                pass

        return grid, broadened_spectrum
       
    '''Main workflow'''
    # Find all .out files in the specified directory
    files = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x]

    # Create a Pandas DataFrame to store the results
    df = pd.DataFrame()

    for file in files:

        # Extract IR spectrum from the file
        state, energy, intensity = extract_TDDFT_spectrum(os.path.join(directory, file), unit, scale_freq, fixed_shift)

        # Broaden the spectrum
        energy_grid, broadened_spectrum = broaden_spectrum(energy, intensity, fwhm, lower_bound, upper_bound, step_size, normalization)

        # Add the broadened spectrum to the DataFrame with the file title as the column header
        df[os.path.basename(file)] = broadened_spectrum

    # Add the wavenumber column to the DataFrame
    df['Wavenumber'] = energy_grid

    # Reorder columns
    df = df[['Wavenumber'] + [col for col in df.columns if col != 'Wavenumber']]

    # Write the DataFrame to an Excel file
    if os.path.isfile(output_excel_file):
        i = 2
        while os.path.isfile(output_excel_file):
            output_excel_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq)}_unit_{unit}_v{i}.xlsx')
            i += 1    
    try:
        df.to_excel(output_excel_file, index=False)
    except (PermissionError):
        i = 2

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} A file with the same name is open in Excel. Writing data to a new file name.')
        output_excel_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq)}_unit_{unit}_v{i}.xlsx')
        
        while os.path.isfile(output_excel_file):
            output_excel_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq)}_unit_{unit}_v{i}.xlsx')
            i += 1
        
        df.to_excel(output_excel_file, index=False)

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The new output file is {os.path.basename(output_excel_file)}')
        
    stime = time.time()
    ftime = time.time()
    
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {len(files)} .out files have been processed in {np.round(ftime - stime, 2)} seconds.')

if __name__ == "__main__":
    # Specify the directory containing .out files
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Manual\Appendix_A\Testing\TDDFT_plot'

    # Specify the FWHM for broadening
    fwhm = 100

    # Specify the lower and upper bounds
    lower_bound = 25000
    upper_bound = 67000
    step_size = 100

    #harmonic scaling factor
    scale_freq = 1610

    #fixed shift in eV if nm is being used
    fixed_shift = True
    if fixed_shift:
       scale_freq = 0.2 

    #option to normalize all calculated spectra to a maximum intensity of 1
    normalization = True

    unit = 'cm-1' #eV, nm, or cm-1

    # Specify the output Excel file
    output_excel_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq)}_unit_{unit}.xlsx')

    # Process the .out files and write the results to an Excel file
    extract_TDDFT_VG(directory, fwhm, lower_bound, upper_bound, step_size, normalization, unit, scale_freq, fixed_shift, output_excel_file)
