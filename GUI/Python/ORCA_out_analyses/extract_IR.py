import os, re, time
import numpy as np
import pandas as pd
from datetime import datetime
from PyQt6.QtWidgets import QApplication

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('QtAgg') #Use matplotlib backend that is compatible w/ PyQt6 to prevent competition w/ GUI event loop

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams.update({'lines.linewidth': 1.5})

def extract_IR_spectra(directory, fwhm, lower_bound, upper_bound, step_size, IR_output_basename, plotting, scale_freq = 1.):

    def extract_ir_freqs(file, scale_freq):
        '''Extracts out the IR spectrum block from an ORCA .out file. Function returns the freq in cm**-1 and the intensity in km/mol as arrays.'''
        
        #open the .out file and read its contents
        
        with open(file, 'r') as opf:
            file_content = opf.read()
        
        #Find the start and end indices of the IR spectrum section
        start_idx = file_content.find('-----------\nIR SPECTRUM\n-----------')
        end_idx = file_content.find('--------------------------\nTHERMOCHEMISTRY AT')

        if start_idx == -1 or end_idx == -1:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} IR data could not be found in {os.path.basename(file)}. Processing of this file will be skipped.')
            return None, None
        
        else:
            #Extract the IR spectrum block
            ir_block = file_content[start_idx:end_idx].strip()

            #Split the block into lines and extract frequency and intensity
            lines = ir_block.split('\n')
            data_lines = [line.split() for line in lines[7:-6]]
            try:
                freq = [float(line[1]) * scale_freq for line in data_lines]
                intensity = [float(line[3]) for line in data_lines]
                return freq, intensity
            
            except Exception as e:
                #print(f'{datetime.now().strftime("[ %H:%M:%S ]")} IR data from {os.path.basename(file)} could not be processed: {e}.')
                return None, None

    def extract_gibbs_energy(file):
        '''Extracts the Gibbs energy and temperature from the ORCA .out file.'''
        
        gibbs, temp = None, None

        with open(file, 'r') as f:
            for line in f:
                if 'Final Gibbs free energy' in line:
                    try:
                        gibbs = float(line.split()[-2])
                    except ValueError:
                        break
                
                if 'THERMOCHEMISTRY AT' in line:
                    try:
                        temp = float(line.split()[2].replace('K', ''))
                    except (ValueError, IndexError):
                        break
            
            if gibbs is not None and temp is not None:
                return gibbs, temp
            else:
                return None, None
            
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No Gibbs free energy found in {os.path.basename(file)}.')
        
        return None
                
    def broaden_spectrum(file, freq, intensity, fwhm, lower_bound, upper_bound, step_size):
        '''Broadens each absoprtion in an IR stick strectrum with a user specified Gaussian FWHM.
        Returns vib freq., intensity, and normalized intensity'''
        #Create a grid from the lower_bound to the upper_bound with 1 cm**-1 spacing
        freq_grid = np.arange(lower_bound, upper_bound + step_size, step_size)

        #define sigma for gaussian function
        w = fwhm / (2 * np.sqrt(np.log(4)))

        #Create an empty array to store the broadened spectrum
        broadened_spectrum = np.zeros_like(freq_grid, dtype=float)

        #Broaden the spectrum for each frequency/intensity pair
        for f, i in zip(freq, intensity):

            broadened_intensity = i * np.exp(-0.5 * np.square((freq_grid - f) / w))

            #Add the broadened intensity to the overall spectrum
            broadened_spectrum += broadened_intensity        

        #Get the maximum value of the spectrum to prepare for normalization
        max_value = np.max(broadened_spectrum)
        
        if np.isclose(max_value, 0, atol=1E-3):
           print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The extracted spectrum has a maximum value of zero. Please check the IR block {os.path.basename(file)} for correctness.')
           return None, None

        #normalize spectrum if requested - but need to return both for Boltzmann weighting
        norm_broadened_spectrum = broadened_spectrum / max_value 
        return freq_grid, broadened_spectrum, norm_broadened_spectrum
       
    '''Main workflow'''
    
    stime = time.time()
    df_norm = pd.DataFrame()  # For normalized spectra
    df_unorm = pd.DataFrame()  # For unnormalized spectra

    #get list of filenames, and check if the directory does not contain any .out files
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x.lower()]

    if len(filenames) == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} There are no .out files in the provided directory.')
        return

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting extraction of IR spectra from {len(filenames)} ORCA .out files...')
    QApplication.processEvents()

    # Initialize lists to store data only for successfully processed files
    successful_filenames = []
    successful_energies = []
    successful_temps = []
    failed_files = []

    for file in filenames:

        #Extract IR spectrum and gibbs energy from the file - extract_ir_freqs returns None for freq, intensity if an .out file does not contain any IR freq data
        freq, intensity = extract_ir_freqs(os.path.join(directory, file),scale_freq)
        energy, temp = extract_gibbs_energy(os.path.join(directory, file))

        #Broaden the spectrum only if IR data was sucessfuly extracted from the file
        if freq is not None and intensity is not None and energy is not None:
            freq_grid, unorm_spectrum, norm_spectrum = broaden_spectrum(file, freq, intensity, fwhm, lower_bound, upper_bound, step_size)

            #Add the broadened spectrum to the DataFrame with the file title as the column header - broaden_spectrum returns None for freq_grid, broadened_spectrum if the max intensity of all IR freqs is zero. 
            if freq_grid is not None:
                df_unorm[os.path.basename(file)] = unorm_spectrum
                df_norm[os.path.basename(file)] = norm_spectrum
                successful_filenames.append(os.path.basename(file))
                successful_energies.append(energy)
                successful_temps.append(temp)

            else:
                failed_files.append(file)
                continue
        
        #If IR data cannot be found and/or the maximum intensity of the IR freqs is zero, skip to the next file
        else:
            failed_files.append(file)
            continue

    #Compute the Boltzmann Weights
    emin = min(successful_energies)
    tmin = min(successful_temps) #use minimum temp from all temps. These should all be consistent anyways
    rel_energies = [(E - emin) * 2625.5 for E in successful_energies]
    populations = [np.exp(-E / (8.3145E-3 * tmin)) for E in rel_energies]
    total_population = sum(populations)
    relative_populations = [pop / total_population for pop in populations]
   
   # Collect data for the populations sheet - need to use "successful" lists to ensure that each array is the same size. 
    population_data = {
        'Filename': [f for f in successful_filenames],
        'Gibbs Energy (Eh)': successful_energies,
        'Relative Energy (kJ/mol)': rel_energies,
        'Population': populations,
        'Relative Population': relative_populations
    }
    df_populations = pd.DataFrame(population_data)

    # Calculate Boltzmann-weighted spectrum
    bw_spectrum = np.zeros_like(freq_grid, dtype=float)
    for i, col in enumerate(df_unorm.columns[1:]):  # Skip first column, which is wavenumber
        bw_spectrum += relative_populations[i] * df_unorm[col].to_numpy()

    norm_bw_spectrum = bw_spectrum / np.max(bw_spectrum)
    
    # Generate unique output file name to prevent overwriting
    i = 2 
    output_file = os.path.join(directory, f'{IR_output_basename}_scl_{scale_freq}_fwhm_{fwhm}cm.xlsx')
    while os.path.isfile(output_file):
        output_file = os.path.join(directory, f'{IR_output_basename}_scl_{scale_freq}_fwhm_{fwhm}cm_v{i}.xlsx')
        i += 1     

    # Add Boltzmann-weighted spectrum to a new sheet in the Excel file
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        
        #write unnormalizaed spectra
        df_unorm['Wavenumber'] = freq_grid
        df_unorm = df_unorm[['Wavenumber'] + [col for col in df_unorm.columns if col != 'Wavenumber']]
        df_unorm.to_excel(writer, sheet_name='Unnormalized Spectra', index=False)

        #write normalized spectra
        df_norm['Wavenumber'] = freq_grid
        df_norm = df_norm[['Wavenumber'] + [col for col in df_norm.columns if col != 'Wavenumber']]
        df_norm.to_excel(writer, sheet_name='Normalized Spectra', index=False)

        #write boltzmann weighted spectra (norm and unorm)
        pd.DataFrame({'Wavenumber': freq_grid, 'BW Spectrum': bw_spectrum, 
                      'Normalized BW Spectrum': norm_bw_spectrum}).to_excel(writer, sheet_name='Boltzmann Weighted', index=False)

        #write populations of each isomer
        df_populations.to_excel(writer, sheet_name=f'Populations @ {int(tmin)}K', index=False)

    #Plotting
    if plotting:

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Extracted IR spectra will be plotted momentarily...')
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

        #initialize list to store filenames for generating plot legend
        labels = []

        # Plot unnormalized spectra in the first subplot
        for column in df_unorm.columns:
            if column != 'Wavenumber':
                ax1.plot(df_unorm['Wavenumber'], df_unorm[column], label=column)

                # Generate a short label for the legend if needed
                if len(column) > 20:
                    match = re.search(r"(\d+)(?!.*\d)", column)  # Regex to find the last number
                    basename = column.split('_')[0]
                    
                    if match and basename:
                        new_label = f'{basename}_{match.group(1)}'
                    elif match and not basename:
                        new_label = f'File_{match.group(1)}'
                    else:
                        new_label = column[:int(np.ceil(len(column) * 0.25))]  # Use 25% of filename
                        
                else:
                    new_label = column

                labels.append(new_label)

        # Customize the first subplot
        ax1.set_ylabel('Intensity (Unnormalized)')
        ax1.set_title('Unnormalized IR Spectra')
        ax1.legend(labels=labels, loc='upper left', bbox_to_anchor=(1, 1))

        # Plot the normalized Boltzmann-weighted spectrum in the second subplot
        ax2.plot(freq_grid, norm_bw_spectrum, color='purple', linewidth=2, label='Normalized BW Spectrum')
        
        # Customize the second subplot
        ax2.set_xlabel('Wavenumber (cm⁻¹)')
        ax2.set_ylabel('Intensity (Normalized)')
        ax2.legend(loc='upper left')

        # Adjust layout and save the plot
        plt.tight_layout()
        plt_file = output_file.replace('.xlsx', '.png')

        #save the plot and show it to the user
        try:
            plt.savefig(plt_file)
        except FileExistsError:
            print(f'{datetime.now().strftime("[ %H:%M:%S ]")} {os.path.basename(plt_file)} already exists within the output directory. You can save the plot manually from the plot that shows.')
            pass

        #show the plot
        plt.tight_layout()
        plt.show()

    if len(failed_files) > 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} IR spectra could not be processed from {len(failed_files)} files. These have been omitted from the IR extraction:\n{",\n".join(failed_files)}')

    ftime = time.time()
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} IR spectra from {len(filenames)} .out files have been processed in {np.round(ftime - stime, 2)} seconds.')

#external testing
if __name__ == "__main__":
    #Specify the directory containing .out files
    directory = r'E:\OneDrive - University of Waterloo\Waterloo\Manuscripts\2024\Leipzig_PIC\Calcs_Iteration2\Met-OMe\RS_Met-OMe'

    #Specify the FWHM for broadening
    fwhm_value = 3

    #Specify the lower and upper bounds
    lower_bound_value = 500 #must be an integer
    upper_bound_value = 3600 #must be an integer
    step_size = 1 #must be an integer

    #harmonic scaling factor
    scale_freq = 1.00

    #option to plot spectra
    plotting = False

    #basename for output file
    out_basename = 'test'

    #Process the .out files and write the results to an Excel file
    extract_IR_spectra(directory, fwhm_value, lower_bound_value, upper_bound_value, step_size, out_basename, plotting, scale_freq)
