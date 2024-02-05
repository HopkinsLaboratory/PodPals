import os, time
import numpy as np
import pandas as pd
from datetime import datetime

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('QtAgg') #Use matplotlib backend that is compatible w/ PyQt6 to prevent competition w/ GUI event loop

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams.update({'lines.linewidth': 1.5})

def extract_IR_spectra(directory, fwhm, lower_bound, upper_bound, step_size, normalization, plotting, scale_freq = 1.):

    def extract_ir_freqs(file, scale_freq):
        '''Extracts out the IR spectrum block from an ORCA .out file. Function returns the freq in cm**-1 and the intensity in km/mol as arrays.'''
        
        #open the .out file and read its contents
        with open(file, 'r') as f:
            file_content = f.read()
        
        # Find the start and end indices of the IR spectrum section
        start_idx = file_content.find('-----------\nIR SPECTRUM\n-----------')
        end_idx = file_content.find('--------------------------\nTHERMOCHEMISTRY AT')

        # Extract the IR spectrum block
        ir_block = file_content[start_idx:end_idx].strip()

        # Split the block into lines and extract frequency and intensity
        lines = ir_block.split('\n')
        data_lines = [line.split() for line in lines[7:-6]]
        freq = [float(line[1]) * scale_freq for line in data_lines]
        intensity = [float(line[3]) for line in data_lines]

        return freq, intensity

    def broaden_spectrum(freq, intensity, fwhm, lower_bound, upper_bound, step_size, normalization):
        '''Broadens each absoprtion in an IR stick strectrum with a user specified Gaussian FWHM.'''
        # Create a grid from the lower_bound to the upper_bound with 1 cm**-1 spacing
        freq_grid = np.arange(lower_bound, upper_bound + step_size, step_size)

        #define sigma for gaussian function
        w = fwhm / (2 * np.sqrt(np.log(4)))

        # Create an empty array to store the broadened spectrum
        broadened_spectrum = np.zeros_like(freq_grid, dtype=float)

        # Broaden the spectrum for each frequency/intensity pair
        for f, i in zip(freq, intensity):

            broadened_intensity = i * np.exp(-0.5 * np.square((freq_grid - f) / w))

            # Add the broadened intensity to the overall spectrum
            broadened_spectrum += broadened_intensity        

        # Get the maximum value of the spectrum to prepare for normalization
        max_value = np.max(broadened_spectrum)
        
        if max_value == 0:
            raise ValueError('The extracted spectrum has a maximum value of zero. Something went wrong with the data extraction. Please check your output file.')

        if normalization:
            broadened_spectrum /= max_value #normalize spectrum

        return freq_grid, broadened_spectrum
       
    '''Main workflow'''
    stime = time.time()
    
    # Find all .out files in the specified directory
    files = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x]

    # Create a Pandas DataFrame to store the results
    df = pd.DataFrame()

    for file in files:

        # Extract IR spectrum from the file
        freq, intensity = extract_ir_freqs(os.path.join(directory, file),scale_freq)

        # Broaden the spectrum
        freq_grid, broadened_spectrum = broaden_spectrum(freq, intensity, fwhm, lower_bound, upper_bound, step_size, normalization)

        # Add the broadened spectrum to the DataFrame with the file title as the column header
        df[os.path.basename(file)] = broadened_spectrum

    # Add the wavenumber column to the DataFrame
    df['Wavenumber'] = freq_grid

    # Reorder columns
    df = df[['Wavenumber'] + [col for col in df.columns if col != 'Wavenumber']]

    # Write the DataFrame to an Excel file, unsuring that it does not overwrite previous data
    output_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}cm.xlsx')
    
    i = 2 
    while os.path.isfile(output_file):
        output_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}cm_v{i}.xlsx')
        i += 1

    #write data to file
    df.to_excel(output_file, index=False)

    #Plotting
    if plotting:
        plt.figure(figsize=(10, 6))

        # Loop through each spectrum in the DataFrame (excluding 'Wavenumber')
        for column in df.columns:
            if column != 'Wavenumber':
                plt.plot(df['Wavenumber'], df[column], label=column)

        # Customizing the plot
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel('Intensity')
        plt.title('Extracted IR Spectra')
        plt.legend()
        
        #Write plot to file, ensuring that it doesn't overwrite a previous version
        plt_file = os.path.join(directory, f'IR_spectra_plot_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}.png')

        i = 2 
        while os.path.isfile(plt_file):
            plt_file = os.path.join(directory, f'IR_spectra_plot_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}_v{i}.png')
            i += 1

        #save the plot and show it to the user
        plt.savefig(plt_file)
        plt.show()

    ftime = time.time()
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} IR spectra from {len(files)} .out files have been processed in {np.round(ftime - stime, 2)} seconds.')

if __name__ == "__main__":
    # Specify the directory containing .out files
    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\MobCal-MPI\MobCal-MPI\Manual\Appendix_A\DFT'

    # Specify the FWHM for broadening
    fwhm_value = 8

    # Specify the lower and upper bounds
    lower_bound_value = 0 #must be an integer
    upper_bound_value = 4000 #must be an integer
    step_size = 2 # must 

    #harmonic scaling factor
    scale_freq = 0.9875

    #option to normalize all calculated spectra to a maximum intensity of 1
    normalization = True

    # Process the .out files and write the results to an Excel file
    extract_IR_spectra(directory, fwhm_value, lower_bound_value, upper_bound_value, step_size, normalization, scale_freq)
