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

def extract_IR_spectra(directory, fwhm, lower_bound, upper_bound, step_size, normalization, plotting, scale_freq = 1.):

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
            
            except (IndexError, ValueError, TypeError) as e:
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} IR data from {os.path.basename(file)} is in a different format than expect. An error was encountered when trying to parse the freq and intensity: {e}.')
                print(f'{datetime.now().strftime("[ %H:%M:%S ]")} The first and second lines of the problematic IR block are:\n{data_lines[0]}\n{data_lines[1]}\nProcessing of this file will be skipped.')
                return None, None
            
    def broaden_spectrum(file, freq, intensity, fwhm, lower_bound, upper_bound, step_size, normalization):
        '''Broadens each absoprtion in an IR stick strectrum with a user specified Gaussian FWHM.'''
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

        #normalize spectrum if requested
        if normalization:
            broadened_spectrum /= max_value 

        return freq_grid, broadened_spectrum
       
    '''Main workflow'''
    
    #get list of filenames, and check if the directory does not contain any .out files
    filenames = [x for x in os.listdir(directory) if x.lower().endswith('.out') and '_atom' not in x.lower()]

    if len(filenames) == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} There are no .out files in the provided directory.')
        return

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Starting extraction of IR spectra from {len(filenames)} ORCA .out files...')
    QApplication.processEvents()

    stime = time.time()

    #Create a Pandas DataFrame to store the results
    df = pd.DataFrame()

    for file in filenames:

        #Extract IR spectrum from the file - extract_ir_freqs returns None for freq, intensity if an .out file does not contain any IR freq data
        freq, intensity = extract_ir_freqs(os.path.join(directory, file),scale_freq)

        #Broaden the spectrum only if IR data was sucessfuly extracted from the file
        if freq is not None and intensity is not None:
            freq_grid, broadened_spectrum = broaden_spectrum(file, freq, intensity, fwhm, lower_bound, upper_bound, step_size, normalization)

            #Add the broadened spectrum to the DataFrame with the file title as the column header - broaden_spectrum returns None for freq_grid, broadened_spectrum if the max intensity of all IR freqs is zero. 
            if freq_grid is not None and broadened_spectrum is not None:
                df[os.path.basename(file)] = broadened_spectrum

        #If IR data cannot be found and/or the maximum intensity of the IR freqs is zero, skip to the next file
        else:
            continue

    #Add the wavenumber column to the DataFrame
    df['Wavenumber'] = freq_grid

    #Reorder columns
    df = df[['Wavenumber'] + [col for col in df.columns if col != 'Wavenumber']]

    #Write the DataFrame to an Excel file, unsuring that it does not overwrite previous data
    if normalization:
        output_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}cm_Norm.xlsx')

    else:
        output_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}cm.xlsx')

    i = 2 
    while os.path.isfile(output_file):
        if normalization:
            output_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}cm_Norm_v{i}.xlsx')
        else:
            output_file = os.path.join(directory, f'output_spectrum_scaled_{str(scale_freq).replace(".", "-")}_fwhm_{fwhm}cm_v{i}.xlsx')
        i += 1

    #write data to file
    df.to_excel(output_file, index=False)

    ftime = time.time()
    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} IR spectra from {len(filenames)} .out files have been processed in {np.round(ftime - stime, 2)} seconds.')

    #Plotting
    if plotting:

        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Extracted IR spectra will be plotted momentarily...')
        plt.figure(figsize=(16, 10))

        #initialize list to store filenames for generating plot legend
        labels = []

        #Loop through each spectrum in the DataFrame (excluding 'Wavenumber')
        for column in df.columns:
            if column != 'Wavenumber':
                plt.plot(df['Wavenumber'], df[column], label=column)
                
                #Check if the label length exceeds a maximum length so that the legend doesn't get hella cluttered
                if len(column) > 20:
                    match = re.search(r"(\d+)(?!.*\d)", column) #regex to look for the last number in a filename
                    basename = column.split('_')[0]
                    
                    if match and basename:
                        new_label = f'{basename}_{match.group(1)}'  #Use first text preceeded by an underscore and the last number as the legend entry

                    elif match and not basename:
                        new_label = f'File_{match.group(1)}'  #If no underscores in fileanme, use the 'File' prefix and the last number

                    else:
                        new_label = column[(int(-1 * np.ceil(len(column) * 0.25))):-4] #Take first 25% of characters in the filename and remove the file extension [:-4]... users can adjust this number if they don't like it
                
                else:
                    new_label = column
                
                labels.append(new_label)

        #Customizing the plot
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel('Intensity')
        plt.title('Extracted IR Spectra')
        
        #bbox_to_anchor to place the legend outside, adjusting 'loc' and 'bbox_to_anchor' as needed
        plt.legend(labels=labels, loc='upper left', bbox_to_anchor=(1, 1))
        
        #Write plot to file with the same basename as the corresponding .xlsx
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

#external testing
if __name__ == "__main__":
    #Specify the directory containing .out files
    directory = r'D:\OneDrive\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Sample_Files\T9_ORCA_Vib_Spectrum_Analyzer'

    #Specify the FWHM for broadening
    fwhm_value = 8

    #Specify the lower and upper bounds
    lower_bound_value = 500 #must be an integer
    upper_bound_value = 3600 #must be an integer
    step_size = 1 #must be an integer

    #harmonic scaling factor
    scale_freq = 0.9875

    #option to normalize all calculated spectra to a maximum intensity of 1
    normalization = False

    #option to plot spectra
    plotting = False

    #Process the .out files and write the results to an Excel file
    extract_IR_spectra(directory, fwhm_value, lower_bound_value, upper_bound_value, step_size, normalization, plotting, scale_freq)
