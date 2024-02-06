import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
matplotlib.use('QtAgg') #Use matplotlib backend that is compatible w/ PyQt6 to prevent competition w/ GUI event loop

from PyQt6.QtWidgets import QApplication

from datetime import datetime
import os

def ORCA_opt_plt(file, save_plot = False):

    def extract_info(lines):

        #empty lists to store data
        energy, e_change, e_threshold = [], [], []
        max_step, max_step_threshold = [], []
        max_grad, max_grad_threshold = [], []
        rms_step, rms_step_threshold = [], []
        rms_grad, rms_grad_threshold = [], []
        e_1st_scf_cycle, change_1st_scf_cycle = [], []

        for index, line in enumerate(lines):
            if line.startswith("FINAL SINGLE POINT ENERGY"):  # ENERGY (1)
                energy.append(float(line.split()[4]))

            elif line.startswith("          Energy change"):  # ENERGY CHANGE (2)
                e_change.append(float(line.split()[2]))
                e_threshold.append(float(line.split()[3]))

            elif line.startswith("          MAX gradient"):  # MAX GRADIENT (3)
                max_grad.append(float(line.split()[2]))
                max_grad_threshold.append(float(line.split()[3]))

            elif line.startswith("          MAX step"):  # MAX STEP (4)
                max_step.append(float(line.split()[2]))
                max_step_threshold.append(float(line.split()[3]))

            elif line.startswith("          RMS gradient"):  # RMS GRADIENT (5)
                rms_grad.append(float(line.split()[2]))
                rms_grad_threshold.append(float(line.split()[3]))

            elif line.startswith("          RMS step"):  # RMS STEP (6)
                rms_step.append(float(line.split()[2]))
                rms_step_threshold.append(float(line.split()[3]))

            elif line.startswith("                         !        ITERATION"):  # 1st SCF cycle (7)
                e_1st_scf_cycle.append(float(lines[index + 1].split()[3]))
                change_1st_scf_cycle.append(float(lines[index + 2].split()[3]))

        return (energy, e_change, e_threshold, max_step, max_step_threshold,
                max_grad, max_grad_threshold, rms_step, rms_step_threshold,
                rms_grad, rms_grad_threshold, e_1st_scf_cycle, change_1st_scf_cycle)

    #If there is only one SCF cycle, enter plot option 1
    def plot_option_1(e_1st_scf_cycle, change_1st_scf_cycle):
        del change_1st_scf_cycle[0]

        plt.figure()

        ax1 = plt.subplot(211)
        plt.title("1st SCF CYCLE Energy")
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax1.plot(range(1, len(e_1st_scf_cycle) + 1), e_1st_scf_cycle, color='red')
        ax1.grid(True)

        ax2 = plt.subplot(212)
        plt.title("1st SCF CYCLE Energy Change")
        ax2.plot(range(1, len(change_1st_scf_cycle) + 1), change_1st_scf_cycle, color='blue')
        ax2.grid(True)

    #if the calculation made it past multiple scf cycles, enter option 2
    def plot_option_2(e_change, e_threshold, energy, max_grad, max_grad_threshold,
                    max_step, max_step_threshold, rms_grad, rms_grad_threshold,
                    rms_step, rms_step_threshold):
        plt.figure(figsize=(20, 10))

        axes = []  # Collect axes for setting grid
        
        ax1 = plt.subplot(321)
        axes.append(ax1)
        plt.title("Energy Change")
        ax1.plot(range(1, len(e_change) + 1), e_change, color='red', label='Energy Change')
        ax1.plot(range(1, len(e_threshold) + 1), e_threshold, color='blue', linestyle='dotted', label='Threshold')
        plt.legend()

        ax2 = plt.subplot(322)
        axes.append(ax2)
        plt.title("SCF Energy")
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))
        ax2.scatter(range(len(energy)), energy, color='red', marker='o')

        ax3 = plt.subplot(323)
        axes.append(ax3)
        plt.title("Max Gradient")
        ax3.plot(range(len(max_grad)), max_grad, color='green', label='Max Gradient')
        ax3.plot(range(len(max_grad_threshold)), max_grad_threshold, color='blue', linestyle='dotted', label='Threshold')
        plt.legend()

        ax4 = plt.subplot(324)
        axes.append(ax4)
        plt.title("Max Step")
        ax4.scatter(range(len(max_step)), max_step, color='green', label='Max Step')
        ax4.plot(range(len(max_step_threshold)), max_step_threshold, color='blue', linestyle='dotted', label='Threshold')
        plt.legend()

        ax5 = plt.subplot(325)
        axes.append(ax5)
        plt.title("RMS Gradient")
        ax5.plot(range(len(rms_grad)), rms_grad, color='black', label='RMS Gradient')
        ax5.plot(range(len(rms_grad_threshold)), rms_grad_threshold, color='blue', linestyle='dotted', label='Threshold')
        plt.legend()

        ax6 = plt.subplot(326)
        axes.append(ax6)
        plt.title("RMS Step")
        ax6.scatter(range(len(rms_step)), rms_step, color='black', label='RMS Step')
        ax6.plot(range(len(rms_step_threshold)), rms_step_threshold, color='blue', linestyle='dotted', label='Threshold')
        plt.legend()

        # Add grid to all axes
        for ax in axes:
            ax.grid(True)

    print(f'{datetime.now().strftime("[ %H:%M:%S ]")} Analyzing the optimization progress in {os.path.basename(file)}. A plot will appear momentarily...')
    QApplication.processEvents()

    try:                        
        with open(file, 'r') as opf:
            lines = opf.readlines()

    except IOError as e:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error was encountered when trying to open {os.path.basename(file)}: {e}.\n')
        return
    
    #extract pertinent information from the ORCA .out file
    try:
        (energy, e_change, e_threshold, max_step, max_step_threshold,
        max_grad, max_grad_threshold, rms_step, rms_step_threshold,
        rms_grad, rms_grad_threshold, e_1st_scf_cycle, change_1st_scf_cycle) = extract_info(lines)
    
    except Exception as e:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error was encountered when trying extracting energies from {os.path.basename(file)}: {e}')
        return

    #if no SCF cycles have been completed, the job probably didn't start correctly. Throw a value error to the user
    if len(energy) == 0:
        print(f'{datetime.now().strftime("[ %H:%M:%S ]")} No SCF cycles have finished in {os.path.basename(file)}! Did it start correctly?')
        return
    
    #if one SCF cycle has been completed, enter plot routine 1 (e.g., a single point energy)    
    if not e_change:
        plot_option_1(e_1st_scf_cycle, change_1st_scf_cycle)
    
    #if multiple scf cycles have been completed, enter plot routine 2 (e.g., an optimization)    
    else:
        plot_option_2(e_change, e_threshold, energy, max_grad, max_grad_threshold,
                      max_step, max_step_threshold, rms_grad, rms_grad_threshold,
                      rms_step, rms_step_threshold)

    plt.tight_layout()
    plt.show()

    if save_plot:
        #Write plot to file, ensuring that it doesn't overwrite a previous version
        plt_file = os.path.join(os.path.dirname(file), f'{str(file).split('.out')[0]}.png')

        i = 2 
        while os.path.isfile(plt_file):
            plt_file = os.path.join(os.path.dirname(file), f'{str(file).split('.out')[0]}_v{1}.png')
            i += 1

        #save the plot and show it to the user
        try:
            plt.savefig(plt_file)
        
        except Exception as e:
            raise Exception(f'{datetime.now().strftime("[ %H:%M:%S ]")} An error occured when trying to save the plot to {plt_file}: {e}.')

#external testing
if __name__ == "__main__":
    input_path = r'E:\Hopkins_Laboratory\Fentanyl_Analogs\Fentanyl\N-Prot\DFT\Fentanyl_NProt_1.out'

