import os
import numpy as np
import matplotlib.pyplot as plt
import csv


def LED_Analysis(fragment1, fragment2, parent, output_csv):

    # Constants
    J2Eh = 1 / (4.3597447222071e-18)  # Joule to Hartree conversion
    N_Av = 6.02214076e23  # Avogadro's number

    def parse_Fi_file(file):
        with open(file) as f:
            lines = f.readlines()

        E_tot = E_C_tot = E_C_T = 0.0
        for line in lines:
            if line.startswith('Triples Correction (T)'):
                E_C_T = float(line.split()[-1])
            elif line.startswith('Final correlation energy '):
                E_C_tot = float(line.split()[-1])
            elif line.startswith('FINAL SINGLE POINT ENERGY'):
                E_tot = float(line.split()[-1])

        E_C_CCSD = E_C_tot - E_C_T
        E_ref = E_tot - E_C_tot
        return E_ref, E_C_CCSD, E_C_T

    def parse_LED_file(file):
        with open(file) as f:
            lines = f.readlines()

        E_tot = E_elstat = E_exch = E_disp = E_nondisp = E_XvY_ref = E_XvY_C = 0.0
        for line in lines:
            if line.startswith('Total energy                                ='):
                E_tot = float(line.split()[-1])
            if line.startswith('Total energy                                ='):
                E_tot = float(line.split()[-1])
            if line.startswith('Intra REF. energy'):
                E_XvY_ref = float(line.split()[-1])
            if line.startswith('Intra Correlation energy'):
                E_XvY_C = float(line.split()[-1])
            if line.startswith('Electrostatics (REF.)'):
                E_elstat += float(line.split()[-1])
            if line.startswith('Exchange (REF.)'):
                E_exch += float(line.split()[-1])
            if line.startswith('Dispersion (strong pairs)'):
                E_disp += float(line.split()[-1])
            if line.startswith('Dispersion (weak pairs)'):
                E_disp += float(line.split()[-1])
            if line.startswith('Non dispersion (strong pairs)'):
                E_nondisp += float(line.split()[-1])
            if line.startswith('Non dispersion (weak pairs)'):
                E_nondisp += float(line.split()[-1])

        return E_tot, E_XvY_ref, E_XvY_C, E_elstat, E_exch, E_disp, E_nondisp

    def print_summary(E_int, dE_elprep_ref, E_elstat, E_exch, E_disp, dE_nondisp_CCSD, dE_C_T, filename='LED_summary.csv'):
        summary = [
            ('E_int', E_int),
            ('E_int (check)', (dE_elprep_ref + E_elstat + E_exch + E_disp + dE_nondisp_CCSD + dE_C_T) / J2Eh * N_Av * 1e-3),
            ('E_elprep(HF)', dE_elprep_ref / J2Eh * N_Av * 1e-3),
            ('E_elstat(HF)', E_elstat / J2Eh * N_Av * 1e-3),
            ('E_exch(HF)', E_exch / J2Eh * N_Av * 1e-3),
            ('E_disp(C)', E_disp / J2Eh * N_Av * 1e-3),
            ('E_nondisp(C)', dE_nondisp_CCSD / J2Eh * N_Av * 1e-3),
            ('E_trip', dE_C_T / J2Eh * N_Av * 1e-3)
        ]

        print('\n== LED summary ==')
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Description', 'Value (kJ/mol)'])  # Header row

            for desc, value in summary:
                print(f'{desc:15s} = {value:+.8f} kJ/mol')
                writer.writerow([desc, f'{value:+.8f}'])

        print(f"\nSummary written to {filename}")

    def plot_summary(Es, label):
        x1 = np.arange(1,8)
        plt.figure()
        off = np.max(Es) * 0.05
        plt.title('LED for ' + label)

        plt.bar(x1, Es, color='C0')
        for i, E in enumerate(Es):
            plt.text(x1[i], E + off * np.sign(E), '%+.1f' % E, ha='center', va='center')

        XXTicks = [np.arange(1,8), [r'$E_{elprep}$', r'$E_{elstat}$', r'$E_{exch}$', r'$E_{disp}$', r'$E_{nondisp}$', r'$E_{trip}$', r'$E_{tot}$']]

        plt.hlines(0.0, 0, 8, colors='k', alpha=0.5, lw=0.5)
        plt.xlim(0, 8)
        plt.ylabel(r'interaction energy in kJ/mol')
        plt.xticks(*XXTicks)
        plt.tight_layout()
        plt.savefig(label + '_LED.png', format='png', dpi=300, bbox_inches='tight')
        plt.show()

    def run_LED_analysis(F1_file, F2_file, LED_file, output_file):
        F1_E_ref, F1_E_C_CCSD, F1_E_C_T = parse_Fi_file(F1_file)
        F2_E_ref, F2_E_C_CCSD, F2_E_C_T = parse_Fi_file(F2_file)
        Ad_E_ref, Ad_E_C_CCSD, Ad_E_C_T = parse_Fi_file(LED_file)
        E_tot, E_XvY_ref, E_XvY_C, E_elstat, E_exch, E_disp, E_nondisp = parse_LED_file(LED_file)

        # Calculations
        F1_E_tot = F1_E_ref + F1_E_C_CCSD + F1_E_C_T
        F2_E_tot = F2_E_ref + F2_E_C_CCSD + F2_E_C_T
        Ad_E_tot = Ad_E_ref + Ad_E_C_CCSD + Ad_E_C_T

        E_int = E_tot - (F1_E_tot + F2_E_tot)
        dE_elprep_ref = E_XvY_ref - (F1_E_ref + F2_E_ref)
        dE_nondisp_CCSD = E_nondisp - (F1_E_C_CCSD + F2_E_C_CCSD)
        dE_C_T = Ad_E_C_T - (F1_E_C_T + F2_E_C_T)

        # Print Summary
        print_summary(E_int, dE_elprep_ref, E_elstat, E_exch, E_disp, dE_nondisp_CCSD, dE_C_T, output_file)

        # Plot Summary
        Es = np.array([dE_elprep_ref, E_elstat, E_exch, E_disp, dE_nondisp_CCSD, dE_C_T, E_int]) / J2Eh * N_Av * 1e-3  # Convert to kJ/mol
        label = LED_file.rsplit('/', 1)[-1][:-8]
        plot_summary(Es, label)

    # Main execution
    run_LED_analysis(fragment1, fragment2, parent, output_csv)

if __name__ == '__main__':
    directory = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Python\Special_Analyses'
    parent_LED = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Python\Special_Analyses\HMX_NO3_C0_LED.out'
    fragment1_Ehigh = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Python\Special_Analyses\NO3_Ehigh.out'
    fragment2_Ehigh = r'C:\Users\Chris\OneDrive - University of Waterloo\Waterloo\GitHub\ORCA_Analysis_GUI\Python\Special_Analyses\HMX_neut_C0_Ehigh.out'
    output_csv = os.path.join(directory, 'LED_analysis.csv')

    LED_Analysis(fragment1_Ehigh, fragment2_Ehigh, parent_LED, output_csv)
