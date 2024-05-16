import numpy as np

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

'''#Functions for unit conversions'''

def wavelength_to_wavelength(energy):
    '''Wavelength (nm) to Wavelength (nm); same input - redunency needed for the way the dictionary is set up that calls the unit conversion functions'''
    return energy

def wavelength_to_wavenumber(energy):
    '''Convert wavelength (nm) to wavenumber (cm^-1)
    1E-7 cm per nm, then inversed'''

    return 1 / (energy * 1E-7) 

def wavelength_to_eV(energy):
    '''Convert wavelength (nm) to eV
    Wavelength (nm) to Joules via (h*c)/(lambda * 1E-9), then J to eV via multiplication by electric constant (eV2J)'''

    return ((c['h_SI'] * c['c_SI']) / (energy * 1.E-9)) / c['eV2J'] 

def wavenumber_to_wavelength(energy):
    '''Convert wavenumber (cm^-1) to wavelength (nm)
    1E-7 cm per nm, then inversed'''
    
    return 1 / (energy * 1E-7) # 1E-7 cm per nm

def wavenumber_to_wavenumber(energy):
    '''Wavenumber (cm**-1) to Wavenumber (cm**-1); same input - redunency needed for the way the dictionary is set up that calls the unit conversion functions'''
    
    return energy

def wavenumber_to_eV(energy):
    '''# Convert wavenumber (cm^-1) to eV
    Wavenumber (cm**-1) to Joules via (h*c) * wavenumber * 100cm/m, then J to eV via multiplication by electric constant (eV2J)'''
    
    return (c['h_SI'] * c['c_SI'] * energy * 100.) / c['eV2J']

def eV_to_wavelength(energy):
    '''Convert eV to wavelength (nm)
    First, convert eV to Joules via multiplication by electric constant, then Joules to wavelength via (hc/lambda) * 1E9 nm/m '''
    return ((c['h_SI'] * c['c_SI']) / (energy * c['eV2J'])) * 1.E9

def eV_to_wavenumber(energy):
    '''# Convert eV to wavenumber (cm^-1)
    First, convert eV to Joules via multiplication by electric constant, then Joules to wavenumber (cm**-1) via (E/hc) * (1m / 100cm) '''
    return ((energy * c['eV2J']) / (c['h_SI'] * c['c_SI'])) / 100.

def eV_to_eV(energy):
    return energy


def convert_energy(energy, input_unit, output_unit):
    '''A handy dictionary to call energy conversion functions from. Currently supported units are nm, eV, and cm**-1. 
    Usage: (energy, input_unit, output_unit)
    Energy: Energy to be converted
    input_unit: Unit of energy provided (nm, eV, or cm**-1)
    output_unit: Unit to convert the energy into (nm, eV, or cm**-1)
    
    Returns the conversion function that will do the energy conversion. '''

    conversion_function = {
        ('nm', 'nm'): wavelength_to_wavelength,
        ('nm', 'cm**-1'): wavelength_to_wavenumber,
        ('nm', 'eV'): wavelength_to_eV,
        ('cm**-1', 'nm'): wavenumber_to_wavelength,
        ('cm**-1', 'cm**-1'): wavenumber_to_wavenumber,
        ('cm**-1', 'eV'): wavenumber_to_eV,
        ('eV', 'nm'): eV_to_wavelength,
        ('eV', 'cm**-1'): eV_to_wavenumber,
        ('eV', 'eV'): eV_to_eV
    }
    
    conversion_key = (input_unit, output_unit)
    if conversion_key in conversion_function:
        return conversion_function[conversion_key](energy)
    else:
        raise ValueError(f'Conversion from {input_unit} to {output_unit} is not supported.')

#external testing
if __name__ == "__main__":
    input_unit = 'nm'
    output_unit = 'nm'

    input_energy = 300 #nm
    
    shift_unit = 'eV'
    shift_by = -0.25

    inp_to_shift = convert_energy(input_energy, input_unit, shift_unit)
    print(f'Input to shift unit: {input_energy} {input_unit} was converted to {inp_to_shift} {shift_unit}')

    shifted_inp = inp_to_shift + shift_by
    print(f' {inp_to_shift} {shift_unit} has been shifted to {shifted_inp} {shift_unit}')

    conv_shifted_energy = convert_energy(shifted_inp, shift_unit, output_unit)
    print(f'{input_energy} {input_unit} was shifted by {shift_by} {shift_unit} and coverted to {output_unit}, yielding a value of {conv_shifted_energy} {output_unit}.')

