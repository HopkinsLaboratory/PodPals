import numpy as np

## fundamental constants
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

## fundamental constants in other units
kB_Eh = kB*J2Eh # hartree/kelvin
h_cm = h_SI/(h_SI*c_SI*100.) # cm**-1 * s, Planck constant in wavenumber*s


## fundamental constants
print(f"pi4eps: {pi4eps}")
print(f"kB: {kB}")
print(f"N_Av: {N_Av}")
print(f"h_SI: {h_SI}")
print(f"h_bar: {h_bar}")
print(f"c_SI: {c_SI}")
print(f"elC: {elC}")
print(f"m_el: {m_el}")
print(f"r_Bohr: {r_Bohr}")

## unit conversions
print(f"bohr2m: {bohr2m}")
print(f"bohr2A: {bohr2A}")
print(f"cm2in: {cm2in}")
print(f"J2Eh: {J2Eh}")
print(f"amu2kg: {amu2kg}")
print(f"D2Cm: {D2Cm}")
print(f"eigval2s2: {eigval2s2}")

## fundamental constants in other units
print(f"kB_Eh: {kB_Eh}")
print(f"h_cm: {h_cm}")


