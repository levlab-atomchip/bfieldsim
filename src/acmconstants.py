# -*- coding: utf-8 -*-
"""
Created on Thu May 02 21:43:22 2013

@author: Will
"""

from math import pi
from Window import window

# Physical Constants, from Steck Rb87 Alkali D Line Data
C = 2.99792e8 # m/s (exact)
ALPHA = 7.29735e-3 # from wikipedia
HBAR = 1.05457e-34 # J*sec
MU_B = 9.274e-24 # J/T
A_0 = 5.29e-11 #m, Bohr radius
MU_0 = 4*pi*1e-7
K_B = 1.38e-23 #J/K, boltzmann constant
E = 1.602e-19 # C, elementary charge


## Imaging Transition Properties
#THESE ARE WRONG!!!
#G1 = 1
#D12 = 1.73135e-29 #C*m
#D12 = 2.042*A_0
D_eff = 2.989*E*A_0
#G2 = 1
OMEGA_RES = 2*pi*384.230e12 #Hz
LINEWIDTH_RES = 2*pi*6.0666e6 #Hz
WAVELENGTH_RES = 780e-9 #m
I_SAT = 16.7 #W/m^2
SIGMA_0 = 2.907e-13 #m^2, resonant x section

# Rb87 Properties
A = 5.2e-9 #m, scattering length
M = 87 * 1.7e-27 #kg
G = (4 * pi * HBAR**2 * A) / M

# Dy Properties
#A = 100 * A_0
#M = 162.5 * 1.7e-27
#MU_B = 10* 9.274e-24

CLOUD_THICKNESS = window.cell_size
BG_NOISE_CURRENT = 30 #e/p/s

# Camera Properties
# PIXIS 1024BR
#PIXEL_SIZE = 13e-6 #m
PIXEL_SIZE = 100e-6 / 1024
NUM_PIXELS = 1024
DARK_CURRENT = .07 #e/p/s