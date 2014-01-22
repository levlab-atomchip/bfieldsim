# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:45:57 2013

@author: Will
"""

## Trap parameters
# The script needs:
#       -the bias field (B_ybias) in Gauss
#       -the central wire length (l) in microns
#       -the trap height (h) in microns
#       -the wire current (I) in amps
#
# Formulae are provided for deriving the necessary bias field to realize a
# particular trap height. Trap length must be provided, plus any two of
# (B_ybias, I, and h).


# Axes are defined so that z is up (out of the chip plane), x is along the
# central axis of the z trap ( = ODT axis) and y follows the right-hand
# rule.

#The z ground plane is the top of the MACOR surface.


#physical constants
mu_B = 9.274e-24 #J/Tesla
m_Rb87 = 87*1.66053886e-27 #kg
hbar = 1.05457e-34 #J*s
g = 9.813 #m/s^2
k_B = 1.38e-23 #J/K
rho_Cu = 16.78e-9 #ohms * meters
rho_Ag = 15.87e-9 #ohm * meter

#Wire Geometry Parameters
l_trap = 13.5e3#microns central and y-bias wire lengths
l_ax = 10e3 #microns longitudinal trapping wire lengths
l_med = 2e3
l_micro = 5e3
l_arms = 4.5e3
a = 2.5e3 #microns wire spacing, between central wire and integrated bias wires
a_ax = 12.5e3 #microns wire spacing, between submerged longitudinal bias wires
#sub = 1.7e3 #microns wire submersion, for 
#longitudinal bias wires from ground plane to middle
w_rad = 1587 #wire width, microns; 1/16"
w_ax = 1191 #microns; 3/64"
w_dimple = 1191 #microns; 3/64"
w_micro = 100
w_micro_lead = 500
w_micro_thlead = 1000
h_rad = 1000
h_ax = 1000 #microns, wire height
h_dimple = 1000
h_micro = 5
sp = 10
sub = 1.5*h_rad + 200 #microns assumes 200 um vertical spacing...
chip_height = 500
n = 10 #number of subwires per side, used to simulate wide wires.

#Lead Geometry Parameters, meters
mount_width = 30e-3
mount_length = 50e-3

#Trap parameters
h = 1000#microns
I_in = 25.4#amp (a guess)
beta = .5*((a/(h+h_rad/2))**2 + 1)
I_out = .80*beta * I_in #Use a fudge factor to tune trap height
#I_out = 37.4
I_ax =11.7 #amps, guess
I_dimple =0
I_micro_trap = 0
I_micro_arm = 0
I_micro_dimple = 0
# I_out = 1.02*[((h_rad/2 + h)^2 + a^2)/(2*h + h_rad)]*[I_in/(h + h_rad/2) 
#+ (I_micro_trap)/(h - chip_height)] #Use when microwires are included

#These should be altered to match the output of the script, and iterated
#until input and output converge.
x_trap = -1.25e-6 #meters
y_trap = 1.25e-6 #meters

#x_trap = 13.75e-6 #meters
#y_trap = -2.5e-6 #meters

# Set bias field to ensure trap height
# Bias fields (in Gauss)
B_xbias=0
B_zbias=0
B_ybias = 2.06 #integrated bias
#B_ybias = (2000*I_micro_trap/(h-chip_height)) #accurate for h << l
#B_ybias = 1000*I*(l/h)*(h^2 +.25*l^2)^(-.5) #more accurate

# set the trap min for calculating the slice in 3D
#If height is a fixed parameter:
z=h
#If B_zbias is a fixed parameter:
#z = abs(2000*I/B_ybias)
#alpha = 1000*I/B_ybias
#z = (l/sqrt(8))*((1+64*alpha^2/l^2)^.5 -1)^.5 # more accurate

## Unit Conversions (MKS from here on out!)
# convert microns to meters
a = a*1e-6
a_ax = a_ax*1e-6
l_trap = l_trap*1e-6
l_ax = l_ax*1e-6
l_med = l_med*1e-6
l_micro = l_micro*1e-6
l_arms = l_arms*1e-6
w_rad = w_rad*1e-6
w_ax = w_ax*1e-6
w_dimple = w_dimple*1e-6
w_micro = w_micro*1e-6
h_rad = h_rad*1e-6
h_ax = h_ax*1e-6
h_dimple = h_dimple*1e-6
h_micro = h_micro*1e-6
sp = sp*1e-6
sub = sub*1e-6
chip_height = chip_height*1e-6
z = z*1e-6


# convert fields to Tesla
B_xbias=B_xbias*1e-4
B_ybias=B_ybias*1e-4
B_zbias=B_zbias*1e-4
