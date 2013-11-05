# -*- coding: utf-8 -*-
"""
Created on 2013-07-12

@author: Will
"""

from acwires import HWire, VWire, NWire, HThinWire, VThinWire

#from AtomChip import *

n=4 #number of subwires

I_central = -2.5
I_GU1 = 0
I_GU2 = 0
I_GU3 = 0.5
I_GU4 = 0
I_GU5 = 0
I_GL1 = 0
I_GL2 = 0
I_GL3 = 0.5
I_GL4 = 0
I_GL5 = 0

I_XBias = 3
I_YBias = 0
I_ZBias = 0

I_Macro_Bias_1 = 21
I_Macro_Bias_2 = 21
I_Macro_Central = 0
I_Macro_Axial_1 = 0
I_Macro_Axial_2 = 0
I_Macro_Dimple = 0

# Field Calibration
# Taken from ACM Science Chamber Field Calibrations on wiki
xbiascal = 1054e-7 #T/A
xgradcal = 0.06211e-2 #T/Am
ybiascal = 223.67 #T/A
zbiascal = 1039.7 #T/A

# Bias fields
B_xbias = I_XBias * xbiascal # T
B_ybias = I_YBias * ybiascal # T
B_zbias = I_ZBias * zbiascal # T

# Height Targeting bias selection
# height_target = 500e-6
# B_ybias = 2e-3*I_central*1e-4 / height_target
print("By Bias: %2.2f G"%(B_ybias*1e4))

# Define chip geometry

mww = 100e-6
mwh = 5e-6
mwz = 0

macro_left_H = -2.5e-2
macro_left_ax1 = -0.685e-2
macro_left_ax2 = 0.565e-2
macro_left_dimple = -0.0595e-2
macro_length_H = 5e-2
macro_length_V = 0.119e-2

macro_bottom_H = -0.61e-2
macro_bottom_V = -0.731e-2
macro_height_H = 0.09e-2
macro_height_V = 0.1e-2

macro_low_bias1 = 0.17e-2
macro_low_bias2 = -0.33e-2
macro_low_central = -0.0785e-2
macro_low_V = -2.5e-2
macro_width_H = 0.157e-2
macro_width_V = 3e-2

## Wire definition

# Horizontal Wires
hwires = []
vwires = []
nwires = []
#(name, length, width, height, current, xl, y0, z0, subwires = 1):
hwires.append(HWire('WG',               # name
                    13.6e-3,            # length / m
                    mww,                # width / m
                    mwh,                # height / m
                    I_central,          # current / A
                     -6.8e-3,          # xl / m
                     0,            # y0 / m
                     mwz,
                     n))              # z0 / m
                     
vwires.append(VWire('GU1',               # name
                    3.3e-3,              # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GU1,               # current / A
                     -3e-3,              # x0 / m
                     0,              # yd / m
                     mwz,
                     1))               # z0 / m
                     
vwires.append(VWire('GU2',               # name
                    3.8e-3,              # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GU2,               # current / A
                     -1e-3,              # x0 / m
                     0,              # yd / m
                     mwz))               # z0 / m
                     
vwires.append(VWire('GU3',               # name
                    4.4e-3,              # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GU3,               # current / A
                     0e-6,             # x0 / m
                     0,              # yd / m
                     mwz,
                     n))               # z0 / m
                     
vwires.append(VWire('GU4',               # name
                    4.4e-3,              # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GU4,               # current / A
                     1e-3,             # x0 / m
                     0,              # yd / m
                     mwz))               # z0 / m
                     
vwires.append(VWire('GU5',               # name
                    3.8e-3,              # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GU5,               # current / A
                     3e-3,             # x0 / m
                     0,              # yd / m
                     mwz))               # z0 / m
                     
vwires.append(VWire('GL1',               # name
                    5e-3,                # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GL1,               # current / A
                     -3e-3,              # x0 / m
                     -5e-3,              # yd / m
                     mwz))               # z0 / m
                     
vwires.append(VWire('GL2',               # name
                    5e-3,                # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GL2,               # current / A
                     -1e-3,              # x0 / m
                     -5e-3,              # yd / m
                     mwz))               # z0 / m
                     
vwires.append(VWire('GL3',               # name
                    5e-3,                # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GL3,               # current / A
                     0,                  # x0 / m
                     -5e-3,              # yd / m
                     mwz,
                     n))               # z0 / m
                     
vwires.append(VWire('GL4',               # name
                    5e-3,                # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GL4,               # current / A
                     1e-3,             # x0 / m
                     -5e-3,              # yd / m
                     mwz))               # z0 / m
                     
vwires.append(VWire('GL5',               # name
                    5e-3,                # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GL5,               # current / A
                     3e-3,             # x0 / m
                     -5e-3,              # yd / m
                     mwz))               # z0 / m
                     
hwires.append(HWire('MacroBias1',               # name
                    macro_length_H,            # length / m
                    macro_width_H,                # width / m
                    macro_bottom_H,                # height / m
                    I_Macro_Bias_1,          # current / A
                    macro_left_H,          # xl / m
                    macro_low_bias1,            # y0 / m
                    macro_bottom_H,
                    n))              # z0 / m

hwires.append(HWire('MacroBias2',               # name
                    macro_length_H,            # length / m
                    macro_width_H,                # width / m
                    macro_bottom_H,                # height / m
                    I_Macro_Bias_2,          # current / A
                    macro_left_H,          # xl / m
                    macro_low_bias2,            # y0 / m
                    macro_bottom_H,
                    n))              # z0 / m

hwires.append(HWire('MacroCentral',               # name
                    macro_length_H,            # length / m
                    macro_width_H,                # width / m
                    macro_bottom_H,                # height / m
                    I_Macro_Central,          # current / A
                    macro_left_H,          # xl / m
                    macro_low_central,            # y0 / m
                    macro_bottom_H))              # z0 / m                    
                     
vwires.append(VWire('MacroAxial1',               # name
                    macro_length_V,                # length / m
                    macro_width_V,                 # width / m
                    macro_height_V,                 # height / m
                    I_Macro_Axial_1,               # current / A
                    macro_left_ax1,                  # x0 / m
                    macro_low_V,              # yd / m
                    macro_bottom_V))               # z0 / m
                    
vwires.append(VWire('MacroAxial2',               # name
                    macro_length_V,                # length / m
                    macro_width_V,                 # width / m
                    macro_height_V,                 # height / m
                    I_Macro_Axial_2,               # current / A
                    macro_left_ax2,                  # x0 / m
                    macro_low_V,              # yd / m
                    macro_bottom_V))               # z0 / m

vwires.append(VWire('MacroDimple',               # name
                    macro_length_V,                # length / m
                    macro_width_V,                 # width / m
                    macro_height_V,                 # height / m
                    I_Macro_Dimple,               # current / A
                    macro_left_dimple,                  # x0 / m
                    macro_low_V,              # yd / m
                    macro_bottom_V))               # z0 / m                    
                     
allwires = hwires + vwires + nwires

atomchip_v4 = {'wirespecs' : allwires,
               'B_xbias' : B_xbias,
               'B_ybias' : B_ybias,
               'B_zbias' : B_zbias}
# 
if __name__ == '__main__':
    import bfsimulator
    b_f_sim = bfsimulator.BFieldSimulator()
    b_f_sim.set_chip(atomchip_v4)
    b_f_sim.calc_trap_height()
    b_f_sim.plot_z()
#    b_f_sim.calc_xy()
#    b_f_sim.plot_xy()
    sim_results = b_f_sim.find_trap_freq()
    print 'x_trap : %2.0f um \ny_trap : %2.0f um \nz_trap : %2.0f um'%(b_f_sim.x_trap*1e6, b_f_sim.y_trap*1e6, b_f_sim.z_trap*1e6)
    print 'f_long : %2.0f Hz \nf_trans : %2.0f Hz \nf_z : %2.0f Hz'%(sim_results['f_long'], sim_results['f_trans'], sim_results['f_z'])
#    b_f_sim.set_chip(atomchip_v3)
    b_f_sim.plot_xy()