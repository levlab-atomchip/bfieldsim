# -*- coding: utf-8 -*-
"""
Created on 2013-07-12

@author: Will
"""

from acwires import HWire, VWire, NWire, HThinWire, VThinWire

#from AtomChip import *

I_central = 2
I_GU1 = 0
I_GU2 = 0
I_GU3 = 0
I_GU4 = 0
I_GU5 = 2
I_GL1 = 2
I_GL2 = 0
I_GL3 = 0
I_GL4 = 0
I_GL5 = 0

B_xbias = 0 # T
#B_ybias = 20e-4 # T
height_target = 500e-6
B_ybias = 2e-3*I_central*1e-4 / height_target
print B_ybias
B_zbias = 0 # T

mww = 100e-6
mwh = 5e-6
mwz = 0

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
                     mwz))              # z0 / m
                     
vwires.append(VWire('GU1',               # name
                    3.3e-3,              # length / m
                    mww,                 # width / m
                    mwh,                 # height / m
                    I_GU1,               # current / A
                     -3e-3,              # x0 / m
                     0,              # yd / m
                     mwz))               # z0 / m
                     
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
                     mwz))               # z0 / m
                     
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
                     mwz))               # z0 / m
                     
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
                     

allwires = hwires + vwires + nwires

atomchip_v3 = {'wirespecs' : allwires,
               'B_xbias' : B_xbias,
               'B_ybias' : B_ybias,
               'B_zbias' : B_zbias}
# 
if __name__ == '__main__':
    import bfsimulator
    b_f_sim = bfsimulator.BFieldSimulator()
    b_f_sim.set_chip(atomchip_v3)
    b_f_sim.calc_trap_height()
    b_f_sim.plot_z()
#    b_f_sim.calc_xy()
#    b_f_sim.plot_xy()
    sim_results = b_f_sim.find_trap_freq()
    print 'x_trap : %2.0f um \ny_trap : %2.0f um \nz_trap : %2.0f um'%(b_f_sim.x_trap*1e6, b_f_sim.y_trap*1e6, b_f_sim.z_trap*1e6)
    print 'f_long : %2.0f Hz \nf_trans : %2.0f Hz \nf_z : %2.0f Hz'%(sim_results['f_long'], sim_results['f_trans'], sim_results['f_z'])
#    b_f_sim.set_chip(atomchip_v3)
    b_f_sim.plot_xy()