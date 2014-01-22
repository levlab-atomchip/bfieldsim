# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:46:12 2013

@author: Will
"""

from acwires import HWire, VWire, NWire, HThinWire, VThinWire

from AtomChip import *

## Wire definition

# Horizontal Wires
hwires = []
vwires = []
nwires = []
#(name, length, width, height, current, xl, y0, z0, subwires = 1):
hwires.append(HWire('Central Test Wire', 1e-2, 100e-6, 5e-6, 1, 
                     -0.5e-2, -50e-6, 0))

vwires.append(VWire('Arm 1', 0.5e-2, 100e-6, 5e-6, 1, -0.5e-2, -0.5e-2, 0))
vwires.append(VWire('Arm 2', 0.5e-2, 100e-6, 5e-6, 1,  0.5e-2, 0, 0))

allwires = hwires + vwires + nwires

# 