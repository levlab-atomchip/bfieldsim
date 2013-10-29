# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 20:28:42 2013

@author: Will
"""

from biot_savart import *
from CurrentSlab import *
from math import pi
from Window import Window
import numpy as np
from acmconstants import MU_0
import matplotlib.pyplot as plt

rez = 1e-2
th = 1e-3

y_0 = 0
z_0 = 0.1
slab = Slab(1, 10, th, rez)

wire_center = 0
wire_width = 4*rez
wire_J = 1e-2

def wire_wrapper(x_0, w, J):
    return lambda xx, yy: J_perp_wire(xx, yy, x_0, w, J)

def J_perp_wire(xx, yy, x_0, w, J):
    zz1 = xx < (x_0 + 0.5*w)
    zz2 = xx >= (x_0 - 0.5*w)
    return J * np.logical_and(zz1, zz2)
    
def theta_perp_wire(xx, yy):
    zz = 0.5 * pi * np.ones(xx.shape)
    return zz
    
window = Window(-0.1, 0.1, 100)    

current_slab = CurrentSlab(slab, wire_wrapper(wire_center, wire_width, wire_J), theta_perp_wire)

B = biot_savart(window, y_0, z_0, current_slab)

B.plot_mag()

prefactor = MU_0 * (wire_width * th * wire_J)  * 0.5 / pi
r = np.sqrt(window.window**2 + y_0**2 + z_0**2)
B_inf = prefactor / r
plt.plot(window.window, B_inf)
plt.show()
Bratio = B_inf / B.get_mag()
plt.plot(window.window, Bratio - 1)
current_slab.plot_J()