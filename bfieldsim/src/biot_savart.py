# -*- coding: utf-8 -*-
"""
Created on Fri Jun 07 21:15:30 2013

@author: Will
"""

from acmconstants import MU_0
from math import pi
import numpy as np
from BField import BField, BField_2D
from Window import Point

def biot_savart(window, y_0, z_0, current_slab):
    delx = current_slab.get_delx()
    th = current_slab.get_th()
    Y = current_slab.get_Y()
    X = current_slab.get_X()
    dx = current_slab.get_dx()
    dy = current_slab.get_dy()
    J = current_slab.get_J()

    B = BField(window)    
    prefactor = MU_0 * delx * th / (4 * pi)
    r2_yz = (y_0 - Y)**2 + z_0**2
    bx = z_0 * dy
    by = -1*z_0 * dx
    bz1 = (y_0 - Y)*dx
    bz2 = X*dy
    
    for x in window.window:
        bz3 = -1*dy*x
        bz = bz1 + bz2 + bz3
        r2_x = (x - X)**2
        denom = (r2_x + r2_yz)**1.5
        Bx = prefactor * np.sum(J * bx / denom)
        By = prefactor * np.sum(J * by / denom)
        Bz = prefactor * np.sum(J * bz / denom)
        B.set_field(x, [Bx, By, Bz])
    
    return B

def biot_savart_2D(window, z_0, current_slab):
    delx = current_slab.get_delx()
    th = current_slab.get_th()
    Y = current_slab.get_Y()
    X = current_slab.get_X()
    dx = current_slab.get_dx()
    dy = current_slab.get_dy()
    J = current_slab.get_J()

    B = BField_2D(window)    
    prefactor = MU_0 * delx * th / (4 * pi)
    r2_z = z_0**2
    bx = z_0 * dy
    by = -1*z_0 * dx
    bz1 = -1*Y*dx
    bz2 = X*dy

    for x in np.nditer(window.X):
        for y in np.nditer(window.Y):
#            x = window.X
#            y = window.Y
            bz3 = -1*dy*x
            bz4 = dx*y
            bz = bz1 + bz2 + bz3 + bz4
            r2_x = (x - X)**2
            r2_y = (y - Y)**2
            denom = (r2_x + r2_y + r2_z)**1.5
            Bx = prefactor * np.sum(J * bx / denom)
            By = prefactor * np.sum(J * by / denom)
            Bz = prefactor * np.sum(J * bz / denom)
            B.set_field(Point(x,y), [Bx, By, Bz])
    
    return B
    
if __name__ == '__main__':
    import biot_savart_test