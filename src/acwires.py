# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:54:24 2013

@author: Will
"""

from math import pi, sqrt
import numpy as np
from acmconstants import MU_0

class Wire():
    def __init__(self, name, length, width, height, current, subwires):
        self.name = name
        self.length = length
        self.width = width
        self.height = height
        self.current = current
        self.subwires = subwires
        
class HWire(Wire):
    def __init__(self, name, length, width, height, current, xl, y0, z0, subwires = 1):
        Wire.__init__(self, name, length, width, height, current, subwires)
        self.xl = xl
        self.y0 = y0
        self.z0 = z0
        
    def bfieldcalc(self,x,y,z):
        #loops over subwires, calculating magnetic field at point (x,y,z)
#        print 'executing bfieldcalc'
        XL = self.xl
        XR = self.xl + self.length
        suby = np.linspace(self.y0, self.y0 + self.width, self.subwires)
        subz = np.linspace(self.z0, self.z0 + self.height, self.subwires)
        nn = self.subwires**2
        const_G=MU_0*self.current/(4*pi*nn)
        for i in range(self.subwires): #loop over y
            for j in range(self.subwires): #loop over z
                beta = (z-subz[j])**2 + (y-suby[i])**2
                B_G=const_G*((x-XL)/(beta*sqrt(beta+
                    (x-XL)**2))-(x-XR)/(beta*sqrt(beta+(x-XR)**2)))
#                print 'B_G = %f'%B_G
                B_Gy=B_G*(subz[j]-z)
                B_Gz=B_G*(y-suby[i])
        return np.array((0, B_Gy, B_Gz))
        
class HThinWire(HWire):        
    def bfieldcalc(self,x,y,z):
        '''Calculate field analytically for finite length infinitely thin wire'''
        XL = self.xl
        XR = self.xl + self.length
        const_G=MU_0*self.current/(4*pi)
        beta = (z-self.z0)**2 + (y-self.y0)**2
        B_G=const_G*((x-XL)/(beta*sqrt(beta+
                    (x-XL)**2))-(x-XR)/(beta*sqrt(beta+(x-XR)**2)))
#                print 'B_G = %f'%B_G
        B_Gy=B_G*(self.z0-z)
        B_Gz=B_G*(y-self.y0)
        return np.array((0, B_Gy, B_Gz))
        
class HFlatWire(HWire):
    def bfieldcalc(self, x, y, z):
        ''' Calculate field analytically for finite length,
        finite width, infinitely flat wire'''
        
        
                
class VWire(Wire):
    def __init__(self, name, length, width, height, current, x0, yd, z0, subwires = 1):
        Wire.__init__(self, name, length, width, height, current, subwires)
        self.x0 = x0
        self.yd = yd
        self.z0 = z0
        
    def bfieldcalc(self,x,y,z):
        #loops over subwires, calculating magnetic field at point (x,y,z)
        YD = self.yd
        YU = self.yd + self.length
        subx = np.linspace(self.x0, self.x0 + self.width, self.subwires)
        subz = np.linspace(self.z0, self.z0 + self.height, self.subwires)
        nn = self.subwires**2
        const_G=MU_0*self.current/(4*pi*nn)
        for i in range(self.subwires): #loop over x
            for j in range(self.subwires): #loop over z
                beta = (z-subz[j])**2 + (x-subx[i])**2
                B_G=const_G*((y - YD)/(beta*sqrt(beta+
                    (y - YD)**2))-(y - YU)/(beta*sqrt(beta+(y - YU)**2)))
                B_Gx=B_G*(z - subz[j])
                B_Gz=B_G*(subx[i] - x)
        return np.array((B_Gx, 0, B_Gz))
        
class VThinWire(VWire):
    def bfieldcalc(self,x,y,z):
        #loops over subwires, calculating magnetic field at point (x,y,z)
        YD = self.yd
        YU = self.yd + self.length
        const_G=MU_0*self.current/(4*pi)
        beta = (z-self.z0)**2 + (x-self.x0)**2
        B_G=const_G*((y - YD)/(beta*sqrt(beta+
                    (y - YD)**2))-(y - YU)/(beta*sqrt(beta+(y - YU)**2)))
        B_Gx=B_G*(z - self.z0)
        B_Gz=B_G*(self.x0 - x)
        return np.array((B_Gx, 0, B_Gz))
        
        
class NWire(Wire):
    def __init__(self, name, length, width, height, current, x0, y0, zd, subwires = 1):
        Wire.__init__(self, name, length, width, height, current, subwires)
        self.x0 = x0
        self.y0 = y0
        self.zd = zd
        
    def bfieldcalc(self,x,y,z):
        #loops over subwires, calculating magnetic field at point (x,y,z)
        ZD = self.zd
        ZU = self.zd + self.length
        subx = np.linspace(self.x0, self.x0 + self.width, self.subwires)
        suby = np.linspace(self.y0, self.y0 + self.height, self.subwires)
        nn = self.subwires**2
        const_G=MU_0*self.current/(4*pi*nn)
        for i in range(self.subwires): #loop over x
            for j in range(self.subwires): #loop over y
                beta = (x-subx[i])**2 + (y-suby[j])**2
                B_G=const_G*((z - ZD)/(beta*sqrt(beta+
                    (z - ZD)**2))-(z - ZU)/(beta*sqrt(beta+(z - ZU)**2)))
                B_Gx=B_G*(suby[j] - y)
                B_Gy=B_G*(x - subx[i])
        return np.array((B_Gx, B_Gy, 0))