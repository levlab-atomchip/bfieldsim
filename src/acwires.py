# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:54:24 2013

@author: Will
"""

from math import pi, sqrt, atan, atanh
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
        
# class HFlatWire(HWire):
    # def bfieldcalc(self, x, y, z):
        # ''' Calculate field analytically for finite length,
        # finite width, infinitely flat wire'''
        
        
                
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
        B_G=const_G*((y - YD)/(beta*sqrt(beta+(y - YD)**2))
                    -(y - YU)/(beta*sqrt(beta+(y - YU)**2)))
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
        
class ThickFinWire(Wire):
    '''Wire of finite length and width, no height'''
    def __init__(self, name, length, width, height, current, x0,y0,z0):
        Wire.__init__(self,name,length,width,height,current,1)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
       
    def bfieldcalc(self, x, y, z):
        ''' Set up for a wire in the x, y plane parallel to x axis of length 2*L at y = R'''
        #aliases
        self.L = self.length / 2
        self.R = self.y0
        w = self.width
        I = self.current
        
        G = (mu_0 * I / (4*pi*w))
        B_y = -1*G * (__AuxBy__(x, y, z, w) + __AuxBy__(-1*x, y, z, w) - __AuxBy__(x, y, z, -1*w) - __AuxBy__(-1*x, y, z, -1*w))
        B_z = G * (__AuxBz__(x, y, z, w) + __AuxBz__(-1*x, y, z, w) - __AuxBz__(x, y, z, -1*w) - __AuxBz__(-1*x, y, z, -1*w))
        return np.array((0, B_y, B_z))
    
    def __AuxBy__(self, x, y, z, w):
        x0 = self.L + x
        y0 = self.R - y
        num = x0*(y0 + 0.5*w)
        den = z * sqrt(x0**2 + y0**2 + z**2)
        return atan(num / den)
        
    def __AuxBz__(self, x, y, z, w):
        x0 = self.L + x
        y0 = self.R - y
        num = sqrt(x0**2 + (y0 - 0.5*w)**2 + z**2)
        return atanh(num / x0)
        
class ThickFinHWire(ThickFinWire):
    def __init__(self, name, length, width, height, current, x0,y0,z0)
        super().__init__(self, name, length, width, height, current, x0,y0,z0)
    def bfieldcalc(self, x, y, z):
        return super().bfieldcalc(self, x - self.x0, y, z - self.z0)

class ThickFinVWire(ThickFinWire):
    '''stored in the rotated frame; never access properties directly, write appropriate getters'''
    def __init__(self, name, length, width, height, current, x0,y0,z0)
        super().__init__(self, name, length, width, height, current, y0,-1*x0,z0)
    def bfieldcalc(self, x, y, z):
        '''Need to rotate the x, y axes for this...'''
        rot_frame_field = super().bfieldcalc(self, x - self.x0, y, z - self.z0)
        return np.array((-1*rot_frame_field[1], 0, rot_frame_field[2]))