# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:54:24 2013

@author: Will
"""

from math import pi, sqrt, atan, atanh, log
import numpy as np
MU_0 = 4 * pi * 1e-7  # Permeability of Free Space


class Wire():

    '''Basic class describing a wire of rectangular cross section
    with uniform current density.'''
    # subwires is still here only for backwards compatibility

    def __init__(self, name, length, width, height, current, subwires=1):
        self.name = name
        self.length = length
        self.width = width
        self.height = height
        self.current = current
        # This pre-factor shows up in all the field formulae
        self.const_G = MU_0 * current / (4 * pi)

    def bfieldcalc(self, x, y, z):
        '''Based on a finite thin wire parallel to x axis
        and centered on the y-axis in the x-y plane'''
        xL = x
        xR = x - self.length

        beta = z ** 2 + y ** 2
        B_G = self.const_G * (xL / (beta * sqrt(beta +
                                                xL ** 2)) - xR
                                                / (beta * sqrt(beta + xR ** 2)))
        return (0, -1 * B_G * z, B_G * y)


class HThinWire(Wire):

    '''Define a wire with current parallel to the x axis.'''

    def __init__(self, name, length, width, height,
                 current, xl, y0, z0, subwires=1):
        Wire.__init__(self, name, length, width, height, current, subwires)
        self.xl = xl
        self.y0 = y0
        self.z0 = z0

    def bfieldcalc(self, x, y, z):
        return Wire.bfieldcalc(self, x - self.xl, y - self.y0, z - self.z0)


class VThinWire(Wire):

    '''Define a wire with current parallel to the y axis.'''

    def __init__(self, name, length, width, height,
                 current, x0, y0, z0, subwires=1):
        Wire.__init__(self, name, length, width, height, current, subwires)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0

    def bfieldcalc(self, x, y, z):
        '''Rotate the coordinates to use the basic field calculator,
        then rotate the field vector back.'''
        rot_frame_field = Wire.bfieldcalc(
            self,
            y -
            self.y0,
            self.x0 -
            x,
            z -
            self.z0)
        return (-1 * rot_frame_field[1], 0, rot_frame_field[2])


class ThickFinWire(Wire):

    '''Wire of finite length and width, no height'''

    def __init__(self, name, length, width, height, current, x0, y0, z0):
        Wire.__init__(self, name, length, width, height, current)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.L = self.length / 2
        self.R = self.y0 + self.width / 2

    def bfieldcalc(self, x, y, z):
        ''' Set up for a wire in the x, y plane parallel to x axis of
        length 2*L at y = R. Formula is taken from Extravour thesis appendix.'''
        w = self.width
        I = self.current

        G = (MU_0 * I / (4 * pi * w))
        B_y = -1 * G * (self.__AuxBy__(x, y, z, w) +
                        self.__AuxBy__(-1 * x, y, z, w)
                        - self.__AuxBy__(x, y, z, -1 * w) -
                        self.__AuxBy__(-1 * x, y, z, -1 * w))
        # B_z = G * (self.__AuxBz__(x, y, z, w) + self.__AuxBz__(-1*x, y, z, w)
        #- self.__AuxBz__(x, y, z, -1*w) - self.__AuxBz__(-1*x, y, z, -1*w))
        # reorganized to avoid imaginary sub-results
        B_z = G * \
            (self.__AuxBz2__(x, y, z, w) + self.__AuxBz2__(-1 * x, y, z, w))
        return (0, B_y, B_z)

    def __AuxBy__(self, x, y, z, w):
        x0 = x + self.L
        y0 = self.R - y
        num = x0 * (y0 + 0.5 * w)
        den = z * sqrt(x0 ** 2 + y0 ** 2 + z ** 2)
        if den == 0:
            # m^2, correct length scale for 'small' distances (sub-um)
            den = 1e-14
        return atan(num / den)

    def __AuxBz__(self, x, y, z, w):
        x0 = x + self.L
        y0 = self.R - y
        num = sqrt(x0 ** 2 + (y0 + 0.5 * w) ** 2 + z ** 2)
        # scipy.optimize.minimize won't catch the ZeroDivisionError so I added
        # this.
        if x0 == 0:
            x0 = 1e-7  # m, a reasonable 'small' distance
        try:
            return (num / x0)
        except ZeroDivisionError:
            # print('Zero Division Error')
            return num / 1e-7

    def __AuxBz2__(self, x, y, z, w):
        u = self.__AuxBz__(x, y, z, w)
        v = self.__AuxBz__(x, y, z, -1 * w)
        arg = (u - v) / (1 - u * v)
        return atanh(arg)


class HThickFinWire(ThickFinWire):

    def __init__(self, name, length, width, height, current, x0, y0, z0):
        ThickFinWire.__init__(
            self,
            name,
            length,
            width,
            height,
            current,
            x0,
            y0,
            z0)

    def bfieldcalc(self, x, y, z):
        '''Translate the coordinates appropriately, then calculate the field'''
        return (
            ThickFinWire.bfieldcalc(
                self,
                x - (self.x0 + self.L),
                y,
                z - self.z0)
        )


class VThickFinWire(ThickFinWire):

    def __init__(self, name, length, width, height, current, x0, y0, z0):
        ThickFinWire.__init__(
            self,
            name,
            length,
            width,
            height,
            current,
            x0,
            y0,
            z0)
        self.R = -1 * (self.x0 + self.width / 2)

    def bfieldcalc(self, x, y, z):
        '''Translate and rotate the coordinates to use the basic field
        calculator, then rotate the field vector back.'''
        rot_frame_field = ThickFinWire.bfieldcalc(
            self,
            y - (self.y0 + self.L),
            - x,
            z - self.z0)
        return (-1 * rot_frame_field[1], 0, rot_frame_field[2])


class RectWire(Wire):

    '''Wire of finite length, width, and height'''

    def __init__(self, name, length, width, height, current, x0, y0, z0):
        Wire.__init__(self, name, length, width, height, current)

        self.xlims = [x0, x0 + self.length]
        self.ylims = [y0, y0 + self.width]
        self.zlims = [z0, z0 + self.height]

        self.const = (MU_0 * self.current /
                      (4 * pi * self.width * self.height))

    def bfieldcalc(self, x, y, z):
        '''Set up for a wire with current in the x direction; formula from
        Treutlein thesis appendix'''
        b = [0, 1]
        indices_set = [(i, j, k) for i in b for j in b for k in b]
        y_terms = [(
            (-1) ** sum(indices)) * self.__AuxFn__(x - self.xlims[indices[0]],
                                                   y -
                                                   self.ylims[
                                                       indices[1]],
                                                   z - self.zlims[indices[2]])
            for indices in indices_set]
        z_terms = [(
            (-1) ** sum(indices)) * self.__AuxFn__(x - self.xlims[indices[0]],
                                                   z -
                                                   self.zlims[
                                                       indices[2]],
                                                   y - self.ylims[indices[1]])
            for indices in indices_set]
        By = -1 * self.const * sum(y_terms)
        Bz = self.const * sum(z_terms)
        return (0, By, Bz)

    def __AuxFn__(self, x, y, z):
        r = sqrt(x ** 2 + y ** 2 + z ** 2)
        if z * r == 0:
            t1 = z * atan(x * y / 1e-14)
        else:
            t1 = z * atan(x * y / (z * r))
        t2 = -1 * x * log(y + r)
        t3 = -1 * y * log(x + r)
        return (t1 + t2 + t3)


class HRectWire(RectWire):

    def __init__(self, name, length, width, height,
                 current, x0, y0, z0, subwires=1):
        RectWire.__init__(
            self,
            name,
            length,
            width,
            height,
            current,
            x0,
            y0,
            z0)

    def bfieldcalc(self, x, y, z):
        return RectWire.bfieldcalc(self, x, y, z)


class VRectWire(RectWire):

    def __init__(self, name, length, width, height,
                 current, x0, y0, z0, subwires=1):
        '''Limits are stored as in rotated frame'''
        RectWire.__init__(
            self,
            name,
            length,
            width,
            height,
            current,
            x0,
            y0,
            z0)
        self.xlims = [y0, y0 + self.length]
        self.ylims = [-1 * x0, -1 * (x0 + self.width)]

    def bfieldcalc(self, x, y, z):
        '''Rotate the coordinates to use the basic field calculator,
        then rotate the field vector back.'''
        rot_frame_field = RectWire.bfieldcalc(self, y, -1 * x, z)
        return (-1 * rot_frame_field[1], 0, rot_frame_field[2])
