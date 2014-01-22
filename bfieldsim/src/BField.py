# -*- coding: utf-8 -*-
"""
Created on Fri Jun 07 22:31:28 2013

1D Bfield class

@author: Will
"""

import numpy as np
import matplotlib.pyplot as plt

class BField():
    def __init__(self, window):
        self.window = window.window
        self.Bx = np.zeros(self.window.shape)
        self.By = np.zeros(self.window.shape)
        self.Bz = np.zeros(self.window.shape)
    def set_field(self, x, vec):
        nearest_x_ind = np.argmin(np.absolute(self.window - x))
        self.Bx[nearest_x_ind] = vec[0]
        self.By[nearest_x_ind] = vec[1]
        self.Bz[nearest_x_ind] = vec[2]
    def plot_mag(self):
        B = np.sqrt(self.Bx**2 + self.By**2 + self.Bz**2)
        plt.plot(self.window, 1e4*B)
        plt.xlabel('X')
        plt.ylabel('Field (G)')
#        plt.show()
    def get_mag(self):
        B = np.sqrt(self.Bx**2 + self.By**2 + self.Bz**2)
        return B
    def add_bias(self, bias):
        for x in xrange(len(self.window)):
            self.Bx[x] += bias[0]
            self.By[x] += bias[1]
            self.Bz[x] += bias[2]
            
class BField_2D():
    def __init__(self, window):
        self.window = window.window
        print 'window shape'
        print self.window.shape
        self.X = window.X
        self.Y = window.Y
        self.Bx = np.zeros(self.window.shape)
        self.By = np.zeros(self.window.shape)
        self.Bz = np.zeros(self.window.shape)
    def set_field(self, pt, vec):
        A = np.sqrt((self.X - pt.x)**2 + (self.Y - pt.y)**2)
        nearest_pt_ind = np.unravel_index(A.argmin(), A.shape)
        self.Bx[nearest_pt_ind] = vec[0]
        self.By[nearest_pt_ind] = vec[1]
        self.Bz[nearest_pt_ind] = vec[2]
    def plot_mag(self):
        B = np.sqrt(self.Bx**2 + self.By**2 + self.Bz**2)
        B.reshape(self.window.shape)
        print B.shape
        plt.contourf(self.X, self.Y, 1e4*B)
        plt.colorbar()
        plt.xlabel('X')
        plt.ylabel('Field (G)')
#        plt.show()
    def get_mag(self):
        B = np.sqrt(self.Bx**2 + self.By**2 + self.Bz**2)
        return B
    def add_bias(self, bias):
        for ind in np.ndenumerate(self.window):
            self.Bx[ind] += bias[0]
            self.By[ind] += bias[1]
            self.Bz[ind] += bias[2]