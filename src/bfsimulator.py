# -*- coding: utf-8 -*-
"""
Created on 2013-07-10

Magnetic Trap Simulator, Python edition
This is a Python port of the most up-to-date simulator I had written in 
MatLab as of 7/4/10. It simulates
the combined fields from the macrowires and microwires.

This is a class-based implementation, unlike BFieldSim

@author: Will
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, atan2, ceil
from acmconstants import K_B, M
import logging
from scipy.optimize import minimize
import numdifftools as nd

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

clear = "\n"*100
um = 1e-6
mm = 1e-3
CONV_THRESH = 0.01 #percent change allowed in f_trans to declare convergence

class BFieldSimulator():
    def __init__(self, chip = None):
        self.resolution=.000025 # resolution of the plots
        self.plotleft=-1*mm # boundary of plots, made it smaller for small arm trap on 11/6/13
        self.plotright=1*mm
        self.nhorz = ceil((self.plotright - self.plotleft) / self.resolution)
        self.plottop = 1*mm
        self.plotbottom = -1*mm
        self.nvert = ceil((self.plottop - self.plotbottom) / self.resolution)
        self.x, self.y =np.meshgrid(np.arange(self.plotleft, self.plotright, self.resolution), 
                                    np.arange(self.plotbottom, self.plottop, self.resolution))
        self.plothigh = 5*mm
        self.plotlow = 1*um
        self.nz = ceil((self.plothigh - self.plotlow) / self.resolution)
        self.z_range = np.linspace(self.plothigh,self.plotlow,self.nz) #meters
        self.z_spacing = self.z_range[1]-self.z_range[0] #meters
        
        if chip is not None:
            self.wirespecs = chip['wirespecs'] #this is a bad kludge
            self.wirespecs = [wire for wire in self.wirespecs if wire.current != 0]
            self.B_ybias = chip['B_ybias']
            self.B_bias = np.array((chip['B_xbias'], chip['B_ybias'], chip['B_zbias']))
            self.x_trap = 0
            self.y_trap = 0
            self.calc_trap_height() #initialize B_tot_z
            self.find_z_min() #initialize z_trap
            self.calc_xy()# initialize B_tot_xy
            self.find_xy_min() #initialize x_trap, y_trap
    
        
    def set_chip(self, chip):
        '''Choose a chip configuration, defined in some other file via wirespecs'''
        self.__init__(chip)
        
    def zoom(self, zoom_factor, center_pt = [0, 0, 200e-6]):
        '''Zoom in on the trap center'''
        logging.debug('Zoom in on: %2.0f, %2.0f, %2.0f'%(self.x_trap*1e6, self.y_trap*1e6, self.z_trap*1e6))
        self.resolution = self.resolution / zoom_factor
        x_cen = self.x_trap
        y_cen = self.y_trap
        z_cen = self.z_trap
        self.plotleft = x_cen - 0.5*self.nhorz*self.resolution
        self.plotright = x_cen + 0.5*self.nhorz*self.resolution
        self.plottop = y_cen + 0.5*self.nvert*self.resolution
        self.plotbottom = y_cen - 0.5*self.nvert*self.resolution
        self.plothigh = z_cen + 0.5*self.nz*self.resolution
        self.plotlow = max(z_cen -0.5*self.nz*self.resolution, 1e-6)

        self.x, self.y =np.meshgrid(
                np.linspace(self.plotleft, self.plotright, self.nhorz), 
                np.linspace(self.plotbottom, self.plottop, self.nvert))
        self.z_range = np.linspace(self.plotlow,self.plothigh,self.nz) #meters
        self.z_spacing = self.resolution #meters
 
    def calc_trap_height(self):
        '''Find B field magnitude along z axis through trap center'''
        self.B_tot_z = np.zeros(len(self.z_range))
        for ii in xrange(len(self.z_range)):
            self.B_tot_z[ii] = self.calc_field_mag(self.x_trap, self.y_trap, self.z_range[ii])

    def find_z_min(self):
        '''Find minimum value of B field magnitude along z axis through trap center'''
        self.min_index = np.argmin(self.B_tot_z)           
        self.z_trap = self.z_range[self.min_index]

    def plot_z(self):
        '''Plot field against z through trap center'''
        self.calc_trap_height()
        plt.plot(self.z_range*1e3, self.B_tot_z*1e4)
        plt.xlabel('Z axis (mm)') #Standard axis labelling
        plt.ylabel('Effective |B|, (G)')
        plt.ylim((1e4*min(self.B_tot_z),1e4*min(self.B_tot_z[0], self.B_tot_z[-1])))
        plt.show()

    def calc_xz(self):
        '''Calculate B field magnitude in xz plane through trap center'''
        self.B_tot = np.zeros(self.x.shape)
        for coords in np.ndenumerate(self.x):    
            self.B_tot[coords[0]] = self.calc_field_mag(self.x[coords[0]], 
                                                  self.y_trap, 
                                                    self.z[coords[0]])
            
    def plot_xz(self):
        '''Plot field magnitude in xz plane through trap center'''
        self.calc_xz()
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_surface(self.x*1e3,self.z*1e3,self.B_tot_xy*1e4, rstride=8, cstride=8, alpha=0.3)
        plt.xlabel('X axis (mm)')
        plt.ylabel('Z axis (mm)') #standard axis labelling
        ax.set_zlabel('B field (G)')
        plt.show()
        
    def calc_xy(self):
        '''Calculate field magnitude in xy plane through trap center'''
        self.B_tot_xy = np.zeros(self.x.shape)
        for coords in np.ndenumerate(self.x):    
            self.B_tot_xy[coords[0]] = self.calc_field_mag(self.x[coords[0]], 
                                                  self.y[coords[0]], 
                                                    self.z_trap)
    def find_xy_min(self):
        '''Find minimum value of B field magnitude in xy plane through trap center'''
        min_ind = np.unravel_index(self.B_tot_xy.argmin(), self.B_tot_xy.shape)
        x_ind = min_ind[0]
        y_ind = min_ind[1]
        self.x_trap = self.x[0, y_ind]
        self.y_trap = self.y[x_ind, 0]
        logging.debug('x_trap: %f'%self.x_trap)
        logging.debug('y_trap: %f'%self.y_trap)
    
    def plot_xy(self):
        '''plot field magnitude in xy plane through trap center'''
        self.calc_xy()
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_surface(self.x*1e3,self.y*1e3,self.B_tot_xy*1e4, 
                        rstride=8, 
                        cstride=8, 
                        alpha=0.3)
        plt.xlabel('X axis (mm)')
        plt.ylabel('Y axis (mm)') #standard axis labelling
        ax.set_zlabel('B field (G)')
        plt.show()
    
    def calc_field(self, x, y, z):
        '''Calculate total B field due to all wires at a point'''
        tot_field = [0,0,0]
        for wire in self.wirespecs:
            this_field = wire.bfieldcalc(x, y, z)
            tot_field[0] += this_field[0]
            tot_field[1] += this_field[1]
            tot_field[2] += this_field[2]
        return np.array(tot_field)
    
    def calc_field_mag(self, x, y, z):
        '''Calculate total B field magnitude due to all wires at a point'''
        tot_field = self.calc_field(x, y, z)
        return np.linalg.norm(tot_field + self.B_bias)

    def calc_field_dir(self, x, y, z):
        '''Calculate direction in xy plane of B field due to all wires, at a point'''
        tot_field = self.calc_field(x, y, z)
        return atan2(tot_field[1], tot_field[0])
                    
    def find_trap_cen(self, method='3D', min_method = 'Nelder-Mead'):
        '''Find the location of the minimum B field magnitude'''
        if method == '1+2D':
            field_mag = lambda x: self.calc_field_mag(self.x_trap, self.y_trap, x)
            results = minimize(field_mag, self.z_trap, method = min_method, options = {'disp':True})
            [self.z_trap] = results.x
            field_mag = lambda x: self.calc_field_mag(x[0], x[1], self.z_trap)
            results = minimize(field_mag, (self.x_trap, self.y_trap), method = min_method, options = {'disp':True})
            [self.x_trap, self.y_trap] = results.x
        elif method == '1D':
            field_mag = lambda x: self.calc_field_mag(self.x_trap, self.y_trap, x)
            results = minimize(field_mag, self.z_trap, method = min_method, options = {'disp':True})
            [self.z_trap] = results.x
        else:
            field_mag = lambda x: self.calc_field_mag(x[0], x[1], x[2])
            results = minimize(field_mag, (self.x_trap, self.y_trap, self.z_trap), method = min_method, options = {'disp':True})
            [self.x_trap, self.y_trap, self.z_trap] = results.x

    def analyze_trap(self, method = '3D'):
        '''Extract trap frequencies'''
        trap_params = {}
       
        if method == '2D':
            field_mag_xy = lambda x: self.calc_field_mag(x[0], x[1], self.z_trap)
            trap_loc = [self.x_trap, self.y_trap]
        else:
            field_mag_xy = lambda x: self.calc_field_mag(x[0], x[1], x[2])
            trap_loc = [self.x_trap, self.y_trap, self.z_trap]
        Hxy = nd.Hessian(field_mag_xy)
        try:
            Principal_Matrix = Hxy(trap_loc)
            self.values, self.vectors = np.linalg.eig(Principal_Matrix)
            freqs = [(1.2754)*sqrt(abs(val)) for val in self.values]
        except IndexError:
            print('Frequency Finding Failure')
            freqs = [-1,-1,-1]
            
        
        self.f_transverse = max(freqs)
        self.f_longitudinal = min(freqs)
        self.omega_transverse = 2*pi*self.f_transverse
        self.omega_longitudinal = 2*pi*self.f_longitudinal
        
        if method == '2D':
            field_mag_z = lambda z: self.calc_field_mag(self.x_trap, self.y_trap, z)
            Hz = nd.Hessian(field_mag_z)
            ddBdzz = Hz([self.z_trap])[0][0]
            self.f_z = 1.2754*sqrt(abs(ddBdzz))
        else:
            self.f_z = sorted(freqs)[1]
        self.omega_z = 2*pi*self.f_z
        
        trap_params['h'] = self.z_trap
        trap_params['f_long'] = self.f_longitudinal
        trap_params['f_trans'] = self.f_transverse
        trap_params['f_z'] = self.f_z
        return trap_params
                       
    def plot_xy_dir(self):
        '''plot field direction in xy plane through trap center'''
        fig = plt.figure()
        plt.hsv

        plt.contourf(self.x*1e3,self.y*1e3,self.B_dir, 16, cmap=cm.get_cmap('binary'))
        plt.xlabel('X axis (mm)')
        plt.ylabel('Y axis (mm)') #standard axis labelling
        plt.colorbar()
        plt.show()
        
    def plot_xy_coupling(self):
        '''plot field coupling in xy plane through trap center'''
        fig = plt.figure()
        plt.hsv

        plt.contourf(self.x*1e3,self.y*1e3,np.cos(self.B_dir), 16, cmap=cm.get_cmap('binary'))
        plt.xlabel('X axis (mm)')
        plt.ylabel('Y axis (mm)') #standard axis labelling
        plt.colorbar()
        plt.show()
        
        
    def find_trap_freq(self, method = '3D', 
                            analyze_method = '3D', 
                            trap_find_method = '3D', 
                            min_method = 'Nelder-Mead', 
                            convcheck = False, 
                            debug = False):
        '''Iterate trap analysis until transverse frequency converges'''
        logging.debug('\nxtrap: %e\nytrap: %e'%(self.x_trap, self.y_trap))
        sim_results = self.analyze_trap()
        
        if method == '1D':
            error = abs(sim_results['f_z'] - sim_results['f_trans']) / sim_results['f_z']
            logging.debug('f_trans and f_z differ by: %2.1f %%'%(error*100))

        else:
            f_trans_prev = 0.1 # An unreasonably low frequency to start
            error = abs(sim_results['f_trans'] - f_trans_prev) / sim_results['f_trans']
            f_trans_prev = sim_results['f_trans']
            logging.debug('f_trans_prev and f_trans differ by: %2.1f %%'%(error*100))

            
        if debug:
            self.plot_z()
            self.plot_xy()
        
        logging.debug('x_trap : %2.0f um \n\ty_trap : %2.0f um \n\tz_trap : %2.0f um'%(self.x_trap*1e6, self.y_trap*1e6, self.z_trap*1e6))
        logging.debug('f_long : %2.0f Hz \n\tf_trans : %2.0f Hz \n\tf_z : %2.0f Hz'%(sim_results['f_long'], sim_results['f_trans'], sim_results['f_z']))
        n_tries = 1
        while error > CONV_THRESH and n_tries < 10:
            self.zoom(2) #used to be 4; 11/6/13
            self.find_trap_cen(method = trap_find_method, min_method = min_method)

            sim_results = self.analyze_trap()
            last_error = error
            if method == '1D':
                error = abs(sim_results['f_z'] - sim_results['f_trans']) / sim_results['f_z']
                logging.debug('f_trans and f_z differ by: %2.1f %%'%(error*100))

            else:
                error = abs(sim_results['f_trans'] - f_trans_prev) / sim_results['f_trans']
                f_trans_prev = sim_results['f_trans']
                logging.debug('f_trans_prev and f_trans differ by: %2.1f %%'%(error*100))

                
            if debug:
                self.plot_z()
                self.plot_xy()
            
            logging.debug('x_trap : %2.0f um \n\ty_trap : %2.0f um \n\tz_trap : %2.0f um'%(self.x_trap*1e6, self.y_trap*1e6, self.z_trap*1e6))
            logging.debug('f_long : %2.0f Hz \n\tf_trans : %2.0f Hz \n\tf_z : %2.0f Hz'%(sim_results['f_long'], sim_results['f_trans'], sim_results['f_z']))
            
            if convcheck:
                if abs(last_error - error) < .005:
                    print('Not Converging')
                    break
            n_tries += 1       
        return sim_results
        
if __name__ == '__main__':
    import atomchip_v3
    b_f_sim = BFieldSimulator()
    b_f_sim.set_chip(atomchip_v3.atomchip_v3)
    b_f_sim.calc_trap_height()
    b_f_sim.plot_z()
    sim_results = b_f_sim.find_trap_freq()
    print 'x_trap : %2.0f um \ny_trap : %2.0f um \nz_trap : %2.0f um'%(b_f_sim.x_trap*1e6, b_f_sim.y_trap*1e6, b_f_sim.z_trap*1e6)
    print '\nf_long : %2.0f Hz \nf_trans : %2.0f Hz \nf_z : %2.0f Hz'%(sim_results['f_long'], sim_results['f_trans'], sim_results['f_z'])
    b_f_sim.set_chip(atomchip_v3.atomchip_v3)
    b_f_sim.plot_xy()