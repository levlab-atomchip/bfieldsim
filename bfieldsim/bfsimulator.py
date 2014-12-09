# -*- coding: utf-8 -*-
"""
Created on 2013-07-10

Magnetic Trap Simulator, Python edition
This is a Python port of the most up-to-date simulator I had written in
MatLab as of 7/4/10. It simulates
the combined fields from the macrowires and microwires.

This is a class-based implementation, unlike the old BFieldSim

@author: Will
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, atan2, ceil
from .acmconstants import K_B, M, G, MU_B, M_F, G_F
# Boltzmann constant, mass of Rb-87, Gravitational acceleration,
# Bohr Magneton, m_f quantum number, g-factor
import logging
from scipy.optimize import minimize
import numdifftools as nd

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

clear = "\n" * 100
um = 1e-6
mm = 1e-3
CONV_THRESH = 0.01  # percent change allowed in f_trans to declare convergence


class BFieldSimulator():

    def __init__(self, chip=None):
        '''Set up variables to define window and initialize
        the fields and trap location.'''
        self.resolution = .0001  # resolution of the plots; was 25 um
        self.plotleft = -5 * mm
        # boundary of plots, made it smaller for small arm trap on 11/6/13
        self.plotright = 5 * mm
        self.nhorz = ceil((self.plotright - self.plotleft) / self.resolution)
        self.plottop = 5 * mm
        self.plotbottom = -5 * mm
        self.nvert = ceil((self.plottop - self.plotbottom) / self.resolution)
        self.x, self.y = np.meshgrid(np.arange(self.plotleft,
                                               self.plotright, self.resolution),
                                     np.arange(self.plotbottom,
                                               self.plottop, self.resolution))
        self.plothigh = 5 * mm
        self.plotlow = 1 * um
        self.nz = ceil((self.plothigh - self.plotlow) / self.resolution)
        self.z_range = np.linspace(
            self.plothigh,
            self.plotlow,
            self.nz)  # meters
        self.z_spacing = self.z_range[1] - self.z_range[0]  # meters

        if chip is not None:
            self.wirespecs = chip['wirespecs']  # this is a bad kludge
            self.wirespecs = [wire for wire in self.wirespecs
                              if wire.current != 0]
            self.B_ybias = chip['B_ybias']
            self.B_bias = np.array((chip['B_xbias'],
                                    chip['B_ybias'],
                                    chip['B_zbias']))
            self.x_trap = 0
            self.y_trap = 0
            self.calc_trap_height()  # initialize B_tot_z
            self.find_z_min()  # initialize z_trap
            self.calc_xy()  # initialize B_tot_xy
            self.find_xy_min()  # initialize x_trap, y_trap

    def set_chip(self, chip):
        '''Choose a chip configuration, defined
        in some other file via wirespecs.'''
        self.__init__(chip)

    def zoom(self, zoom_factor, center_pt=[0, 0, 200e-6]):
        '''Zoom in on the trap center.

        Args:
            zoom_factor (float): Resolution will increase by this factor.
            center_pt ([float, float, float]): The window center
                                            will move to this point.

        Effects:
            Update all window attributes.

        '''
        logging.debug('Zoom in on: %2.0f, %2.0f, %2.0f' % (self.x_trap * 1e6,
                                                           self.y_trap * 1e6,
                                                           self.z_trap * 1e6))
        self.resolution = self.resolution / zoom_factor
        x_cen = self.x_trap
        y_cen = self.y_trap
        z_cen = self.z_trap
        self.plotleft = x_cen - 0.5 * self.nhorz * self.resolution
        self.plotright = x_cen + 0.5 * self.nhorz * self.resolution
        self.plottop = y_cen + 0.5 * self.nvert * self.resolution
        self.plotbottom = y_cen - 0.5 * self.nvert * self.resolution
        self.plothigh = z_cen + 0.5 * self.nz * self.resolution
        self.plotlow = max(z_cen - 0.5 * self.nz * self.resolution, 1e-6)

        self.x, self.y = np.meshgrid(
            np.linspace(self.plotleft, self.plotright, self.nhorz),
            np.linspace(self.plotbottom, self.plottop, self.nvert))
        self.z_range = np.linspace(
            self.plotlow,
            self.plothigh,
            self.nz)  # meters
        self.z_spacing = self.resolution  # meters

    def calc_trap_height(self):
        '''Find B field magnitude along z axis through trap center.

        Effects:
            Update B_tot_z attribute.

        '''
        self.B_tot_z = np.zeros(len(self.z_range))
        for ii in xrange(len(self.z_range)):
            self.B_tot_z[ii] = self.calc_eff_bmag(self.x_trap,
                                                  self.y_trap,
                                                  self.z_range[ii])

    def find_z_min(self):
        '''Find minimum value of B field magnitude
        along z axis through trap center.

        Effects:
            Update z_trap attribute.
        '''
        self.min_index = np.argmin(self.B_tot_z)
        self.z_trap = self.z_range[self.min_index]

    def plot_z(self):
        '''Plot field against z through trap center.'''
        self.calc_trap_height()
        plt.plot(self.z_range * 1e3, self.B_tot_z * 1e4)
        # Scales length to mm and field to Gauss
        plt.xlabel('Z axis (mm)')  # Standard axis labelling
        plt.ylabel('Effective |B|, (G)')
        plt.show()

    def calc_xz(self):
        '''Calculate B field magnitude in xz plane through trap center

        Effects:
            Update B_tot_xz attribute.
        '''
        self.B_tot_xz = np.zeros(self.x.shape)
        for coords in np.ndenumerate(self.x):
            self.B_tot_xz[coords[0]] = self.calc_eff_bmag(self.x[coords[0]],
                                                          self.y_trap,
                                                          self.z[coords[0]])

    def plot_xz(self):
        '''Plot field magnitude in xz plane through trap center'''
        self.calc_xz()
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_surface(
            self.x * 1e3, self.z * 1e3, self.B_tot_xy * 1e4, rstride=8,
            cstride=8,
            alpha=0.3)
        plt.xlabel('X axis (mm)')
        plt.ylabel('Z axis (mm)')  # standard axis labelling
        ax.set_zlabel('B field (G)')
        plt.show()

    def calc_xy(self):
        '''Calculate field magnitude in xy plane through trap center.

        Effects:
            Updates B_tot_xy attribute.
        '''
        self.B_tot_xy = np.zeros(self.x.shape)
        for coords in np.ndenumerate(self.x):
            self.B_tot_xy[coords[0]] = self.calc_eff_bmag(self.x[coords[0]],
                                                          self.y[coords[0]],
                                                          self.z_trap)

    def find_xy_min(self):
        '''Find minimum value of B field magnitude in
        xy plane through trap center.

        This method relies on B_tot_xy being up to date.
        It just finds the minimum by brute force.

        Effects:
            Updates x_trap, y_trap attributes.
        '''
        min_ind = np.unravel_index(self.B_tot_xy.argmin(), self.B_tot_xy.shape)
        x_ind = min_ind[0]
        y_ind = min_ind[1]
        self.x_trap = self.x[0, y_ind]
        self.y_trap = self.y[x_ind, 0]
        logging.debug('x_trap: %f' % self.x_trap)
        logging.debug('y_trap: %f' % self.y_trap)

    def plot_xy(self):
        '''plot field magnitude in xy plane through trap center'''
        self.calc_xy()
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_surface(self.x * 1e3, self.y * 1e3, self.B_tot_xy * 1e4,
                        rstride=8,
                        cstride=8,
                        alpha=0.3)
        plt.xlabel('X axis (mm)')
        plt.ylabel('Y axis (mm)')  # standard axis labelling
        ax.set_zlabel('B field (G)')
        plt.show()

    def calc_field(self, x, y, z):
        '''Calculate total B field due to all wires
        and bias fields at a point.'''
        tot_field = [0, 0, 0]
        for wire in self.wirespecs:
            this_field = wire.bfieldcalc(x, y, z)
            tot_field[0] += this_field[0]
            tot_field[1] += this_field[1]
            tot_field[2] += this_field[2]
        return (np.array(tot_field) + self.B_bias)

    def calc_field_mag(self, x, y, z):
        '''Calculate total B field magnitude due to
        all wires and bias fields at a point.'''
        tot_field = self.calc_field(x, y, z)
        return np.linalg.norm(tot_field)

    def calc_eff_bmag(self, x, y, z):
        '''Calculate total potential, expressed as an effective field,
        due to magnetic field and gravity (referenced to z = 0)'''
        field_mag = self.calc_field_mag(x, y, z)
        field_pot = M_F * G_F * MU_B * field_mag
        # Because g_f = 1/2, m_f = 2 for the 2, 2 state
        gravity_pot = -1 * M * G * z  # -1 because z axis is 'upside-down'
        total_pot = field_pot + gravity_pot
        return total_pot / (M_F * G_F * MU_B)

    def calc_field_dir(self, x, y, z):
        '''Calculate direction in xy plane of B field
        due to all wires, at a point.'''
        tot_field = self.calc_field(x, y, z)
        return atan2(tot_field[1], tot_field[0])

    def find_trap_cen(self, method='3D', min_method='Nelder-Mead'):
        '''Find the location of the minimum B field magnitude.

        This method uses a numerical optimization algorithm to find the trap.
        It is generally faster and more
        accurate than the old brute force method.

        Args:
            method (str): Sets the trap finding method from
            ('1D', '1+2D', '3D'). Defaults to '3D'.
                '1D' - Constrain the search to the z-axis running through
                        the pre-existing x_trap, y_trap.
                        This was mostly used during testing,
                        no practical application.
                '1+2D' - First find the optimum height,
                        using pre-existing x_trap and y_trap.
                        Then find the minimum within the xy plane
                        at that height.
                        Similar to the old trap-finding method.
                        Potentially faster but less accurate.
                '3D' - Find the trap, no constraints.

            min_method (str): Sets the optimization method used by
                        scipy.optimize.minimize. Defaults to 'Nelder-Mead'.
                'Nelder-Mead' - A fast, simplex-type algorithm.
                                Doesn't need to calculate any derivatives.
                                Not so accurate.
                'BFGS' - A slow algorithm that numerically calculates
                        the jacobian and hessian. The default for
                        scipy.optimize.minimize. Very accurate,
                        but the first run can take a while.

                See the scipy.optimize.minimize documentation for
                the other available methods (which I haven't tried)

        Effects:
            Updates x_trap, y_trap, z_trap attributes.
        '''
        if method == '1+2D':
            eff_bmag = lambda x: self.calc_eff_bmag(self.x_trap,
                                                    self.y_trap,
                                                    x)
            results = minimize(eff_bmag, self.z_trap, method=min_method,
                               options={'disp': True})
            [self.z_trap] = results.x
            eff_bmag = lambda x: self.calc_eff_bmag(x[0], x[1], self.z_trap)
            results = minimize(eff_bmag, (self.x_trap, self.y_trap),
                               method=min_method, options={'disp': True})
            [self.x_trap, self.y_trap] = results.x
        elif method == '1D':
            eff_bmag = lambda x: self.calc_eff_bmag(
                self.x_trap,
                self.y_trap,
                x)
            results = minimize(eff_bmag, self.z_trap, method=min_method,
                               options={'disp': True})
            [self.z_trap] = results.x
        else:
            eff_bmag = lambda x: self.calc_eff_bmag(x[0], x[1], x[2])
            results = minimize(eff_bmag, (self.x_trap, self.y_trap,
                                          self.z_trap), method=min_method,
                                                        options={'disp': True})
            [self.x_trap, self.y_trap, self.z_trap] = results.x

    def analyze_trap(self, method='3D'):
        '''Extract trap frequencies.

        This method uses a numerical differentiation module to calculate the
        Hessian matrix, and then extracts
        the principal axes and frequencies from that. It is faster and more
        reliable than the naive numerical
        differentiation I used in old versions of this code.

        Args:
            method (str): Sets the method for calculating trap frequencies.
                        Defaults to '3D'.
                '2D' - Calculate the frequencies in the x-y plane
                        and z-axis separately.
                        Not guaranteed to produce the true trap frequencies,
                        especially if there is
                        a vertical tilt. But it works fine for most traps.
                '3D' - Calculate the true principal frequencies.

        Returns:
            trap_params: A dictionary of calculated trap parameters.
                        This format made it easy to test the module.
                h - The trap height in microns.
                f_long - The smallest of the calculated frequencies, in Hz.
                f_trans - The largest of the calculated frequencies, in Hz.
                f_z - The middle of the calculated frequencies, in Hz.

        '''
        trap_params = {}

        freq_prefactor = sqrt((G_F * M_F * MU_B) / M) / (2 * pi)
        # convert eff_bmag to frequency

        if method == '2D':
            eff_bmag_xy = lambda x: self.calc_eff_bmag(x[0], x[1], self.z_trap)
            trap_loc = [self.x_trap, self.y_trap]
        else:
            eff_bmag_xy = lambda x: self.calc_eff_bmag(x[0], x[1], x[2])
            trap_loc = [self.x_trap, self.y_trap, self.z_trap]
        Hxy = nd.Hessian(eff_bmag_xy)
        try:
            Principal_Matrix = Hxy(trap_loc)
            self.values, self.vectors = np.linalg.eig(Principal_Matrix)
            freqs = [freq_prefactor * sqrt(abs(val)) for val in self.values]
        except IndexError:
            print('Frequency Finding Failure')
            freqs = [-1, -1, -1]

        self.f_transverse = max(freqs)
        self.f_longitudinal = min(freqs)
        self.omega_transverse = 2 * pi * self.f_transverse
        self.omega_longitudinal = 2 * pi * self.f_longitudinal

        if method == '2D':
            eff_bmag_z = lambda z: self.calc_eff_bmag(self.x_trap,
                                                      self.y_trap,
                                                      z)
            Hz = nd.Hessian(eff_bmag_z)
            ddBdzz = Hz([self.z_trap])[0][0]
            self.f_z = freq_prefactor * sqrt(abs(ddBdzz))
        else:
            self.f_z = sorted(freqs)[1]
        self.omega_z = 2 * pi * self.f_z

        trap_params['h'] = self.z_trap
        trap_params['f_long'] = self.f_longitudinal
        trap_params['f_trans'] = self.f_transverse
        trap_params['f_z'] = self.f_z
        return trap_params

    def plot_xy_dir(self):
        '''plot field direction in xy plane through trap center'''
        fig = plt.figure()
        plt.hsv

        plt.contourf(self.x * 1e3, self.y * 1e3, self.B_dir, 16,
                     cmap=cm.get_cmap('binary'))
        plt.xlabel('X axis (mm)')
        plt.ylabel('Y axis (mm)')  # standard axis labelling
        plt.colorbar()
        plt.show()

    def plot_xy_coupling(self):
        '''plot field coupling in xy plane through trap center'''
        fig = plt.figure()
        plt.hsv

        plt.contourf(self.x * 1e3, self.y * 1e3, np.cos(self.B_dir), 16,
                     cmap=cm.get_cmap('binary'))
        plt.xlabel('X axis (mm)')
        plt.ylabel('Y axis (mm)')  # standard axis labelling
        plt.colorbar()
        plt.show()

    def find_trap_freq(self, method='3D',
                       analyze_method='3D',
                       trap_find_method='3D',
                       min_method='Nelder-Mead',
                       convcheck=False,
                       debug=False):
        '''Iterate trap analysis until transverse frequency converges.

        Args:
            method (str): Sets the criterion for declaring convergence.
                        Defaults to '3D'.
                '1D' - Declare convergence when f_trans approaches f_z.
                        Appropriate for traps with cylindrical symmetry.
                '3D' - Declare convergence when f_trans approaches
                        a limiting value.
            analyze_method (str): Sets the method for analyze_trap.
                                    Defaults to '3D'.
            trap_find_method (str): Sets the method for find_trap_cen.
                                    Defaults to '3D'.
            min_method (str): Sets the method for scipy.optimize.minimize.
                                    Defaults to 'Nelder-Mead'.
            convcheck (bool): Controls whether the loop terminates
                            if error stops improving. Defaults to False.
            debug (bool): Controls whether plots are displayed
                            while the loop is running. Defaults to False.

        Returns:
            sim_results: Output from analyze_trap.

        Effects:
            All the effects of the subroutines.
        '''
        logging.debug('\nxtrap: %e\nytrap: %e' % (self.x_trap, self.y_trap))
        sim_results = self.analyze_trap()

        if method == '1D':
            error = abs(sim_results['f_z'] - sim_results['f_trans']) \
                / sim_results['f_z']
            logging.debug(
                'f_trans and f_z differ by: %2.1f %%' %
                (error * 100))

        else:
            f_trans_prev = 0.1  # An unreasonably low frequency to start
            error = abs(sim_results['f_trans'] - f_trans_prev) \
                / sim_results['f_trans']
            f_trans_prev = sim_results['f_trans']
            logging.debug('f_trans_prev and f_trans differ by: %2.1f %%'
                          % (error * 100))

        if debug:
            self.plot_z()
            self.plot_xy()

        logging.debug('x_trap : %2.0f um \n\ty_trap : \
                        %2.0f um \n\tz_trap : %2.0f um'
                      % (self.x_trap * 1e6, self.y_trap * 1e6,
                        self.z_trap * 1e6))
        logging.debug('f_long : %2.0f Hz \n\tf_trans : \
                        %2.0f Hz \n\tf_z : %2.0f Hz'
                      % (sim_results['f_long'], sim_results['f_trans'],
                         sim_results['f_z']))
        n_tries = 1
        while error > CONV_THRESH and n_tries < 10:
            self.zoom(2)  # used to be 4; 11/6/13
            self.find_trap_cen(method=trap_find_method,
                               min_method=min_method)

            sim_results = self.analyze_trap(method=analyze_method)
            last_error = error
            if method == '1D':
                error = abs(sim_results['f_z'] - sim_results['f_trans']) \
                    / sim_results['f_z']
                logging.debug(
                    'f_trans and f_z differ by: %2.1f %%' %
                    (error * 100))

            else:
                error = abs(sim_results['f_trans'] - f_trans_prev)\
                    / sim_results['f_trans']
                f_trans_prev = sim_results['f_trans']
                logging.debug('f_trans_prev and f_trans differ by: %2.1f %%'
                              % (error * 100))

            if debug:
                self.plot_z()
                self.plot_xy()

            logging.debug('x_trap : %2.0f um \n\ty_trap : \
                            %2.0f um \n\tz_trap : %2.0f um' %
                          (self.x_trap * 1e6, self.y_trap * 1e6,
                            self.z_trap * 1e6))
            logging.debug('f_long : %2.0f Hz \n\tf_trans : \
                            %2.0f Hz \n\tf_z : %2.0f Hz' %
                          (sim_results['f_long'],
                           sim_results['f_trans'],
                           sim_results['f_z']))

            if convcheck:
                if abs(last_error - error) < .005:
                    print('Not Converging')
                    break
            n_tries += 1
        return sim_results

if __name__ == '__main__':
    from atomchip_rectwires import atomchip
    b_f_sim = BFieldSimulator()
    b_f_sim.set_chip(atomchip)
    sim_results = b_f_sim.find_trap_freq()
    print 'x_trap : %2.0f um \ny_trap : %2.0f um \nz_trap : \
            %2.0f um' % (b_f_sim.x_trap * 1e6,
                         b_f_sim.y_trap * 1e6,
                         b_f_sim.z_trap * 1e6)
    print '\nf_long : %2.0f Hz \nf_trans : %2.0f Hz \nf_z : %2.0f Hz'\
        % (sim_results['f_long'], sim_results['f_trans'], sim_results['f_z'])
    b_f_sim.plot_xy()
