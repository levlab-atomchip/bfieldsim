# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:48:36 2013

@author: Will
"""

import unittest
import bfsimulator
from acmconstants import MU_0
from math import pi, sqrt
from acwires import HWire, VWire, NWire, HThinWire, VThinWire
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

class TestAgainstLit(unittest.TestCase):

    def setUp(self):
        schneider = {'name' : 'Schneider',
                     'f_long' : 18, 
                     'f_trans' : 51,
                     'f_z' : 51,
                     'h' : 3000e-6, 
                     'l' : 6.7e-3, 
                     'I' : 49.7, 
                     'B_zbias' : 0, 
                     'B_xbias' : 0,
                     'B_ybias' : 25.4e-4}
        hansel_I = {'name' : 'Hansel I',
                    'f_long' : 28, 
                     'f_trans' : 220, 
                     'f_z' : 220,
                     'h' : 445e-6, 
                     'l' : 1.95e-3, 
                     'I' : 2, 
                     'B_zbias' : 0, 
                     'B_xbias' : 0,
                     'B_ybias' : 8e-4}
        hansel_II = {'name' : 'Hansel II',
                    'f_long' : 17, 
                     'f_trans' : 6200, 
                     'f_z' : 6200,
                     'h' : 70e-6, 
                     'l' : 1.95e-3, 
                     'I' : 2, 
                     'B_zbias' : 0, 
                     'B_xbias' : 1.9e-4,
                     'B_ybias' : 55e-4}
        wildermuth = {'name' : 'Wildermuth',
                    'f_long' : 39, 
                     'f_trans' : 132, 
                     'f_z' : 132,
                     'h' : 1525e-6, 
                     'l' : 4e-3, 
                     'I' : 60, 
                     'B_zbias' : 0, 
                     'B_xbias' : 0,
                     'B_ybias' : 60e-4}
        treutlein_z1 = {'name' : 'Treutlein Z1',
                    'f_long' : 28, 
                     'f_trans' : 400, 
                     'f_z' : 400,
                     'h' : 350e-6, 
                     'l' : 2e-3, 
                     'I' : 2, 
                     'B_zbias' : 0, 
                     'B_xbias' : 0,
                     'B_ybias' : 10.8e-4}
        treutlein_z2 = {'name' : 'Treutlein Z2',
                    'f_long' : 17, 
                     'f_trans' : 8600, 
                     'f_z' : 8600,
                     'h' : 70e-6, 
                     'l' : 2e-3, 
                     'I' : 2, 
                     'B_zbias' : 0, 
                     'B_xbias' : 0.8e-4,
                     'B_ybias' : 55e-4}
        treutlein_z3 = {'name' : 'Treutlein Z3',
                    'f_long' : 13, 
                     'f_trans' : 3200, 
                     'f_z' : 3200,
                     'h' : 100e-6, 
                     'l' : 2e-3, 
                     'I' : 1, 
                     'B_zbias' : 0, 
                     'B_xbias' : 0.4e-4,
                     'B_ybias' : 20e-4}
        self.lit_chips = [hansel_I, hansel_II, wildermuth, treutlein_z1, treutlein_z2, treutlein_z3]
        self.lit_chips = map(self.make_chip_wires, self.lit_chips)
        self.bfsim = bfsimulator.BFieldSimulator()
        
    def make_chip_wires(self, chip):
        hwires = []
        vwires = []
        #(name, length, width, height, current, xl, y0, z0, subwires = 1):
        hwires.append(HThinWire('Central Test Wire', chip['l'], 100e-6, 5e-6, chip['I'], 
                             -0.5*chip['l'],0, 0))
        
        vwires.append(VThinWire('Arm 1', 10*chip['l'], 100e-6, 5e-6, chip['I'], -0.5*chip['l'], -10*chip['l'], 0))
        vwires.append(VThinWire('Arm 2', 10*chip['l'], 100e-6, 5e-6, chip['I'], 0.5*chip['l'], 0, 0))
        
        chip['wirespecs'] = hwires + vwires
        return chip

    def test_against_lit(self):
        for chip in self.lit_chips:
            print '%s Chip:'%chip['name']
            self.bfsim.set_chip(chip)
#            self.bfsim.calc_trap_height()
#            self.bfsim.calc_xy()
##            logging.debug('\nxtrap: %e\nytrap: %e'%(self.bfsim.x_trap, self.bfsim.y_trap))
##            self.bfsim.plot_z()
##            self.bfsim.plot_xy()
#            sim_results = self.bfsim.analyze_trap()
#            
#            f_z_trans_diff = abs(sim_results['f_z'] - sim_results['f_trans']) / sim_results['f_z']
#            n_tries = 1
#            while f_z_trans_diff > 0.05 and n_tries < 10:
#                logging.debug('f_z and f_trans differ by: %2.1f %%'%(f_z_trans_diff*100))
#                self.bfsim.zoom(4)
#                self.bfsim.calc_xy()
#                sim_results = self.bfsim.analyze_trap()
#                f_z_trans_diff = abs(sim_results['f_z'] - sim_results['f_trans']) / sim_results['f_z']
#                n_tries += 1
#            logging.debug('\nxtrap: %e\nytrap: %e'%(self.bfsim.x_trap, self.bfsim.y_trap))
            sim_results = self.bfsim.find_trap_freq()
#            self.bfsim.plot_xy()
#            sim_results = self.bfsim.analyze_trap()
#            print self.bfsim.B_tot
#            print self.bfsim.B_tot.shape
#            logging.debug('\nxtrap: %e\nytrap: %e'%(self.bfsim.x_trap, self.bfsim.y_trap))
            f_z_error = abs(sim_results['f_z'] - chip['f_z']) / chip['f_z']
            f_trans_error = abs(sim_results['f_trans'] - chip['f_trans']) / chip['f_trans']
            f_long_error = abs(sim_results['f_long'] - chip['f_long']) / chip['f_long']
            height_error = abs(sim_results['h'] - chip['h']) / chip['h']
            self.assertLessEqual(f_z_error, 0.10, 
                'Vertical Frequency:\nTrue: %2.0f Hz\nSim:  %2.0f Hz'%(chip['f_z'], sim_results['f_z']))
            self.assertLessEqual(f_trans_error, 0.10, 
                'Transverse Frequency:\nTrue: %2.0f Hz\nSim:  %2.0f Hz'%(chip['f_trans'], sim_results['f_trans']))
#            self.assertLessEqual(f_long_error, 0.10,
#                'Longitudinal Frequency:\nTrue: %2.0f Hz\nSim:  %2.0f Hz'%(chip['f_long'], sim_results['f_long']))
            self.assertLessEqual(height_error, 0.10,
                'Trap Height:\nTrue: %e\nSim:  %e'%(chip['h'], sim_results['h']))
            print 'OK'

            
        
if __name__ == '__main__':
    unittest.main()