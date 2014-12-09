# -*- coding: utf-8 -*-
"""
Created on Sat May 18 21:25:57 2013

@author: Will
"""

import unittest
import acwires
from acmconstants import MU_0
from math import pi, sqrt

class TestHThinWireFunctions(unittest.TestCase):

    def setUp(self):
        self.wirelength = 1000
        self.testdist = .1
        self.current = 1
        self.hwire = acwires.HThinWire('test hthinwire', 
                                   self.wirelength, 
                                   .001, 
                                   .001, 
                                   self.current, 
                                   1, 
                                   0, 
                                   0, 
                                   0)

    def test_bfield(self):
        bfield = self.hwire.bfieldcalc(self.wirelength / 2, 0, self.testdist)
        truefield = (0, -1*((MU_0 * self.current * self.wirelength) 
            / (4 * pi * self.testdist * 
            sqrt(self.testdist**2 + 0.25*self.wirelength**2))), 0)
        self.assertAlmostEqual(bfield[1], truefield[1], 7,'True field: %f\nSim Field: %f'%(truefield[1], bfield[1]))
        
if __name__ == '__main__':
    unittest.main()