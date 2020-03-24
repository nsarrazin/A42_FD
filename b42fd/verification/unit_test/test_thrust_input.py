#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:28:35 2020

@author: suyiwen
"""

import unittest
from b42fd.data_processing.thrust_input import *

class TestAnalytical(unittest.TestCase):
    
    def test_Mach_zeroVc(self):
        
        """Results are verified by calculations by hand: when velocity is zero"""
        
        input_data={"Vc": 0,  "gamma":1.4, "rho0": 1.125, "p0": 101325, "p": 90000}
        hand_M=0
        
        self.assertEqual( Mach(input_data["Vc"], input_data["gamma"], input_data["rho0"], input_data["p0"], input_data["p"]), hand_M)
        
    def test_Mach_negativeVc(self):
        
        input_data={"Vc": -10 , "gamma":3, "rho0": 1, "p0": 1, "p": 1}
        hand_M=5.773502692
        self.assertAlmostEqual( Mach(input_data["Vc"], input_data["gamma"], input_data["rho0"], input_data["p0"], input_data["p"]), hand_M, places=7)
    
    def test_Mach_largeVc(self):
        
        input_data={"Vc": 1000, "gamma":2, "rho0": 1, "p0": 1, "p": 1000}
        hand_M=17.72843277
        self.assertAlmostEqual( Mach(input_data["Vc"], input_data["gamma"], input_data["rho0"], input_data["p0"], input_data["p"]), hand_M, )
     
    def test_pressure_ambiance(self):
        input_data={"hp": 0, "gamma":1.4, "T0": 288.15, "p0": 101325, "lamb": 9.81, "R": 288.17, "p0": 101325 }
        test_result=1.225
        self.assertEqual( pressure(hp, gamma, T0, lamb, g0, R, p0), test_result)
        
        
if __name__ == '__main__':
    unittest.main()
        