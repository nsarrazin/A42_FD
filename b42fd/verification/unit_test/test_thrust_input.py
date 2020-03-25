#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:28:35 2020

@author: suyiwen
"""

import unittest
from b42fd.data_processing.thrust_input import *
from ambiance import Atmosphere
class TestAnalytical(unittest.TestCase):
    
    def test_Mach_zeroVc(self):
        
        """Results are verified by calculations by hand: when velocity is zero"""
        
        input_data={"Vc": 0,  "gamma":1.4, "rho0": 1.125, "p0": 101325, "p": 90000}
        hand_M=0
        
        self.assertEqual( Mach(input_data["Vc"], input_data["gamma"], input_data["rho0"], input_data["p0"], input_data["p"]), hand_M)
        
    def test_Mach_negativeVc(self):
        
        """Results are verified by calculations by hand: when velocity is negative"""
        
        input_data={"Vc": -10 , "gamma":3, "rho0": 1, "p0": 1, "p": 1}
        hand_M=5.773502692
        self.assertAlmostEqual( Mach(input_data["Vc"], input_data["gamma"], input_data["rho0"], input_data["p0"], input_data["p"]), hand_M, places=5)
    
    def test_Mach_largeVc(self):
        
        input_data={"Vc": 1000, "gamma":2, "rho0": 1, "p0": 1, "p": 1000}
        hand_M=125.73564201434377
        self.assertAlmostEqual( Mach(input_data["Vc"], input_data["gamma"], input_data["rho0"], input_data["p0"], input_data["p"]), hand_M, )
     
    def test_pressure_sealevel(self):
        input_data={"hp": 0, "gamma":1.4, "T0": 288.15, "lamb": -0.0065, "p0": 101325, "g0": 9.81, "R": 288.17,}
        test_result=101325
        self.assertEqual( pressure(input_data["hp"], input_data["gamma"], input_data["T0"], input_data["lamb"], input_data["g0"],input_data["R"],input_data["p0"]), test_result)

    def test_presssure_large_altitude(self):
        
        "compare results with another isa calculator"
        input_data={"hp": 11000, "gamma":1.4, "T0": 288.15, "p0": 101325, "lamb": -0.0065, "R": 287, "g0": 9.81}
        hand_result=22614.2066866830
        self.assertAlmostEqual( pressure(input_data["hp"], input_data["gamma"], input_data["T0"], input_data["lamb"], input_data["g0"],input_data["R"],input_data["p0"]), hand_result)
    
    def test_pressure_negativealtitude(self):
        
        "test for negative altitude"
        input_data={"hp": -10, "gamma":1.4, "T0": 288.15, "p0": 101325, "lamb": -0.0065, "R": 287, "g0": 9.81}
        hand_result=print("Pressure can not be determined for negative altitudes")
        self.assertEqual( pressure(input_data["hp"], input_data["gamma"], input_data["T0"], input_data["lamb"], input_data["g0"],input_data["R"],input_data["p0"]), hand_result)
    
    def test_sound_speed(self):
        input_data={"gamma":1, "R":1, "T":1}
        hand_result=1
        self.assertEqual(sound_speed(input_data["gamma"],input_data["R"], input_data["T"]),hand_result)
        
    def test_sound_speed_nagativeT(self):
        input_data={"gamma":1, "R":1, "T":-1}
        hand_result=print("Sound speed can not be determined for negative temperatures")
        self.assertEqual(sound_speed(input_data["gamma"],input_data["R"], input_data["T"]),hand_result)
   
    def test_sound_speed_0(self):
        input_data={"gamma":1, "R":1, "T":0}
        hand_result=0
        self.assertEqual(sound_speed(input_data["gamma"],input_data["R"], input_data["T"]),hand_result)
        
    def test_corrected_temp_0(self):
        input_data={"gamma":1, "Tm": 0, "M":1}
        hand_result=0
        self.assertEqual(corrected_temp(input_data["Tm"],input_data["M"], input_data["gamma"]),hand_result) 
    
    def test_corrected_temp_negative(self):
        input_data={"gamma":1, "Tm": -10, "M":1}
        hand_result= -10
        self.assertEqual(corrected_temp(input_data["Tm"],input_data["M"], input_data["gamma"]),hand_result) 
        
    def test_corrected_temp_large(self):
        input_data={"gamma":1, "Tm": 10000, "M":0.5}
        hand_result= 10000
        self.assertEqual(corrected_temp(input_data["Tm"],input_data["M"], input_data["gamma"]),hand_result) 
    
    def test_true_airspeed_0M(self):
        "test when mach is zero"
        input_data={"a":1, "M":0}
        hand_result= 0
        self.assertEqual(true_airspeed(input_data["a"],input_data["M"]),hand_result) 
        
    def test_true_airspeed_0a(self):
        "test when speed of sound is zero"
        input_data={"a":0, "M":1}
        hand_result= 0
        self.assertEqual(true_airspeed(input_data["a"],input_data["M"]),hand_result) 
    
    def test_true_airspeed_negativeMa(self):
        
        "test when speed of sound and mach are negative"
        input_data={"a":-1, "M":-1}
        hand_result= None
        self.assertEqual(true_airspeed(input_data["a"],input_data["M"]),hand_result) 
    
    def test_true_airspeed_negativeM(self):
        
        "test when speed of sound and mach are negative"
        input_data={"a":1, "M":-1}
        hand_result= None
        self.assertEqual(true_airspeed(input_data["a"],input_data["M"]),hand_result)  
    
    def test_true_airspeed_negativea(self):
        
        "test when speed of sound and mach are negative"
        input_data={"a":-1, "M":1}
        hand_result= None
        self.assertEqual(true_airspeed(input_data["a"],input_data["M"]),hand_result)      
        
if __name__ == '__main__':
    unittest.main()
        