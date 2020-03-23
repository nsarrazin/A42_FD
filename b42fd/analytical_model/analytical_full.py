#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 14:31:52 2020

@author: SUYIWEN
"""

from math import pi, sin, cos, sqrt, log
import numpy as np
from b42fd.numerical_model.Cit_par import *
from b42fd.analytical_model.time import TimeTool
from b42fd.helpers import load_data
from pathlib import Path

class Analytical_Eigenmotion:
    def __init__(self,motion_type,t=0, M_u=0):
        
        if motion_type =="short period motion":
            self.eigvalues, self.properties=self.get_eignval_spm(t, M_u) 
        
        if motion_type=="phugoid oscillation":
            self.eigvalues, self.properties=self.get_eigval_phu(t, M_u)
            
        if motion_type=="dutch roll":
            self.eigvalues, self.properties=self.get_eigval_dutchroll(t, M_u)
            
        if motion_type=="aperiodic roll":
            self.eigvalues, self.properties=self.get_eigval_ape_roll(t, M_u)
            
        if motion_type=="aperiodic spiral":
            self.eigvalues, self.properties=self.get_eigval_ape_spiral(t, M_u)
    
    @staticmethod
    def cal_eigenvalues(A,B,C,V,c):
        return (complex(-B/2/A*V/c,np.sqrt(4*A*C-B**2)/2/A*V/c),complex(-B/2/A*V/c,-np.sqrt(4*A*C-B**2)/2/A*V/c) )
   
    @staticmethod
    def eigenvalue_properties(eigenvalue):       
        
        eta=abs(eigenvalue.imag)     #imaginary part of eigenvalue
        zeta=abs(eigenvalue.real)   #real part of eigenvalue
        
        half_t = -log(1 / 2) /zeta

        if eigenvalue.imag == 0:
            P = None
            damping_ratio= None
            
        else:
            P = (2*pi)/eta
            damping_ratio=-eigenvalue.real/sqrt(eigenvalue.real**2+eigenvalue.imag**2)

        return damping_ratio, P, half_t.real
    
    def get_eignval_spm(self,t,M_u):
        params=TimeTool(t,M_u)
        muc, V0=params.muc, params.true_airspeed
        
        A=2*muc*KY2*(2*muc-CZadot)
        B=-2*muc*KY2*CZa-(2*muc+CZq)*Cmadot-(2*muc-CZadot)*Cmq
        C=CZa*Cmq-(2*muc+CZq)*Cma
        
        eigval=self.cal_eigenvalues(A, B, C,V0, c)[0]
        properties= self. eigenvalue_properties(eigval)
        
        return  eigval, properties
    
    def get_eigval_phu(self,t, M_u):
        
        params=TimeTool(t, M_u)
        muc, V0, CZ0=params.muc, params.true_airspeed, params.CZ0
        
        A=2*muc*(CZa*Cmq-2*muc*Cma)
        B=2*muc*(CXu*Cma-Cmu*CXa)+Cmq*(CZu*CXa-CXu*CZa)
        C=CZ0*(Cmu*CZa-Cma*CZu)
        
        eigval=self.cal_eigenvalues(A, B, C,V0, c)[0]
        properties= self.eigenvalue_properties(eigval)
        
        return eigval, properties
    
    def get_eigval_dutchroll(self, t, M_u):
        
        params=TimeTool(t, M_u)
        mub, V0=params.mub, params.true_airspeed
    
        A=8*mub**2*KZ2
        B=-2*mub*(Cnr+2*KZ2*CYb)
        C=4*mub*Cnb+CYb*Cnr
        
        eigval=self.cal_eigenvalues(A, B, C,V0, b)[0]
        properties= self.eigenvalue_properties(eigval)
        
        return eigval, properties
    
    def get_eigval_ape_roll(self,  t, M_u):
        params=TimeTool(t, M_u)
        mub, V0=params.mub, params.true_airspeed
        
        eig1=complex(Clp/(4*mub*KX2)*V0/b,0)
        
        return eig1, self.eigenvalue_properties(eig1)
    
    def get_eigval_ape_spiral(self, t, M_u):
        
        params=TimeTool(t, M_u)
        mub, V0, CL=params.mub, params.true_airspeed, params.CL
        eig2=complex((2*CL*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))*V0/b,0)
        
        return eig2, self.eigenvalue_properties(eig2)
 
if __name__ == "__main__":
    
    t_phugoid=53*60
    t_spm    =58*60
    t_dutchroll=60*60   #there is another dutch roll in yawning direction. not sure if we had to use it. Time=61*60
    
    t_ape_roll=57*60
    t_ape_spiral=62*60
    
    print("===========================FOR FLIGHT DATA========================")
    
    print("--------------------- Symmetric motions ---------------------")
    
    print("\n                    Short Period Motion")
    m_pax=np.array([95,102,89,82,66,81,69,85,96])    #passenger weights in kg
    M_e=9165*0.453592                     #empty aircraft weight in kg
    M_u=2640*0.453592                     #mass of fuel

    short_period= Analytical_Eigenmotion("short period motion", t=t_spm)
    print("eigenvalues are :", short_period.eigvalues)
    print("eigenvalues properties (damping ratio, period, halving time):", short_period.properties)
    
    print("\n                    Phugoid Motion")
    
    phugoid= Analytical_Eigenmotion("phugoid oscillation", t=t_phugoid)
    print("eigenvalues:", phugoid.eigvalues)
    print("eigenvalues properties (damping ratio, period, halving time):", phugoid.properties)
    
    print("\n-------------------- Asymmetric motions -------------")
    
    print("\n                     Dutch Roll")
    
    dutch_roll= Analytical_Eigenmotion("dutch roll", t=t_phugoid)
    print("eigenvalues for Dutch Roll:", dutch_roll.eigvalues)
    print("eigenvalues properties for Dutch Roll (damping ratio, period, halving time):", dutch_roll.properties)
    
    print("\n                     Aperiodic Roll")
    
    aperiodic_roll= Analytical_Eigenmotion("aperiodic roll", t=t_ape_roll)
    print("eigenvalues for aperiodic roll:", aperiodic_roll.eigvalues)
    print("eigenvalues properties for aperiodic roll(damping ratio, period, halving time):", aperiodic_roll.properties)
    
    print("\n                     Aperiodic Spiral")
     
    aperiodic_spiral=Analytical_Eigenmotion("aperiodic spiral", t=t_ape_spiral)
    print("eigenvalues :", aperiodic_spiral.eigvalues)
    print("eigenvalues properties (damping ratio, period, halving time):", aperiodic_spiral.properties)
    
    
    data=load_data("data/ref_data/ref_data.json")
    m_pax=np.array([95,92,74,66,61,75,78,86,68])
    M_e=9165*0.453592 
    M_u=4050*0.453592
    
    t_phugoid=53*60+57
    t_spm    =60*60+35
    t_dutchroll=60*60   #there is another dutch roll in yawning direction. not sure if we had to use it. Time=61*60
    
    t_ape_roll=59*60+10
    t_ape_spiral=62*60+20
    
    print("\n===========================FOR REFERENCE DATA========================")
    
    print("\n--------------------- Symmetric motions ---------------------")
    
    
    short_period= Analytical_Eigenmotion("short period motion", t=t_spm)
    print("eigenvalues are :", short_period.eigvalues)
    print("eigenvalues properties (damping ratio, period, halving time):", short_period.properties)
    
    print("\n                    Phugoid Motion")
    
    phugoid= Analytical_Eigenmotion("phugoid oscillation", t=t_phugoid)
    print("eigenvalues:", phugoid.eigvalues)
    print("eigenvalues properties (damping ratio, period, halving time):", phugoid.properties)
    
    print("\n-------------------- Asymmetric motions -------------")
    
    print("\n                     Dutch Roll")
    
    dutch_roll= Analytical_Eigenmotion("dutch roll", t=t_phugoid)
    print("eigenvalues for Dutch Roll:", dutch_roll.eigvalues)
    print("eigenvalues properties for Dutch Roll (damping ratio, period, halving time):", dutch_roll.properties)
    
    print("\n                     Aperiodic Roll")
    
    aperiodic_roll= Analytical_Eigenmotion("aperiodic roll", t=t_ape_roll)
    print("eigenvalues for aperiodic roll:", aperiodic_roll.eigvalues)
    print("eigenvalues properties for aperiodic roll(damping ratio, period, halving time):", aperiodic_roll.properties)
    
    print("\n                     Aperiodic Spiral")
     
    aperiodic_spiral=Analytical_Eigenmotion("aperiodic spiral", t=t_ape_spiral)
    print("eigenvalues :", aperiodic_spiral.eigvalues)
    print("eigenvalues properties (damping ratio, period, halving time):", aperiodic_spiral.properties)