#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:45:39 2020
 
@author=suyiwen
"""

from b42fd.validation.fuelmass import data
import numpy as np
from b42fd.numerical_model.Cit_par import *

#from b42fd.helpers import load_data
#from pathlib import Path

class TimeTool:
    
    def __init__(self, t):
        self.data=data    #flight data
        #sefl.data=load_data("data/ref_data/ref_data.json")
        
        self.time=data["time"]["data"]
        self.altitude=data["Dadc1_alt"]["data"]
        self.rh_FU=data["rh_engine_FU"]["data"]
        self.lf_FU=data["lh_engine_FU"]["data"]
        self.V_TAS=data['Dadc1_tas']["data"]
        self.alpha=data['vane_AOA']["data"]
        self.theta=data["Ahrs1_Pitch"]["data"]
        
        self.m_pax=np.array([95,92,74,66,61,75,78,86,68])    #passenger weights in kg
        self.M_e=9165*0.453592                     #empty aircraft weight in kg
        self.M_u=4050*0.453592                     #mass of fuel
                        
        if t !=0:
            self.altitude, self.true_airspeed, self.angle_of_attack, self.theta, self.weight, self.rho,self.mub, self.muc, self.CL, self.CD, self.CX0, self.CZ0 =self.get_flight_conditions(t)
        
    def get_flight_conditions(self,t):
        """
        returns all parameters that are depended on time

        Parameters
        ----------
        t : time in seconds 

        Returns
        -------
        hp0 : pressure altiture in the stationary flight condation in meters
        V0 : true airspeed in the stationary flight condition in m/s
        alpha0 : angle of attack in the stationary flight condition in rad
        th0 : pitch angle in the stationary flight condition in rad
        W : current weight at time t in kg
        rho : air density in kg/m^3
        mub : mass parameter 
        muc : mass parameter
        CD : Lift coefficients 
        CL : Drag coefficient 

        """
        """# Constant values concerning atmosphere and gravity
        rho0   = 1.2250          # air density at sea level [kg/m^3] 
        lambd = -0.0065         # temperature gradient in ISA [K/m]
        Temp0  = 288.15          # temperature at sea level in ISA [K]
        R      = 287.05          # specific gas constant [m^2/sec^2K]
        g      = 9.81            # [m/sec^2] (gravity constant)
        
        # Aircraft geometry
        S      = 30.00	         # wing area [m^2]
        c      = 2.0569	         # mean aerodynamic cord [m]
        b      = 15.911	         # wing span [m]
        A      = b ** 2 / S      # wing aspect ratio [ ]"""
        
        time=self.time
        W0=self.M_e+self.M_u+np.sum(self.m_pax)
        rh_FU=self.rh_FU
        lh_FU=self.lf_FU
        
        weight=np.zeros(len(rh_FU))
        
        #get an array of weight at all instance of time
        for i in range(len(rh_FU)):
                m_fuel=rh_FU[i]+lh_FU[i]
                weight[i]=(W0-m_fuel)*g
        
        for idx, t_i in enumerate(time):
            if time[idx] < t <= time[idx+1]:
                break
        
        hp0=self.altitude[idx]
        th0=self.theta[idx]
    
        V0=self.V_TAS[idx]
        alpha0=self.alpha[idx]
        
        W=weight[idx]
        
        rho = rho0 * pow( ((1+(lambd * hp0 / Temp0))), (-((g / (lambd*R)) + 1)))   
        
        m=W/g
        
        mub=m/(rho*S*b)
        muc=m/(rho*S*c)
        
        CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
        CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]
        
        #Stability derivaties
        CX0    = W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
        CZ0    = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
        
        return hp0, V0, alpha0, th0, W, rho, mub, muc, CL, CD, CX0, CZ0
    
if __name__ == "__main__":
    
    #Flight Data 
    #all in seconds
    t_phugoid=53*60
    t_spm    =58*60
    t_dutchroll=60*60   #there is another dutch roll in yawning direction. not sure if we had to use it. Time=61*60
    
    t_ape_roll=57*60
    t_ape_spiral=62*60
    
    """Reference data
    t_phugoid=53*60+57
    t_spm    =60*60+35
    t_dutchroll=60*60   #there is another dutch roll in yawning direction. not sure if we had to use it. Time=61*60
    
    t_ape_roll=59*60+10
    t_ape_spiral=62*60+20
    """
    
    short_period=TimeTool(t=t_spm)
    phugoid     =TimeTool(t=t_phugoid)
    dutch_roll  =TimeTool(t=t_dutchroll)
    aperiodic_roll =TimeTool(t=t_ape_roll)
    aperiodic_spiral=TimeTool(t=t_ape_spiral)
    
    print("for short period motion:", short_period.altitude, short_period.true_airspeed, short_period.angle_of_attack, short_period.theta, short_period.weight, short_period.rho, short_period.mub, short_period.muc, short_period.CL, short_period.CD, short_period.CX0, short_period.CZ0)
    print("for phugoid oscillation:", phugoid.altitude, phugoid.true_airspeed, phugoid.angle_of_attack, phugoid.theta, phugoid.weight, phugoid.rho, phugoid.mub, phugoid.muc, phugoid.CL, phugoid.CD, phugoid.CX0, phugoid.CZ0 )
    print("for Dutch Roll:", dutch_roll.altitude, dutch_roll.true_airspeed,  dutch_roll.angle_of_attack,  dutch_roll.theta,  dutch_roll.weight,  dutch_roll.rho,  dutch_roll.mub,  dutch_roll.muc,  dutch_roll.CL,  dutch_roll.CD,  dutch_roll.CX0,  dutch_roll.CZ0)        
    print("for aperiodic roll motion", aperiodic_roll.altitude, aperiodic_roll.true_airspeed,  aperiodic_roll.angle_of_attack,  aperiodic_roll.theta,  aperiodic_roll.weight,  aperiodic_roll.rho,  aperiodic_roll.mub,  aperiodic_roll.muc,  aperiodic_roll.CL,  aperiodic_roll.CD,  aperiodic_roll.CX0,  aperiodic_roll.CZ0)
    print("for aperiodic roll motion", aperiodic_spiral.altitude, aperiodic_spiral.true_airspeed,  aperiodic_spiral.angle_of_attack,  aperiodic_spiral.theta,  aperiodic_spiral.weight,  aperiodic_spiral.rho,  aperiodic_spiral.mub,  aperiodic_spiral.muc,  aperiodic_spiral.CL,  aperiodic_spiral.CD,  aperiodic_spiral.CX0,  aperiodic_spiral.CZ0)