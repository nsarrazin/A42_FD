#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:45:39 2020
 
@author=suyiwen
"""

from b42fd.validation.fuelmass import data
import numpy as np

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
        
        #to be changed
        self.e=0.8
        self.CD0=0.04
        self.CLa=5.084
        self.Cma=-0.5626
        self.Cmde=1.1642
        
        
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
        # Constant values concerning atmosphere and gravity
        rho0   = 1.2250          # air density at sea level [kg/m^3] 
        lambd = -0.0065         # temperature gradient in ISA [K/m]
        Temp0  = 288.15          # temperature at sea level in ISA [K]
        R      = 287.05          # specific gas constant [m^2/sec^2K]
        g      = 9.81            # [m/sec^2] (gravity constant)
        
        # Aircraft geometry
        S      = 30.00	         # wing area [m^2]
        c      = 2.0569	         # mean aerodynamic cord [m]
        b      = 15.911	         # wing span [m]
        A      = b ** 2 / S      # wing aspect ratio [ ]
        
        e=self.e
        CD0=self.CD0
        CLa=self.CLa
        
        time=self.time
        W0=self.M_e+self.M_u+np.sum(self.m_pax)
        rh_FU=self.rh_FU
        lh_FU=self.lf_FU
        
        weight=np.zeros(len(rh_FU))
        
        #get an array of weight at all instance of time
        for i in range(len(rh_FU)):
                m_fuel=rh_FU[i]+lh_FU[i]
                weight[i]=(W0-m_fuel)*g
        
        print(len(time))
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
        
        mub=m/(rho*V0*S*b)
        muc=m/(rho*V0*S*c)
        
        CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
        CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]
        
        #Stability derivaties
        CX0    = W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
        CZ0    = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
        
        return hp0, V0, alpha0, th0, W, rho, mub, muc, CL, CD, CX0, CZ0
    
if __name__ == "__main__":
    eig=TimeTool(t=2000)
    
    print(eig.altitude, eig.true_airspeed, eig.angle_of_attack, eig.theta, eig.weight, eig.rho, eig.mub, eig.muc, eig.CL, eig.CD, eig.CX0, eig.CZ0)
            
    
