#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:45:39 2020
 
@author=suyiwen
"""

from b42fd.validation.fuelmass import data
import numpy as np
from b42fd.numerical_model.Cit_par import *
from b42fd.helpers import load_data

M_e=9165*0.453592                     #empty aircraft weight in kg
class TimeTool:
    
    def __init__(self, data, t, M_u, m_pax,CLa, CD0, e):
        self.data=data    #flight data
        #sefl.data=load_data("data/ref_data/ref_data.json")
        
        self.time=data["time"]["data"]
        self.altitude=data["Dadc1_alt"]["data"]
        self.rh_FU=data["rh_engine_FU"]["data"]
        self.lf_FU=data["lh_engine_FU"]["data"]
        self.V_TAS=data['Dadc1_tas']["data"]
        self.alpha=data['vane_AOA']["data"]
        self.theta=data["Ahrs1_Pitch"]["data"]
                        
        if t !=0:
            self.altitude, self.true_airspeed, self.angle_of_attack, self.theta, self.weight, self.rho,self.mub, self.muc, self.CL, self.CD, self.CX0, self.CZ0, self.fuel_mass_used, self.t=self.get_parameters(t, M_u,m_pax, CLa, CD0, e)
        
    def get_parameters(self,t, M_u,m_pax,CLa, CD0, e):
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
        time=self.time
        W0=M_e+M_u+np.sum(m_pax)
        rh_FU=self.rh_FU
        lh_FU=self.lf_FU
        
        weight=np.zeros(len(rh_FU))
        m_fuel=np.zeros(len(rh_FU))
        
        #get an array of weight at all instance of time
        for i in range(len(rh_FU)):
                m_fuel[i]=rh_FU[i]+lh_FU[i]
                weight[i]=(W0-m_fuel[i])*g
        
        for idx, t_i in enumerate(time):
            if time[idx] <= t < time[idx+1]:
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

        return hp0, V0, alpha0, th0, W, rho, mub, muc, CL, CD, CX0, CZ0, m_fuel[idx], time[idx]
    
if __name__ == "__main__":
    
    m_pax=np.array([95,102,89,82,66,81,69,85,96])    #passenger weights in kg
    M_u=2640*0.453592                     #mass of fuel
    
    Cm_de= -1.1941458222011172
    Cma= -0.5402088243290768
    
    CLa=4.371485054942859
    CD0=0.016
    e=0.6
    
    #Flight Data 
    #all in seconds
    t_phugoid=53*60
    t_spm    =58*60
    t_dutchroll=60*60   #there is another dutch roll in yawning direction. not sure if we had to use it. Time=61*60
    
    t_ape_roll=57*60
    t_ape_spiral=62*60
        
    short_period=TimeTool(data,t=t_spm, M_u=M_u,m_pax=m_pax,CLa=CLa, CD0=CD0,  e=e)
    phugoid     =TimeTool(data,t=t_phugoid, M_u=M_u,m_pax=m_pax,CLa=CLa, CD0=CD0,  e=e)
    dutch_roll  =TimeTool(data,t=t_dutchroll, M_u=M_u,m_pax=m_pax,CLa=CLa, CD0=CD0,  e=e)
    aperiodic_roll =TimeTool(data,t=t_ape_roll, M_u=M_u,m_pax=m_pax,CLa=CLa, CD0=CD0,  e=e)
    aperiodic_spiral=TimeTool(data,t=t_ape_spiral, M_u=M_u,m_pax=m_pax, CLa=CLa, CD0=CD0,  e=e)
    
    print("\n---------------------FOR FLIGHT DATA----------------------------")
    
    print("for short period motion:", period.altitude, short_period.true_airspeed, short_period.angle_of_attack, short_period.theta, short_period.weight, short_period.rho, short_period.mub, short_period.muc, short_period.CL, short_period.CD, short_period.CX0, short_period.CZ0, short_period.fuel_mass_used)
    print("\nfor phugoid oscillation:", phugoid.altitude, phugoid.true_airspeed, phugoid.angle_of_attack, phugoid.theta, phugoid.weight, phugoid.rho, phugoid.mub, phugoid.muc, phugoid.CL, phugoid.CD, phugoid.CX0, phugoid.CZ0 )
    print("\nfor Dutch Roll:", dutch_roll.altitude, dutch_roll.true_airspeed,  dutch_roll.angle_of_attack,  dutch_roll.theta,  dutch_roll.weight,  dutch_roll.rho,  dutch_roll.mub,  dutch_roll.muc,  dutch_roll.CL,  dutch_roll.CD,  dutch_roll.CX0,  dutch_roll.CZ0)        
    print("\nfor aperiodic roll motion", aperiodic_roll.altitude, aperiodic_roll.true_airspeed,  aperiodic_roll.angle_of_attack,  aperiodic_roll.theta,  aperiodic_roll.weight,  aperiodic_roll.rho,  aperiodic_roll.mub,  aperiodic_roll.muc,  aperiodic_roll.CL,  aperiodic_roll.CD,  aperiodic_roll.CX0,  aperiodic_roll.CZ0)
    print("\nfor aperiodic spiral motion", aperiodic_spiral.altitude, aperiodic_spiral.true_airspeed,  aperiodic_spiral.angle_of_attack,  aperiodic_spiral.theta,  aperiodic_spiral.weight,  aperiodic_spiral.rho,  aperiodic_spiral.mub,  aperiodic_spiral.muc,  aperiodic_spiral.CL,  aperiodic_spiral.CD,  aperiodic_spiral.CX0,  aperiodic_spiral.CZ0)
    
    data=load_data("data/ref_data/ref_data.json")
    
    m_pax=np.array([95,92,74,66,61,75,78,86,68])
    M_e=9165*0.453592 
    M_u=4050*0.453592
    
    #stationary mesurements results 
    Cm_de= -1.0024724929977598
    Cma_ref= -0.4914080848028234
     
    CLa=4.662367336619402
    CD0=0.016
    e=0.88

    #Reference data
    t_phugoid=53*60+57
    t_spm    =60*60+35
    t_dutchroll=61*60+57   #there is another dutch roll in yawning direction. not sure if we had to use it. Time=62*60+47
    
    t_ape_roll=59*60+10
    t_ape_spiral=65*60+20
    
    short_period=TimeTool(data,t=t_spm, M_u=M_u,m_pax=m_pax, CLa=CLa, CD0=CD0,  e=e)
    phugoid     =TimeTool(data,t=t_phugoid, M_u=M_u,m_pax=m_pax, CLa=CLa, CD0=CD0,  e=e)
    dutch_roll  =TimeTool(data,t=t_dutchroll, M_u=M_u,m_pax=m_pax, CLa=CLa, CD0=CD0,  e=e)
    aperiodic_roll =TimeTool(data,t=t_ape_roll, M_u=M_u,m_pax=m_pax, CLa=CLa, CD0=CD0,  e=e)
    aperiodic_spiral=TimeTool(data,t=t_ape_spiral, M_u=M_u,m_pax=m_pax, CLa=CLa, CD0=CD0,  e=e)
    
    print("\n------------------FOR REFERENCE DATA----------------------------")
    
    print("\nfor short period motion:", short_period.altitude, short_period.true_airspeed, short_period.angle_of_attack, short_period.theta, short_period.weight, short_period.rho, short_period.mub, short_period.muc, short_period.CL, short_period.CD, short_period.CX0, short_period.CZ0)
    print("\nfor phugoid oscillation:", phugoid.altitude, phugoid.true_airspeed, phugoid.angle_of_attack, phugoid.theta, phugoid.weight, phugoid.rho, phugoid.mub, phugoid.muc, phugoid.CL, phugoid.CD, phugoid.CX0, phugoid.CZ0 )
    print("\nfor Dutch Roll:", dutch_roll.altitude, dutch_roll.true_airspeed,  dutch_roll.angle_of_attack,  dutch_roll.theta,  dutch_roll.weight,  dutch_roll.rho,  dutch_roll.mub,  dutch_roll.muc,  dutch_roll.CL,  dutch_roll.CD,  dutch_roll.CX0,  dutch_roll.CZ0)        
    print("\nfor aperiodic roll motion", aperiodic_roll.altitude, aperiodic_roll.true_airspeed,  aperiodic_roll.angle_of_attack,  aperiodic_roll.theta,  aperiodic_roll.weight,  aperiodic_roll.rho,  aperiodic_roll.mub,  aperiodic_roll.muc,  aperiodic_roll.CL,  aperiodic_roll.CD,  aperiodic_roll.CX0,  aperiodic_roll.CZ0)
    print("\nfor aperiodic spiral motion", aperiodic_spiral.altitude, aperiodic_spiral.true_airspeed,  aperiodic_spiral.angle_of_attack,  aperiodic_spiral.theta,  aperiodic_spiral.weight,  aperiodic_spiral.rho,  aperiodic_spiral.mub,  aperiodic_spiral.muc,  aperiodic_spiral.CL,  aperiodic_spiral.CD,  aperiodic_spiral.CX0,  aperiodic_spiral.CZ0)