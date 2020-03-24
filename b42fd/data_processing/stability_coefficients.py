#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 09:33:57 2020

@author: suyi
"""
import numpy as np
from b42fd.analytical_model.time import TimeTool
from b42fd.validation.fuelmass import data

def get_cg_shift(t_start, t_end, time, fuel_mass0, m_pax, x_pax_cg_old, x_pax_cg_new ):
        
    #Weight & Balance Values of Aircraft for CoG Calculations
    M_e=9165    #lbs
    M_e_arm= 291.65 #inches
    
    #taken from the assignment reader
    pax_arm= np.array([131,131,170, 214,214,251,251, 288, 288])
    
    #from figure E2.4 of the assignment reader
    fuel_left=np.linspace(100, 4900, num=49, endpoint=True) #pounds
    
    #moment in pound inches
    moment=np.array([298.16,591.18,879.08,1165.42, 1448.4,1732.53, 2014.8, 2298.84, 2581.92, 2866.3,3150.18, 3434.52, 3718.52, 4003.32, 4287.76, 4572.24, 4856.56, 5141.16, 5425.64, 5709.9, 5994.04, 6278.47, 6562.82, 6846.96, 7131, 7415.33, 7699.6, 7984.34, 8269.06, 8554.05, 8839.04, 9124.8, 9410.62, 9696.97, 9983.40, 10270.08, 10556.84, 10843.87, 11131.0, 11418.2, 11705.5, 11993.31, 12281.18, 12569.04, 12856.86, 13144.73, 13432.48, 13720.56, 14008.46 ])

    m_shift     =m_pax[8]   #kg

    #Aircraft Geometry
    """x_LEMAC = 261.56*2.54E-2
    c = 2.0569	"""
    M_u_kg=fuel_mass0*0.453592
    
    #Be careful with what data is being used. 
    FU1 = TimeTool(data,t_start,M_u_kg,m_pax, CLa, CD0,  e).fuel_mass_used 
    fuel1 = fuel_mass0-FU1*2.20462
    
    FU2= TimeTool(data,t_end,M_u_kg, m_pax, CLa, CD0,  e).fuel_mass_used
    fuel2=fuel_mass0-FU2*2.20462
    
    #get the index of the start time
    for idx_start, t_i in enumerate(time):
        if time[idx_start]< t_start <= time[idx_start+1]:
            break
    
    #get the index of the end time
    for idx_end, t_i in enumerate(time):
        if time[idx_end]< t_end <= time[idx_end+1]:
            break
        
    fuel_moment1=0
    fuel_moment2=0
    
    for i, x_i in enumerate(fuel_left):
        if fuel_left[i] <= fuel1 < fuel_left[i + 1]:
            fuel_moment1=(moment[i+1]-moment[i])/(fuel_left[i+1]-fuel_left[i])*fuel1 +moment[i]
        if fuel_left[i] <= fuel2 < fuel_left[i+1]:
            fuel_moment2 =(moment[i+1]-moment[i])/(fuel_left[i+1]-fuel_left[i])*fuel2 +moment[i]
        
    # Calculate CoG for Passenger Shift
    x_cg_old = (M_e * M_e_arm+ np.dot(m_pax,pax_arm)+ fuel_moment1) / (M_e + fuel1+sum(m_pax))         #in inches
    x_cg_new= (M_e * M_e_arm+ np.dot(m_pax,pax_arm)+ fuel_moment2 - (x_pax_cg_old- x_pax_cg_new)*m_shift) / (M_e + fuel2 + sum(m_pax))   #in inches

    return (x_cg_new-x_cg_old)*0.0254

def get_Cm_de(de, Cn, delta_cg, c):
    """input:     
    de: array of elevator deflections in radians
    Cn: normal force coefficient
    delta_cg: cg shift
    c: mac"""
    delta_de=de[1]-de[0]
    return -1/delta_de*Cn*delta_cg/c

def get_Cm_a(de, a, Cm_de):
    delta_de=min(de)-max(de)
    delta_a=max(a)-min(a)
    
    return -delta_de/delta_a*Cm_de

#wing geometry 
c = 2.0569

m_pax=np.array([95,102,89,82,66,81,69,85,96])*2.20462
t_start=49*60
t_end=52*60
fuel_mass0=2640

x_pax_cg_old=288  #inches
x_pax_cg_new=131  #inches

time=data["time"]["data"]

#stationary mesurements results 
Cm_de= -1.1941458222011172
Cma= -0.5402088243290768

CLa=4.371485054942859
CD0=0.016
e=0.6

delta_cg_flight=get_cg_shift(t_start, t_end, time, fuel_mass0, m_pax, x_pax_cg_old, x_pax_cg_new )
print("-----------------FLIGHT DATA--------------------")
print("\n shift in cg location in meters:", delta_cg_flight)

M_u_kg=fuel_mass0*0.453592

#to determine Cm_de using cg shift data 
de_1=np.radians([-0.2, -0.8])
Cn=TimeTool(data,t_start,M_u_kg, m_pax/2.20462, CLa, CD0,  e).CL
Cm_de= get_Cm_de(de_1, Cn, delta_cg_flight, c)

print("\n Cm_de:", Cm_de)

#to determine Cma using the stationary trim curve data
a=np.radians([5.2, 6.3, 7.5, 4.4, 3.8, 3.3, 5.2])
de=np.radians([-0.3, -0.7, -1.2, 0.1, 0.4, 0.7, -0.2])

Cm_a=get_Cm_a(de, a, Cm_de)

print("\n Cm_a:", Cm_a)
"""
print("\nweight     ", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).weight)
print("rho          ", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).rho)
print("altitude     ", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).altitude)
print("true airspeed", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).true_airspeed)
print("Cn           ", Cn)"""


from b42fd.helpers import load_data


data=load_data("data/ref_data/ref_data.json")

m_pax=np.array([95,92,74,66,61,75,78,86,68])*2.20462
time=data["time"]["data"]

t_start=51*60+2
t_end=52*60+46
fuel_mass0=4050

x_pax_cg_old=288  #inches
x_pax_cg_new=134  #inches

M_u_kg=fuel_mass0*0.453592

#stationary mesurements results 
Cm_de= -1.0024724929977598
Cma_ref= -0.4914080848028234
 
CLa=4.662367336619402
CD0=0.016
e=0.88

delta_cg_ref=get_cg_shift(t_start, t_end, time, fuel_mass0, m_pax, x_pax_cg_old, x_pax_cg_new)

print("\n-----------------REF DATA------------------------")
print("\nshift in cg location in meters:", delta_cg_ref)

#to determine Cm_de using cg shift data 
de_2=np.radians([-0, -0.5])
Cn=TimeTool(data, t_start,M_u_kg, m_pax/2.20462,CLa, CD0,  e).CL
Cm_de= get_Cm_de(de_2, Cn, delta_cg_ref, c)

#to determine Cma using the stationary trim curve data
a=np.radians([5.3, 6.3, 7.3, 8.5, 4.5, 4.1, 3.4])
de=np.radians([0, -0.4, -0.9, -1.5, 0.4, 0.6, 1])

Cm_a=get_Cm_a(de, a, Cm_de)


Cm_a=get_Cm_a(de, a, Cm_de)

print("\n Cm_de:", Cm_de)
print("\n Cm_a:", Cm_a)
"""

print("\nweight     ", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).weight)
print("rho          ", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).rho)
print("altitude     ", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).altitude)
print("true airspeed", TimeTool(data,t_start,M_u_kg, m_pax, CLa=CLa, CD0=CD0,  e=e).true_airspeed)
print("Cn           ", Cn)"""



