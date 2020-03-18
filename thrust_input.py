from math import *
import numpy as np
from ambiance import Atmosphere

gamma=1.4   #specific heat ratio
R=287       #gas constant
rho0=1.225  #sea level air density
p0=101325   #sea level pressure
T0=288.15   #sea level temperature
g0=9.81     #gravitational acceleration
lamb=-0.0065       #bypass ratio

#coversions
ft=0.3048  #ft to m
kt=0.51444 #kts to m/s
lbshr=0.000125998 #lbs/hr to kg/s
#stationary measurement Cl-Cd
h=np.array([18000, 17990, 18020, 17990,18000, 18010])*ft
V=np.array([161,131, 222,200, 182, 114])*kt
TAT=np.array([-9.5,-11.5,-4.8,-6.8,-8.2,-12.8])+273.15
MFl=np.array([450, 378, 668, 548, 488, 480])*lbshr
MFr=np.array([538, 569, 634, 665, 688, 729])*lbshr
    
def Mach(Vc, hp, gamma, rho0,p0, p):
    """input: Vc: caliberated airspeed
              hp: pressure altitude
       output: a :speed of sound"""
    C=1+(gamma-1)/2/gamma*rho0/p0*Vc**2
    D=gamma/(gamma-1)
    M=np.sqrt(2/(gamma-1)*((1+p0/p*(C**D-1))**(1/C)-1))
    return M

def corrected_temp(Tm, M, gamma):
    """
    input:
    Tm: measured total air temperature (indicated as TAT in excel file)
    M: Mach number

    output: T: corrected temperature"""
    
    T=Tm/(1+(gamma-1)*M**2)
    return T


def pressure(hp, gamma,T0,lamb,g0,R,p0):
    """
    input:
    hp: pressure altitude

    output: T: pressure"""
    return p0*(1+lamb*hp/T0)**(-g0/lamb/R)

def temp_diff(T,T_isa):
    return T-T_isa

#input data for stationary measurements Cl-Cd
thrust_input=np.zeros((6,5))
hp=0
Vc=0
Tm=0

for i in range(len(h)):
    Vc=V[i]
    hp=h[i]
    Tm=TAT[i]
    MF1=MFl[i]
    MF2=MFr[i]
    p=pressure(hp, gamma,T0,lamb,g0,R,p0)
    M=Mach(Vc,hp, gamma, rho0,p0, p)
    T=corrected_temp(Tm,M,gamma)
    T_isa=Atmosphere(hp).temperature
    Temp_diff=temp_diff(T,T_isa)
    thrust_input[i,0]=hp
    thrust_input[i,1]=M
    thrust_input[i,2]=Temp_diff
    thrust_input[i,3]=MF1
    thrust_input[i,4]=MF1
    
print(thrust_input)    
    
#stationary measurements elevator trim curve
h=np.array([18060, 18360, 18940, 18350, 18090, 17680, 18360])*ft
V=np.array([156, 147, 134, 168, 176, 186, 156])*kt
TAT=np.array([-10.2, -11.5, -13.5, -10.5, -9.5, -7.8, -11.2])+273.15
MFl=np.array([409,407,404,410,416,420, 480])*lbshr
MFr=np.array([470, 466, 455, 471, 448, 484, 469])*lbshr

#input thrust data elevator trim curve
thrust_input1=np.zeros((7,5))

for i in range(len(h)):
    Vc=V[i]
    hp=h[i]
    Tm=TAT[i]
    MF1=MFl[i]
    MF2=MFr[i]
    p=pressure(hp, gamma,T0,lamb,g0,R,p0)
    M=Mach(Vc,hp, gamma, rho0,p0, p)
    T=corrected_temp(Tm,M,gamma)
    T_isa=Atmosphere(hp).temperature
    Temp_diff=temp_diff(T,T_isa)
    thrust_input1[i,0]=hp
    thrust_input1[i,1]=M
    thrust_input1[i,2]=Temp_diff
    thrust_input1[i,3]=MF1
    thrust_input1[i,4]=MF1
    
print(thrust_input1) 
