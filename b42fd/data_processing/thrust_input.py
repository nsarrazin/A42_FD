from math import *
import numpy as np
from ambiance import Atmosphere   #need to install this package

"""It returns the input data matrix for thrust.exe"""

gamma=1.4   #specific heat ratio
R=287       #gas constant
rho0=1.225  #sea level air density
p0=101325   #sea level pressure
T0=288.15   #sea level temperature
g0=9.81     #gravitational acceleration
lamb=-0.0065       #lapse rate

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

def Input(h,V,TAT,MFl,MFr, gamma,T0,lamb,g0,R,p0,rho0):
    "it returns a matrix of input data for computing thrust"
    
    thrust_input=np.zeros((len(h),5))    
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
        thrust_input[i,4]=MF2
        
    return thrust_input

"""======================================================================================
                             FOR FLIGHT DATA 
======================================================================================"""

#input data for stationary measurements Cl-Cd    
input0= Input(h,V,TAT,MFl,MFr, gamma,T0,lamb,g0,R,p0,rho0)
    
#stationary measurements elevator trim curve
h=np.array([18060, 18360, 18940, 18350, 18090, 17680, 18360])*ft
V=np.array([156, 147, 134, 168, 176, 186, 156])*kt
TAT=np.array([-10.2, -11.5, -13.5, -10.5, -9.5, -7.8, -11.2])+273.15
MFl=np.array([409,407,404,410,416,420, 480])*lbshr
MFr=np.array([470, 466, 455, 471, 448, 484, 469])*lbshr

#input thrust data elevator trim curve
input1= Input(h,V,TAT,MFl,MFr, gamma,T0,lamb,g0,R,p0,rho0)

"""======================================================================================
                             FOR REFERENCE DATA 
======================================================================================"""

#for stationary measurements Cl-Cd 
h=np.array([5010, 5020, 5020, 5030, 5020, 5110])*ft
V=np.array([249,221,192,163,130,118])*kt
TAT=np.array([12.5, 10.5, 8.8, 7.2,6,5.2])+273.15
MFl=np.array([798,673,561,463,443,474])*lbshr
MFr=np.array([813,682,579,484,467,499])*lbshr


input2=Input(h,V,TAT,MFl,MFr, gamma,T0,lamb,g0,R,p0,rho0)

#For elevator trim curve
h=np.array([6060, 6350, 6550, 6880, 6160, 5810, 5310])*ft
V=np.array([161, 150, 140, 130, 173, 179,192])*kt
TAT=np.array([5.5, 4.5, 3.5, 2.5, 5.0, 6.2, 8.2])+273.15
MFl=np.array([462,458,454,449,465,472, 482])*lbshr
MFr=np.array([486,482,477,473,489,496,505])*lbshr


input3=Input(h,V,TAT,MFl,MFr, gamma,T0,lamb,g0,R,p0,rho0)
