import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import M_ramp, c
from b42fd.data_processing.thrust_input import pressure, Mach, corrected_temp,sound_speed, true_airspeed
from b42fd.consts import gamma,T0,lamb,g0,R,p0,rho0
from b42fd.numerical_model.case import Case

print("Flight test data Cl-a")

#fixed values
S = 30                  #[m^2]

#all data from 20200310_V2
h = np.array([18000,17990,18020,17990,18000,18010])                  #altitude 
V = np.array([161,131,222,200,182,114])                              #indicated airspeed
TAT_C = np.array([-9.5,-11.5,-4.8,-6.8,-8.2,-12.8])                  #temp
alpha_deg = np.array([5,8.1,2.0,2.9,3.6,10.7])                       #aoa deg
fuelburnt_lbs = np.array([583,569,634,665,688,729])                  #fuel burnt lbs
mramp_lbs = M_ramp

#conversions 
h_m = h*0.3048
V_ms  = V*0.5144
alpha_rad = np.radians(alpha_deg)
mramp_kg = mramp_lbs*0.45359237
fuelburnt_kg = fuelburnt_lbs*0.45359237
TAT_K = TAT_C+273.15


#density from ambiance package (not a standard package so install)
atmospheres = Atmosphere(h_m)
rho = atmospheres.density

#To find true airspeed
V_TAS=np.zeros(len(h))

for i in range(len(h)):
    Vc_ms=V_ms[i]
    hp_m=h_m[i]
    Tm_K=TAT_K[i]
    p=pressure(hp_m, gamma,T0,lamb,g0,R,p0)
    M=Mach(Vc_ms,hp_m, gamma, rho0,p0, p)
    T=corrected_temp(Tm_K,M,gamma)
    a=sound_speed(gamma,R,T)
    V_TAS[i]=true_airspeed(M,a)

#calculations
mass_1 = mramp_kg - fuelburnt_kg
CL = (mass_1*g0)/(0.5*rho*V_TAS**2*S)

CL_alpha_deg = (CL[-1] - CL[0])/(alpha_deg[-1]-alpha_deg[0])        # in 1/degree
CL_alpha_rad = (CL[-1] - CL[0])/(alpha_rad[-1]-alpha_rad[0])        # in 1/radian

#CL-alpha plot

plt.plot(alpha_deg,CL,'x')
plt.xlabel('angle of attack [deg]')
plt.ylabel('lift coefficient [-]')
plt.show()


#print CL_alpha
print('Cl_a =') 
print(CL_alpha_rad)     # 1.2486831027559244   [1/radian]