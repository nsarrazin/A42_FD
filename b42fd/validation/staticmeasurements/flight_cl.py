import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import c
from b42fd.data_processing.thrust_input import pressure, Mach, corrected_temp,sound_speed, true_airspeed
from b42fd.consts import gamma,T0,lamb,g0,R,p0,rho0
from b42fd.numerical_model.case import Case
from b42fd.validation.staticmeasurements.flight_data import h_m, V_ms, TAT_K, alpha_deg, fuelburnt_kg, mramp_kg

print("Flight test data Cl-a")

#fixed values
S = 30                  #[m^2]
mu = 3.178 * 10**-5 # dynamic viscosity of air (constant between 300-1000 K) kg m−1 s−1
M_ramp = 13491.5

#data
h_m = h_m[:6]                                               #altitude 
V_ms = V_ms[:6]                                             #indicated airspeed
TAT_K = TAT_K[:6]                                           #temp
alpha_deg = alpha_deg[:6]                                   #aoa deg
alpha_rad = np.radians(alpha_deg)
fuelburnt_kg = fuelburnt_kg[:6]                             #fuel burnt lbs

MFl=np.array([392,369,608,508,453,431])/7936.64             #fuel flow left engine (already converted to kg/s)
MFr=np.array([450, 378,668, 548, 488, 480])/7936.64

#density from ambiance package (not a standard package so install)
atmospheres = Atmosphere(h_m)
rho = atmospheres.density

#To find true airspeed
V_TAS=np.zeros(len(h_m))
Mlst = np.zeros(len(h_m))
Relst = np.zeros(len(h_m))

for i in range(len(h_m)):
    Vc_ms=V_ms[i]
    hp_m=h_m[i]
    Tm_K=TAT_K[i]
    p=pressure(hp_m, gamma,T0,lamb,g0,R,p0)
    M=Mach(Vc_ms,gamma, rho0,p0, p)
    print(M)
    T=corrected_temp(Tm_K,M,gamma)
    a=sound_speed(gamma,R,T)
    V_TAS[i]=true_airspeed(M,a)
    Mlst[i] = M
    Re = V_TAS[i]*c*rho[i]/mu
    Relst[i] = Re


Machmin = min(Mlst)
Machmax = max(Mlst)

Remin = min(Relst)
Remax = max(Relst)



#calculations
mass_1 = mramp_kg - fuelburnt_kg
CL = (mass_1*g0)/(0.5*rho*V_TAS**2*S)

CL_alpha_deg = (CL[-1] - CL[0])/(alpha_deg[-1]-alpha_deg[0])        # in 1/degree
CL_alpha_rad = (CL[-1] - CL[0])/(alpha_rad[-1]-alpha_rad[0])        # in 1/radian

#CL-alpha plot
title = 'CL-alpha curve clean configuration with Mach range: M =' + str(round(Machmin,2)) + ' to M =' + str(round(Machmax,2))
subtitle = 'Reynoldsnumber range: Re = ' + str(round(Remin/10**6,1)) + r' * $10^6$ to Re =' + str(round(Remax/10**6,1)) + r' * $10^6$'

s = sorted(zip(alpha_deg,CL))
alpha_degplt,CLplt = map(list, zip(*s))

plt.plot(alpha_degplt,CLplt)
plt.title(title,y=1.07,fontsize=16)
plt.suptitle(subtitle,y=0.92,fontsize=16)
plt.xlabel( r'$\alpha$ - angle of attack [deg]', fontsize = 14)
plt.ylabel(r'$C_L$ [-]', fontsize = 14)
plt.show()


#print CL_alpha
print('Cl_a =') 
print(CL_alpha_rad)     # 1.2486831027559244   [1/radian]