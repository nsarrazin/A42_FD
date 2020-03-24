import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import M_ramp, c
from b42fd.data_processing.thrust_input import pressure, Mach, corrected_temp,sound_speed, true_airspeed
from b42fd.consts import gamma,T0,lamb,g0,R,p0, rho0
from b42fd.numerical_model.case import Case
from b42fd.validation.staticmeasurements.flight_cl import CL
from b42fd.validation.staticmeasurements.flight_data import h_m, V_ms, TAT_K, alpha_deg, fuelburnt_kg

print("Flight test data CD")

#fixed values
b = 15.911
S = 30.
A = b**2/S
mu = 3.178 * 10**-5 # dynamic viscosity of air (constant between 300-1000 K) kg m−1 s−1

#from 20200310_V2
h_m = h_m[:6]                                               #altitude 
V_ms = V_ms[:6]                                             #indicated airspeed
TAT_K = TAT_K[:6]                                           #temp
alpha_deg = alpha_deg[:6]                                   #aoa deg
fuelburnt_kg = fuelburnt_kg[:6]                             #fuel burnt lbs
MFl=np.array([392,369,608,508,453,431])/7936.64             #fuel flow left engine (already converted to kg/s)
MFr=np.array([450, 378,668, 548, 488, 480])/7936.64

mramp_lbs = M_ramp


#density from ambiance package
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
    M=Mach(Vc_ms, gamma, rho0,p0, p)
    T=corrected_temp(Tm_K,M,gamma)
    a=sound_speed(gamma,R,T)
    V_TAS[i]=true_airspeed(M,a)
    Mlst[i] = M
    Re = V_TAS[i]*c*rho[i]/mu
    Relst[i] = Re


thrusts = []
for i in range(6):
    case = Case(h_m[i], V_ms[i], TAT_K[i], MFl[i], MFr[i])
    thrusts.append(case.thrust)

CD = thrusts/(0.5*rho[0:6]*V_TAS**2*S)

plt.plot(alpha_deg[:6],CD,'x')
plt.xlabel('alpha [deg]')
plt.ylabel('CD [-]')
plt.show()

Machmin = min(Mlst)
Machmax = max(Mlst)

Remin = min(Relst)
Remax = max(Relst)

title = 'CL-CD curve with Mach range: M =' + str(round(Machmin,2)) + ' to M =' + str(round(Machmax,2))
subtitle = 'Reynoldsnumber range: Re = ' + str(round(Remin/10**6,1)) + ' * 10^6 to Re =' + str(round(Remax/10**6,1)) + ' * 10^6'

plt.plot(CL,CD,'x',markersize = 10)
plt.title(title,fontsize = 16, y=1.07)
plt.suptitle(subtitle, fontsize = 16, y=0.92)
plt.xlabel('CL [-]',fontsize = 14)
plt.ylabel('CD [-]', fontsize = 14)
plt.show()


#vary these to see get them right

CD0 = 0.04 
e = 0.8

CLsquared = CL**2
CDtheoretical = CD0 + CLsquared/(np.pi*A*e)

plt.plot(CLsquared, CDtheoretical, 'o-')
plt.plot(CLsquared, CD,'x')
plt.xlabel('CL^2 [-]')
plt.ylabel('CD')
plt.show()