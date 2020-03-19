import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

#CL-CD and CL-a curves for flight test 20200306_V2

#these change per test flight
alpha_deg = np.array([1.6,2.4,3.7,5.6,8.5,10.4])                #degrees
V_IAS_kts = np.array([251.,221.,190.,161.,134.,121.])           #kts
hp_ft = np.array([5030,5030,5020,5040,5040,5030])               #ft
fuelburnt_lbs = np.array([387,417,439,469,496,521])             #lbs

#this one should change, but took the 'fixed' value from the report
mramp_lbs = 11623.123   #lbs

#fixed values
rho0 = 1.225            #kg/m^3
S = 30                  #[m^2]

#conversions
alpha_rad = np.radians(alpha_deg)
V_IAS_ms  = V_IAS_kts*0.5144
hp_mt = hp_ft*0.3048
mramp_kg = mramp_lbs*0.45359237
fuelburnt_kg = fuelburnt_lbs*0.45359237

#density from ambiance package (not a standard package so install)
atmospheres = Atmosphere(hp_mt)
rho = atmospheres.density

#calculations
V_TAS_ms = V_IAS_ms * np.sqrt(rho0/rho)
mass = mramp_kg - fuelburnt_kg

CL = mass/(0.5*rho*V_TAS_ms**2*S)

CL_alpha_deg = (CL[-1] - CL[0])/(alpha_deg[-1]-alpha_deg[0])        # in 1/degree
CL_alpha_rad = (CL[-1] - CL[0])/(alpha_rad[-1]-alpha_rad[0])        # in 1/radian

#CL-alpha plot

plt.plot(alpha_deg,CL,'x-')
plt.xlabel('angle of attack [deg]')
plt.ylabel('lift coefficient [-]')
plt.show()

#print CL_alpha
print('Cl_a =') 
print(CL_alpha_rad)     # 0.3522578695587925   [1/radian]



#20200306_V2

de_deg = np.array([0,-0.3,-0.7,-1.3,0.4,0.7,1])
alpha_deg_2 = np.array([5.3,6.1,6.8,8.3,4.2,3.8,3.3]) 

plt.plot(alpha_deg_2,de_deg,'x')
plt.xlabel('angle of attack [deg]')
plt.ylabel('elevator deflection [deg]')
plt.show()

de_da = (max(de_deg)-min(de_deg))/(max(alpha_deg_2)-min(alpha_deg_2))        #0.46 [-]

print('de/dalpha =')
print(de_da)