import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import M_ramp, c
from b42fd.data_processing.thrust_input import pressure, Mach, corrected_temp,sound_speed, true_airspeed
from b42fd.consts import gamma,T0,lamb,g0,R,p0, rho0
from b42fd.numerical_model.case import Case

print("Flight test data trim curve")

#from 20200310_V2
h = np.array([18060,18360,18940,18350,18090,17680,18360]) #altitude 
V = np.array([156,147,134,168,176,186,156])#indicated airspeed
TAT_C = np.array([-10.2,-11.5,-13.5,-10.5,-9.5,-7.8,-11.2]) #temp
alpha_deg = np.array([5.2,6.3,7.5,4.4,3.8,3.3,5.2])  #aoa deg
fuelburnt_lbs = np.array([811,840,865,888,901,912,940]) #fuel burnt lbs
de_deg = np.array([-0.3,-0.7,-1.2,0.1,0.4,0.7,-0.2]) 
mramp_lbs = M_ramp

#conversions 
h_m = h*0.3048
V_ms  = V*0.5144
alpha_rad = np.radians(alpha_deg)
mramp_kg = mramp_lbs*0.45359237
fuelburnt_kg = fuelburnt_lbs*0.45359237
TAT_K = TAT_C+273.15
de_rad_2=np.radians(de_deg)


#plot trim curve
s = sorted(zip(alpha_deg,de_deg))
alpha_deg,de_deg = map(list, zip(*s))

plt.plot(alpha_deg,de_deg,'-')
plt.xlabel('angle of attack [deg]')
plt.ylabel('elevator deflection [deg]')
plt.show()

de_da = -(max(de_deg)-min(de_deg))/(max(alpha_deg)-min(alpha_deg))        # [-]

print('de/dalpha =')
print(de_da)