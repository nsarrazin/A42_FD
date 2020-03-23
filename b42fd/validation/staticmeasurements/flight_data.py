import numpy as np
from b42fd.numerical_model.Cit_par import M_ramp, c

#from 20200310_V2
h = np.array([18000,17990,18020,17990,18000,18010,18060,18360,18940,18350,18090,17680,18360,18369,18550])   #altitude 
V = np.array([161,131,222,200,182,114,156,147,134,168,176,186,156,156,156])                                 #indicated airspeed
TAT_C = np.array([-9.5,-11.5,-4.8,-6.8,-8.2,-12.8,-10.2,-11.5,-13.5,-10.5,-9.5,-7.8,-11.2,-11.2,-11.2])     #temp
alpha_deg = np.array([5,8.1,2.0,2.9,3.6,10.7,5.2,6.3,7.5,4.4,3.8,3.3,5.2,5.2,5.2])                          #aoa deg
fuelburnt_lbs = np.array([583,569,634,665,688,729,811,840,865,888,901,912,940,940,989])                     #fuel burnt lbs
mramp_lbs = M_ramp
MFl =  np.array([392, 369, 608, 508, 453, 431])/7936.64                                                           #fuel flow left engine (already converted to kg/s)
MFr =  np.array([450, 378, 668, 548, 488, 480])/7936.64

print(mramp_lbs)

#conversions 
h_m =   h*0.3048
V_ms  = V*0.5144
alpha_rad = np.radians(alpha_deg)
mramp_kg = mramp_lbs*0.45359237
fuelburnt_kg = fuelburnt_lbs*0.45359237
TAT_K = TAT_C+273.15
