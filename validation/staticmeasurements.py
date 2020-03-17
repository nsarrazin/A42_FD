import numpy as np
import matplotlib.pyplot as plt
import isacalc as isa
#CL-CD and CL-a curves for flight test 20200306_V2

alpha_deg = np.array([1.4,2.3,3.2,5.2,8.3,10.2])
alpha_rad = np.radians(alpha_deg)

V_IAS_kts = np.array([250.,222.,193.,159.,129.,117.])
V_IAS_ms  = V_IAS_kts*0.5144

hp_ft = np.array([9000,9000,8990,8990,9010,8990])
hp_mt = hp_ft*0.3048

rho = np.array([0.933406,0.933406,0.933697,0.933697,0.933115,0.933697])
rho0 = 1.225        #kg/m^3

V_TAS_ms = V_IAS_ms * np.sqrt(rho0/rho)

mramp_lbs = 11623.123 #lbs
mramp_kg = mramp_lbs*0.45359237

fuelburnt_lbs = np.array([387,417,439,469,496,521])
fuelburnt_kg = fuelburnt_lbs*0.45359237

mass = mramp_kg - fuelburnt_kg
S = 30 #[m^2]

CL = mass/(0.5*rho*V_TAS_ms**2*S)

plt.plot(alpha_deg,CL)
plt.show()


#CL-CD and CL-a curves for flight test 20200306_V2

alpha_deg = np.array([1.6,2.4,3.7,5.6,8.5,10.4])
alpha_rad = np.radians(alpha_deg)

V_IAS_kts = np.array([251.,221.,190.,161.,134.,121.])
V_IAS_ms  = V_IAS_kts*0.5144

hp_ft = np.array([5030,5030,5020,5040,5040,5030])
hp_mt = hp_ft*0.3048

rho = np.array([0.933406,0.933406,0.933697,0.933697,0.933115,0.933697])
rho0 = 1.225        #kg/m^3

V_TAS_ms = V_IAS_ms * np.sqrt(rho0/rho)

mramp_lbs = 11623.123 #lbs
mramp_kg = mramp_lbs*0.45359237

fuelburnt_lbs = np.array([387,417,439,469,496,521])
fuelburnt_kg = fuelburnt_lbs*0.45359237

mass = mramp_kg - fuelburnt_kg
S = 30 #[m^2]

CL = mass/(0.5*rho*V_TAS_ms**2*S)

plt.plot(alpha_deg,CL)
plt.show()