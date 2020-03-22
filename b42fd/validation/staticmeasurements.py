import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import M_ramp, c

#CL-CD and CL-a curves for flight test 20200310_V2

#these change per test flight
alpha_deg = np.array([5,8.1,2.0,2.9,3.6,10.7])                #degrees
V_IAS_kts = np.array([161,131,222,200,182,114])           #kts
hp_ft = np.array([18000,17990,18020,17990,18000,18010])               #ft
fuelburnt_lbs = np.array([583,569,634,665,688,729])             #lbs

#this one should change, but took the 'fixed' value from the report
mramp_lbs = M_ramp   #lbs, from flight data 20200310

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



#20200310_V2

de_deg = np.array([-0.3,-0.7,-1.2,0.1,0.4,0.7,-0.2])
alpha_deg_2 = np.array([5.2,6.3,7.5,4.4,3.8,3.3,5.2]) 

plt.plot(alpha_deg_2,de_deg,'x')
plt.xlabel('angle of attack [deg]')
plt.ylabel('elevator deflection [deg]')
plt.show()

de_da = (max(de_deg)-min(de_deg))/(max(alpha_deg_2)-min(alpha_deg_2))        #0.46 [-]

print('de/dalpha =')
print(de_da)

#Cm_delta, Cm_alpha, 20200310
#shift in center of gravity 
hp_ft = np.array([18360,18550])               #ft
V_IAS_kts = np.array([156,156])
de_deg = [-0.2,-0.8]
fuel_used_lbs = np.array([940,989])
xnose_inch = np.array([288,131])

hp_mt = hp_ft*0.3048
V_IAS_ms  = V_IAS_kts*0.5144
fuel_used_kg = fuel_used_lbs*0.45359237
xnose_m = xnose_inch* 0.0254

atmospheres = Atmosphere(hp_mt)
rho = atmospheres.density
V_TAS_ms = V_IAS_ms * np.sqrt(rho0/rho)
mass = mramp_kg - fuel_used_kg
mom = mass *xnose_m   #kgm
x_cg = mom/mass     #m
x_cg_1 = x_cg[0]
x_cg_2 = x_cg[1]
print(x_cg_1,x_cg_2) 

change_de = de_deg[1]-de_deg[0]
change_xcg = x_cg_2-x_cg_1
C_N = (mass[0]*9.81)/(0.5*rho[0]*V_TAS_ms[0]**2*S) #for steady horizontal flight
c_bar = c #MAC

Cm_delta = -1/change_de * C_N * change_xcg/c_bar
# Cm_delta = -1.1642 
Cm_alpha = Cm_delta*de_da
print(Cm_delta,Cm_alpha) 