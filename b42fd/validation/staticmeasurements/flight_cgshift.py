import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import M_ramp, c
from b42fd.data_processing.thrust_input import pressure, Mach, corrected_temp,sound_speed, true_airspeed
from b42fd.consts import gamma,T0,lamb,g0,R,p0, rho0
from b42fd.numerical_model.case import Case
from b42fd.validation.staticmeasurements.flight_trimcurve import de_da
from b42fd.validation.staticmeasurements.flight_data import h_m, V_ms, TAT_K, alpha_deg, fuelburnt_kg, mramp_kg

print("Flight test data cg shift")

#fixed value
Ws = 60500 #[N] table B1 page 11 of assignment
S = 30.

#from 20200310_V2
alpha_deg = alpha_deg[13:15]
de_deg = [-0.2,-0.8]
de_rad = np.radians(de_deg)
h_m = h_m[13:15]
V_ms  = V_ms[13:15]
alpha_rad = np.radians(alpha_deg)
fuelburnt_kg = fuelburnt_kg[13:15]
TAT_K = TAT_K[13:15]

#density from ambiance package (not a standard package so install)
atmospheres = Atmosphere(h_m)
rho = atmospheres.density

#To find true airspeed
V_TAS=np.zeros(len(h_m))

for i in range(len(h_m)):
    Vc_ms=V_ms[i]
    hp_m=h_m[i]
    Tm_K=TAT_K[i]
    p=pressure(hp_m, gamma,T0,lamb,g0,R,p0)
    M=Mach(Vc_ms, gamma, rho0,p0, p)
    T=corrected_temp(Tm_K,M,gamma)
    a=sound_speed(gamma,R,T)
    V_TAS[i]=true_airspeed(M,a)


xnose_inch = np.array([288,131])

#conversion
xnose_m = xnose_inch* 0.0254


#calculations
mass = mramp_kg - fuelburnt_kg
W = mass*g0

x_cg_fuel = 287.58                  #inches
x_cg_fuel_m = x_cg_fuel*0.0254
m_shift = 96                        #kg

def cg_shift(x_cg_fuel_m, xnose_m, mass_3, m_shift):
    mom = []
    mom = np.append(mom, mass_3[0]*x_cg_fuel_m)
    mom = np.append(mom, -m_shift*xnose_m[0])
    mom_tot = np.sum(mom)
    x_cg = []
    # print(mom_tot/(mass_3[0]-m_shift))
    x_cg = np.append(x_cg, mom_tot/(mass_3[0]-m_shift))

    mom = np.append(mom, m_shift*xnose_m[1])
    mom_tot = np.sum(mom)
    x_cg = np.append(x_cg, mom_tot/mass_3[1])
    return x_cg[0],x_cg[1]


V_EAS = V_TAS*np.sqrt(rho/rho0)

# print(V_EAS_ms)
Vej = V_EAS[0]*np.sqrt(W[1]/W[0])

change_de = de_rad[1]-de_rad[0]
x_cg_1,x_cg_2 = cg_shift(x_cg_fuel_m, xnose_m, mass, m_shift)
change_xcg = x_cg_2-x_cg_1

print("change_xcg")
print(change_xcg)
C_N = W[1]/(0.5*rho0*Vej**2*S) #for steady horizontal flight
print("C_N:", C_N)
c_bar = c #MAC

Cm_delta = -C_N/change_de * change_xcg/c_bar    # -1.1642 
Cm_alpha = -Cm_delta*de_da
print("Cm_delta:")
print(Cm_delta)

print("Cm_alpha:")
print(Cm_alpha) 

