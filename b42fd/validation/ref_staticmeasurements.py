import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import c
from b42fd.data_processing.thrust_input import pressure, Mach, corrected_temp,sound_speed, true_airspeed
from b42fd.consts import gamma,T0,lamb,g0,R,p0, rho0
from b42fd.numerical_model.case import Case

print("Ref test data")
"""
mass calculation for ref data
"""
m = np.array([95,92,74,66,61,75,78,86,68]) #[kg] 
m_pax = m*2.2046  #[lbs]
# print(m_pax)

M_payload = np.sum(m_pax)
# print(M_payload)
BEW = 9165 # basic empty weight [lbs]
ZFW = BEW + M_payload
fuel = 4050 # [lbs]
M_ramp = fuel + ZFW
# print(M_ramp)

#from ref data
h = np.array([5010,5020,5020,5030,5020,5110,6060,6350,6550,6880,6160,5810,5310,5730,5790]) #altitude 
V = np.array([249,221,192,163,130,118,161,150,140,130,173,179,192,161,161])#indicated airspeed
TAT_C = np.array([12.5,10.5,8.8,7.2,6,5.2,5.5,4.5,3.5,2.5,5.0,6.2,8.2,5,5]) #temp
alpha_deg = np.array([1.7,2.4,3.6,5.4,8.7,10.6,5.3,6.3,7.3,8.5,4.5,4.1,3.4,5.3,5.3])  #aoa deg
fuelburnt_lbs = np.array([360,412,447,478,532,570,664,694,730,755,798,825,846,881,910]) #fuel burnt lbs
# print(len(h),len(V),len(TAT),len(alpha_deg),len(fuelburnt_lbs)) #, first 6 = CL-CD series, middle 7 = elevator trim curve, last 2 = shift in center of gravity
mramp_lbs = M_ramp  #remember to change for reference data

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

# print(V_TAS)
"""
CL-CD and CL-a curves 
1st set of values
"""

#take values from before
indices = [0,1,2,3,4,5]
alpha_deg_1 = np.take(alpha_deg,indices)
alpha_rad_1 = np.take(alpha_rad,indices)
V_TAS_ms_1 = np.take(V_TAS,indices)          
hp_m_1 = np.take(h_m,indices)
fuelburnt_kg_1 = np.take(fuelburnt_kg,indices)             
rho_1 = np.take(rho,indices)

#fixed values
S = 30                  #[m^2]

#calculations
mass_1 = mramp_kg - fuelburnt_kg_1
CL = (mass_1*g0)/(0.5*rho_1*V_TAS_ms_1**2*S)
# print(CL)
# print(alpha_deg_1)

CL_alpha_deg = (CL[-1] - CL[0])/(alpha_deg_1[-1]-alpha_deg_1[0])        # in 1/degree
CL_alpha_rad = (CL[-1] - CL[0])/(alpha_rad_1[-1]-alpha_rad_1[0])        # in 1/radian

#CL-alpha plot

plt.plot(alpha_deg_1,CL,'x-')
plt.xlabel('angle of attack [deg]')
plt.ylabel('lift coefficient [-]')
plt.show()

#print CL_alpha
print('Cl_a =') 
print(CL_alpha_rad)     # 0.3522578695587925   [1/radian]


"""
trim curve
2nd set of values
"""

#take values from before
indices = [6,7,8,9,10,11,12]
alpha_deg_2 = np.take(alpha_deg,indices)

#take value from excel
de_deg_2 = np.array([0,-0.4,-0.9,-1.5,0.4,0.6,1]) 

#plot trim curve
plt.plot(alpha_deg_2,de_deg_2,'x')
plt.xlabel('angle of attack [deg]')
plt.ylabel('elevator deflection [deg]')
plt.show()

de_da = -(max(de_deg_2)-min(de_deg_2))/(max(alpha_deg_2)-min(alpha_deg_2))        #0.46 [-]

print('de/dalpha =')
print(de_da)

"""
Cm_delta, Cm_alpha
shift in center of gravity = 3rd set of values
"""

indices = [13,14]
hp_m_3 = np.take(h_m,indices)
V_TAS_ms_3 = np.take(V_TAS,indices)
# print(V_TAS_ms_3)
de_deg_3 = [0,-0.5]
fuelburnt_kg_3 = np.take(fuelburnt_kg,indices)
xnose_inch = np.array([288,131])
rho_3 = np.take(rho,indices)

#conversion
xnose_m = xnose_inch* 0.0254

mass_3 = mramp_kg - fuelburnt_kg_3
W = mass_3*g0
# print(W)

x_cg_fuel = 287.58 #inches
x_cg_fuel_m = x_cg_fuel*0.0254
m_shift = 68 #kg

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

V_EAS_ms = V_TAS_ms_3*np.sqrt(rho_3/rho0)
# print(V_EAS_ms)
Vej = V_EAS_ms[0]*np.sqrt(W[1]/W[0])

change_de = de_deg_3[1]-de_deg_3[0]
x_cg_1,x_cg_2 = cg_shift(x_cg_fuel_m, xnose_m, mass_3, m_shift)
change_xcg = x_cg_2-x_cg_1
print("cg change:", change_xcg)
C_N = (W[1])/(0.5*rho0*Vej**2*S) #for steady horizontal flight
print("C_N:", C_N)
c_bar = c #MAC

Cm_delta = -1/change_de * C_N * change_xcg/c_bar    # -1.1642 
Cm_alpha = -Cm_delta*de_da
print("Cm_delta:")
print(Cm_delta)
print("Cm_alpha:")
print(Cm_alpha) 

MFl=np.array([798, 673, 561, 463, 443, 474])/7936.64
MFr=np.array([813, 682, 579, 484, 467, 499])/7936.64

thrusts = []
for i in range(6):
    case = Case(h_m[i], V_ms[i], TAT_K[i], MFl[i], MFr[i])
    thrusts.append(case.thrust)

CD = thrusts/(0.5*rho[0:6]*V_TAS_ms_1[0:6]**2*S)

plt.plot(CL,CD,'x')
plt.xlabel('CL [-]')
plt.ylabel('CD [-]')
plt.show()