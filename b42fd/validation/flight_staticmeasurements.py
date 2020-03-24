import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from b42fd.numerical_model.Cit_par import M_ramp, c
from b42fd.data_processing.thrust_input import pressure, Mach, corrected_temp,sound_speed, true_airspeed
from b42fd.consts import gamma,T0,lamb,g0,R,p0, rho0
from b42fd.numerical_model.case import Case
print("Flight test data")

#from 20200310_V2
h = np.array([18000,17990,18020,17990,18000,18010,18060,18360,18940,18350,18090,17680,18360,18369,18550]) #altitude 
V = np.array([161,131,222,200,182,114,156,147,134,168,176,186,156,156,156])#indicated airspeed
TAT_C = np.array([-9.5,-11.5,-4.8,-6.8,-8.2,-12.8,-10.2,-11.5,-13.5,-10.5,-9.5,-7.8,-11.2,-11.2,-11.2]) #temp
alpha_deg = np.array([5,8.1,2.0,2.9,3.6,10.7,5.2,6.3,7.5,4.4,3.8,3.3,5.2,5.2,5.2])  #aoa deg
fuelburnt_lbs = np.array([583,569,634,665,688,729,811,840,865,888,901,912,940,940,989]) #fuel burnt lbs
# print(len(h),len(V),len(TAT),len(alpha_deg),len(fuelburnt_lbs)) #, first 6 = CL-CD series, middle 7 = elevator trim curve, last 2 = shift in center of gravity
mramp_lbs = M_ramp

#conversions 
h_m = h*0.3048
V_ms  = V*0.5144
alpha_rad = np.radians(alpha_deg)
mramp_kg = mramp_lbs*0.45359237
fuelburnt_kg = fuelburnt_lbs*0.45359237
TAT_K = TAT_C+273.15

# print(TAT_K)

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

plt.plot(alpha_deg_1,CL,'x')
plt.xlabel('angle of attack [deg]')
plt.ylabel('lift coefficient [-]')
plt.show()

#print CL_alpha
print('Cl_a =') 
print(CL_alpha_rad)     # 0.3522578695587925   [1/radian]



trim curve
2nd set of values


#take values from before
indices = [6,7,8,9,10,11,12]
alpha_deg_2 = np.take(alpha_deg,indices)
alpha_rad_2=np.radians(alpha_deg_2)
#take value from excel
de_deg_2 = np.array([-0.3,-0.7,-1.2,0.1,0.4,0.7,-0.2]) 

de_rad_2=np.radians(de_deg_2)
#plot trim curve
plt.plot(alpha_deg_2,de_deg_2,'x')
plt.xlabel('angle of attack [deg]')
plt.ylabel('elevator deflection [deg]')
plt.show()

de_da = -(max(de_deg_2)-min(de_deg_2))/(max(alpha_deg_2)-min(alpha_deg_2))        #0.46 [-]

print('de/dalpha =')
print(de_da)


Cm_delta, Cm_alpha
shift in center of gravity = 3rd set of values


indices = [13,14]
hp_m_3 = np.take(h_m,indices)
V_TAS_ms_3 = np.take(V_TAS,indices)
# print(V_TAS_ms_3)
de_deg_3 = [-0.2,-0.8]
de_rad_3 = np.radians(de_deg_3)

xnose_inch = np.array([288,131])
rho_3 = np.take(rho,indices)
fuelburnt_kg_3 = np.take(fuelburnt_kg,indices)
#conversion
xnose_m = xnose_inch* 0.0254

mass_3 = mramp_kg - fuelburnt_kg_3
W = mass_3*g0
# print(W)

x_cg_fuel = 287.58 #inches
x_cg_fuel_m = x_cg_fuel*0.0254
m_shift = 96 #kg

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

change_de = de_rad_3[1]-de_rad_3[0]
x_cg_1,x_cg_2 = cg_shift(x_cg_fuel_m, xnose_m, mass_3, m_shift)
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


MFl=np.array([392,369,608,508,453,431])/7936.64
MFr=np.array([450, 378,668, 548, 488, 480])/7936.64


print(V_TAS_ms_1)

thrusts = []
for i in range(6):
    case = Case(h_m[i], V_ms[i], TAT_K[i], MFl[i], MFr[i])
    thrusts.append(case.thrust)

CD = thrusts/(0.5*rho[0:6]*V_TAS_ms_1[0:6]**2*S)

plt.plot(alpha_deg[:6],CD,'x')
plt.xlabel('alpha [deg]')
plt.ylabel('CD [-]')
plt.show()

plt.plot([0,0],[0.008,0.02])
plt.plot(CL,CD,'x')
plt.xlabel('CL [-]')
plt.ylabel('CD [-]')
plt.legend()
plt.show()


CD0 = 0.017 #graphically determined

A = 15.911**2/30.
e = 0.9

CDtheoretical = CD0 + CL**2/(np.pi*A*e)

plt.plot(CL**2, CDtheoretical, 'o-')
plt.plot(CL**2, CD,'x')
plt.xlabel('CL^2 [-]')
plt.ylabel('CD')
plt.show()

#e = (CD - CD0)/(CL**2)*(np.pi*A)
#print(e)
"""

Elevator control force curve
flight data
4th set of data
"""

Fe_N = np.array([1,-14, -31, 26, 50, 83, 1])
indices = [6,7,8,9,10,11,12]
V_TAS_ms_4 = np.take(V_TAS,indices)
rho_4 = np.take(rho,indices)
fuelburnt_kg_4 = np.take(fuelburnt_kg,indices)

mass_4 = mramp_kg - fuelburnt_kg_4
V_EAS_ms_4 = V_TAS_ms_4 * np.sqrt(rho_4/rho0)

W_s = 60500 #N
W_4 = mass_4*g0 #N

V_hat_e_ms = V_EAS_ms_4 * np.sqrt(W_s/W_4)

s = sorted(zip(V_hat_e_ms,Fe_N))
V_hat_e_ms,Fe_N = map(list, zip(*s))
# print(V_hat_e_ms,Fe_N)
plt.plot(V_hat_e_ms,Fe_N,'-')
plt.ylim(90,-40)
plt.hlines(0,np.min(V_hat_e_ms),np.max(V_hat_e_ms),linestyles="dashed")
plt.xlabel("Reduced equivalent airspeed [m/s]")
plt.ylabel("Measured elevator control force [N]")
plt.show()

