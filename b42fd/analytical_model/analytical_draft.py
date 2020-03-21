# Citation 550 - Linear simulation

from math import pi, sin, cos
from pathlib import Path
import numpy as np

path = "data/ref_data/ref_data.json"
# data = load_data(path)


# Aircraft mass, post_flight_datasheet 
# FFl=data["lh_engine_FMF"]["data"]
m = [95,92,74,66,61,75,78,86,68] #[kg] 
m_pax = np.array(m, dtype=int)*2.2046  #[lbs]
# print(m_pax)

M_payload = np.sum(m_pax)
# print(M_payload)
BEW = 9165 # basic empty weight [lbs]
ZFW = BEW + M_payload
fuel = 4050 # [lbs]
M_ramp = fuel + ZFW
# print(M_ramp)

# m_flow_l = FFl 
# m_flow_r = FFr
# m_flow = m_flow_l + m_flow_r

m_fuel_used = [360,412,447,478,532,570,664,694,730,755,798,825,846,881,910]

for m_fuel_used_i in m_fuel_used:
    W = M_ramp - m_fuel_used_i
    # print(W)

# cg calculation
x_datum = [131,131,170,214,214,251,251,288,288]     #fixed for ac
mom_tot = 0
mom_pay_tot = 0 

for i in range(len(x_datum)):
    mom = m_pax[i]*x_datum[i]
    mom_pay_tot += mom
    # print(mom_pay_tot) 

nose = 1080     #jackpads
main_r = 4430   #jackpads
main_l = 4405   #jackpads
x_cg_jack = 315.5 - (221.8*nose)/(nose+main_r+main_l)

mom_bew = BEW * x_cg_jack
mom_pay = mom_pay_tot
mom_zfw = mom_bew + mom_pay
x_cg_zfw = mom_zfw/ZFW
x_cg_fuel = 297.58 
mom_fuel = fuel*x_cg_fuel
mom_ramp = mom_fuel + mom_zfw
x_cg_ramp = mom_ramp/M_ramp
# print(x_cg_ramp)


# Stationary flight condition

hp0    =   1          # pressure altitude in the stationary flight condition [m]
V0     =   1          # true airspeed in the stationary flight condition [m/sec]
alpha0 =   1          # angle of attack in the stationary flight condition [rad]
th0    =   1          # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =   6000          # mass [kg]

# aerodynamic properties
e      = 0.8         # Oswald factor [ ]
CD0    = 0.04        # Zero lift drag coefficient [ ]
CLa    = 5.084       # Slope of CL-alpha curve [ ]

# Longitudinal stability (These parameters need to be improved)
Cma    = -0.5626     # longitudinal stabilty [ ]
Cmde   = -1.1642     # elevator effectiveness [ ]

# coefficients above need to be extracted from the flight data


# Aircraft geometry
S      = 30.00	         # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	         # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	         # mean aerodynamic cord [m]
lh_c   = lh / c	         # [ ]
b      = 15.911	         # wing span [m]
bh     = 5.791	         # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	         # [ ]
ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity
rho0   = 1.2250          # air density at sea level [kg/m^3] 
lambd = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]  
rho    = rho0 * pow( ((1+(lambd * hp0 / Temp0))), (-((g / (lambd*R)) + 1)))   
W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.095
CXa    = +0.47966		# Positive! (has been erroneously negative since 1993) 
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939

#Matrix form (vague)
    
import numpy as np

def cal_eigenvalues(A,B,C):
    return (complex(-B/2/A,np.sqrt(4*A*C-B**2)/2/A),complex(-B/2/A,-np.sqrt(4*A*C-B**2)/2/A) )

#Short period
"""Asp=-2*muc*KY2
Bsp=Cmadot+Cmq
Csp=Cma"""

Asp=2*muc*KY2*(2*muc-CZadot)
Bsp=-2*muc*KY2*CZa-(2*muc+CZq)*Cmadot-(2*muc-CZadot)*Cmq
Csp=CZa*Cmq-(2*muc+CZq)*Cma

print("eigenvalues for short period motion are:", cal_eigenvalues(Asp, Bsp, Csp))

#Phugoid oscillation
"""Aph=-4*muc**2
Bph=2*muc*CXu
Cph=-CZu*CZ0"""

Aph=2*muc*(CZa*Cmq-2*muc*Cma)
Bph=2*muc*(CXu*Cma-Cmu*CXa)+Cmq*(CZu*CXa-CXu*CZa)
Cph=CZ0*(Cmu*CZa-Cma*CZu)

print("eigenvalues for phugoid motion are:", cal_eigenvalues(Aph, Bph, Cph))

#Dutch roll
"""Adr=-2*mub*KZ2
Bdr=0.5*Cnr
Cdr=-Cnb"""

Aro=8*mub**2*KZ2
Bro=-2*mub*(Cnr+2*KZ2*CYb)
Cro=4*mub*Cnb+CYb*Cnr

print("eigenvalues for Dutch Roll are:", cal_eigenvalues(Aro, Bro, Cro))

#aperiodic rolling motion (still needs to be changed)

eig1=Clp/(4*mub*KX2)
print("eigenvalues for aperiodic rolling motion are:", eig1)
      
#aperiodic spiral motion (still needs to be changed)

eig2=(2*CL*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))
print("eigenvalues for aperiodic spiral motion are:", eig2)

