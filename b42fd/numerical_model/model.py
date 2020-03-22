import control
import numpy as np
# from b42fd.numerical_model.Cit_par import *
import matplotlib.pyplot as plt
from math import pi, sin, cos
import scipy.signal as signal


# Stationary flight condition

hp0    =   2700          # pressure altitude in the stationary flight condition [m]
V0     =   100          # true airspeed in the stationary flight condition [m/sec]
alpha0 =   .05          # angle of attack in the stationary flight condition [rad]
th0    =   .05          # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =  5000           # mass [kg]

# aerodynamic properties
e      = 0.8         # Oswald factor [ ]
CD0    = 0.04        # Zero lift drag coefficient [ ]
CLa    = 5.084       # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    = -0.5626     # longitudinal stabilty [ ]
Cmde   = -1.1642     # elevator effectiveness [ ]

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
rho    = rho0 * pow(((1+(lambd * hp0 / Temp0))), (-((g / (lambd*R)) + 1)))
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

# symm

# dimensionLess

C1_s_l = np.zeros((4, 4))
C1_s_l[0, 0] = -2 * muc * (c / (V0))
C1_s_l[1, 1] = (CZadot - 2 * muc) * (c / V0)
C1_s_l[2, 2] = -c / V0
C1_s_l[3, 3] = -2 * muc * KY2 * (c / V0)
C1_s_l[3, 1] = Cmadot * c / V0

C2_s_l = np.zeros((4, 4))
C2_s_l[0, 0] = CXu
C2_s_l[0, 1] = CXa
C2_s_l[0, 2] = CZ0
C2_s_l[0, 3] = CXq    #comment this line out to match the FD reader version, but we should include it in our model
C2_s_l[1, 0] = CZu
C2_s_l[1, 1] = CZa
C2_s_l[1, 2] = -CX0
C2_s_l[1, 3] = (CZq + 2 * muc)
C2_s_l[2, 3] = 1
C2_s_l[3, 0] = Cmu
C2_s_l[3, 1] = Cma
C2_s_l[3, 3] = Cmq

C3_s_l = np.array([[CXde], [CZde], [0], [Cmde]])

A_s_l = -1 * np.linalg.inv(C1_s_l) @ C2_s_l
B_s_l = -1 * np.linalg.inv(C1_s_l) @ C3_s_l
# print(A_s)
# print(np.linalg.eig(A_s)[0])
C_s_l = np.eye(4)
D_s_l = np.zeros((4, 1))

# if __name__== '__main__':
#     T = np.linspace(0, 10, 1000)
#     sys_s = control.ss(A_s_l, B_s_l, C_s_l, D_s_l)
#     T, yout = control.impulse_response(sys_s, T)
#
#     print(f'{np.linalg.eig(A_s_l)[0]}')
#     # print(A_s)
#
#     plt.plot(T, yout[1])
#     plt.xlim(0,10)
    # plt.show()


# dimensionHaving

C1_s_h = np.zeros((4, 4))
C1_s_h[0, 0] = -2 * muc * (c / (V0**2))
C1_s_h[1, 1] = (CZadot - 2 * muc) * (c / V0)
C1_s_h[2, 2] = -c / V0
C1_s_h[3, 3] = -2 * muc * KY2 * (c / V0)**2
C1_s_h[3, 1] = Cmadot * c / V0

C2_s_h = np.zeros((4, 4))
C2_s_h[0, 0] = CXu / V0
C2_s_h[0, 1] = CXa
C2_s_h[0, 2] = CZ0
C2_s_h[0, 3] = CXq * (c/V0)    #comment this line out to match the FD reader version, but we should include it in our model
C2_s_h[1, 0] = CZu / V0
C2_s_h[1, 1] = CZa
C2_s_h[1, 2] = -CX0
C2_s_h[1, 3] = (CZq + 2 * muc) * (c/V0)
C2_s_h[2, 3] = 1 * (c/V0)
C2_s_h[3, 0] = Cmu / V0
C2_s_h[3, 1] = Cma
C2_s_h[3, 3] = Cmq * (c/V0)

C3_s_h = np.array([[CXde], [CZde], [0], [Cmde]])

A_s_h = -1 * np.linalg.inv(C1_s_l) @ C2_s_l
B_s_h = -1 * np.linalg.inv(C1_s_l) @ C3_s_l
# print(A_s)
# print(np.linalg.eig(A_s)[0])
C_s_h = np.eye(4)
D_s_h = np.zeros((4, 1))





####################################################


# asymm


# dimensionLess

C1_a_l = np.zeros((4, 4))
C1_a_l[0, 0] = (CYbdot - 2 * mub) * (b / V0)
C1_a_l[1, 1] = -.5 * b / V0
C1_a_l[2, 2] = -4 * mub * KX2 * (b / V0)
C1_a_l[3, 3] = -4 * mub * KZ2 * (b / V0)
C1_a_l[2, 3] = 4 * mub * KXZ * (b / V0)
C1_a_l[3, 2] = 4 * mub * KXZ * (b / V0)
C1_a_l[3, 0] = Cnbdot * (b / V0)

C2_a_l = np.zeros((4, 4))
C2_a_l[0, 0] = CYb
C2_a_l[0, 1] = CL
C2_a_l[0, 2] = CYp
C2_a_l[0, 3] = (CYr - 4 * mub)
C2_a_l[1, 2] = 1
C2_a_l[2, 0] = Clb
C2_a_l[2, 2] = Clp
C2_a_l[2, 3] = Clr
C2_a_l[3, 0] = Cnb
C2_a_l[3, 2] = Cnp
C2_a_l[3, 3] = Cnr

C3_a_l = np.array([[CYda, CYdr], [0, 0], [Clda, Cldr], [Cnda, Cndr]])

A_a_l = -1 * np.linalg.inv(C1_a_l) @ C2_a_l
B_a_l = -1 * np.linalg.inv(C1_a_l) @ C3_a_l
C_a_l = np.eye(4)
D_a_l = np.zeros((4, 2))

# if __name__ == '__main__':
#     sys_a = control.ss(A_a, B_a, C_a, D_a)
#
#     T = np.linspace(0, 10, 1000)
#     T, yout = control.impulse_response(sys_a, T)
#
#     print(f'{np.linalg.eig(A_a)[0]}')
#     # print(A_a)
#
#     plt.plot(T, yout[1])
#     plt.xlim(0, 10)
#     plt.show()


# dimensionHaving

C1_a_h = np.zeros((4, 4))
C1_a_h[0, 0] = (CYbdot - 2 * mub) * (b / V0)
C1_a_h[1, 1] = -.5 * b / V0
C1_a_h[2, 2] = -2 * mub * KX2 * (b / V0)**2
C1_a_h[3, 3] = -2 * mub * KZ2 * (b / V0)**2
C1_a_h[2, 3] = 2 * mub * KXZ * (b / V0)**2
C1_a_h[3, 2] = 2 * mub * KXZ * (b / V0)**2
C1_a_h[3, 0] = Cnbdot * (b / V0)

C2_a_h = np.zeros((4, 4))
C2_a_h[0, 0] = CYb
C2_a_h[0, 1] = CL
C2_a_h[0, 2] = CYp * (b / (2*V0))
C2_a_h[0, 3] = (CYr - 4 * mub) * (b / (2*V0))
C2_a_h[1, 2] = 1 * (b / (2*V0))
C2_a_h[2, 0] = Clb
C2_a_h[2, 2] = Clp * (b / (2*V0))
C2_a_h[2, 3] = Clr * (b / (2*V0))
C2_a_h[3, 0] = Cnb
C2_a_h[3, 2] = Cnp * (b / (2*V0))
C2_a_h[3, 3] = Cnr * (b / (2*V0))

C3_a_h = np.array([[CYda, CYdr], [0, 0], [Clda, Cldr], [Cnda, Cndr]])

A_a_h = -1 * np.linalg.inv(C1_a_h) @ C2_a_h
B_a_h = -1 * np.linalg.inv(C1_a_h) @ C3_a_h
C_a_h = np.eye(4)
D_a_h = np.zeros((4, 2))




###### These are form the FD reader, difference being they assume CXq = 0


# T = np.zeros((4,4))
#
# T[0,0] = (V0/c)*(CXu/(2*muc))
# T[0,1] = (V0/c)*(CXa/(2*muc))
# T[0,2] = (V0/c)*(CZ0/(2*muc))
# T[0,3] = (V0/c)*(CXq/(2*muc))
#
# T[1,0] = (V0/c)*(CZu/(2*muc - CZadot))
# T[1,1] = (V0/c)*(CZa/(2*muc - CZadot))
# T[1,2] = -1*(V0/c)*(CX0/(2*muc-CZadot))
# T[1,3] = (V0/c)*(2*muc+CZq)/(2*muc - CZadot)
#
# T[2,3] = (V0/c)
#
# T[3,0] = (V0/c)*(Cmu + CZu*(Cmadot/(2*muc - CZadot)))/(2*muc*KY2)
# T[3,1] = (V0/c)*(Cma + CZa*(Cmadot/(2*muc - CZadot)))/(2*muc*KY2)
# T[3,2] = -1*(V0/c)*(CX0*(Cmadot/(2*muc - CZadot)))/(2*muc*KY2)
# T[3,3] = (V0/c)*(Cmq + Cmadot*(2*muc+CZq)/(2*muc-CZadot))/(2*muc*KY2)
#
#
#
# P = np.zeros((4,4))
#
# P[0,0] = (V0/b)*(CYb/(2*mub))
# P[0,1] = (V0/b)*(CL/(2*mub))
# P[0,2] = (V0/b)*(CYp/(2*mub))
# P[0,3] = (V0/b)*(CYr - 4*mub)/(2*mub)
#
# P[1,2] = 2*V0/b
#
# P[2,0] = (V0/b)*(Clb*KZ2 + Cnb*KXZ)/(4*mub*(KX2*KZ2 - KXZ**2))
# P[2,2] = (V0/b)*(Clp*KZ2 + Cnp*KXZ)/(4*mub*(KX2*KZ2 - KXZ**2))
# P[2,3] = (V0/b)*(Clr*KZ2 + Cnr*KXZ)/(4*mub*(KX2*KZ2 - KXZ**2))
#
# P[3,0] = (V0/b)*(Clb*KXZ + Cnb*KX2)/(4*mub*(KX2*KZ2 - KXZ**2))
# P[3,2] = (V0/b)*(Clp*KXZ + Cnp*KX2)/(4*mub*(KX2*KZ2 - KXZ**2))
# P[3,3] = (V0/b)*(Clr*KXZ + Cnr*KX2)/(4*mub*(KX2*KZ2 - KXZ**2))
#
# ########
#
# if __name__ == '__main__':
#
#     print(np.linalg.eig(A_s)[0])
#
#     print(np.linalg.eig(T)[0])
#
#     print(T)
#     print('')
#     print(A_s)
#     print('')
#     print(np.where(np.abs(T-A_s)<.001, True, False))
#
#     print(P)
#     print('')
#     print(A_a)
#     print('')
#     print(np.where(np.abs(P-A_a)<.001, True, False))

