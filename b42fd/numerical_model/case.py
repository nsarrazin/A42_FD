import numpy as np
import numba as nb
from b42fd.consts import *
from b42fd.numerical_model.thrust.compute_thrust import compute_thrust
from ambiance import Atmosphere   #need to install this package

class Case:
    def __init__(self, hp, Vc, TAT, MFl, MFr):
        self.hp = hp
        self.Vc = Vc
        self.TAT = TAT
        self.MFl = MFl
        self.MFr = MFr

    @property
    def p(self):
        return p0*(1+lamb*self.hp/T0)**(-g0/lamb/R)
    
    @property
    def M(self):
        C=1+(gamma-1)/2/gamma*rho0/p0*self.Vc**2
        D=gamma/(gamma-1)
        M=np.sqrt(2/(gamma-1)*((1+p0/self.p*(C**D-1))**(1/C)-1))
        return M

    @property
    def T(self):
        return self.TAT/(1+(gamma-1)*self.M**2)

    @property
    def T_isa(self):
        return Atmosphere(self.hp).temperature

    @property
    def deltaT(self):
        return self.T - self.T_isa

    @property
    def thrust_input(self):
        return np.array([self.hp, self.M, self.deltaT, self.MFl, self.MFr], dtype=np.float64)

    @property
    def thrust(self):
        inputs = self.thrust_input
        thrust = compute_thrust(inputs)
        return thrust

if __name__ == "__main__":
    h=np.array([18060, 18360, 18940, 18350, 18090, 17680, 18360])*ft
    V=np.array([156, 147, 134, 168, 176, 186, 156])*kt
    TAT=np.array([-10.2, -11.5, -13.5, -10.5, -9.5, -7.8, -11.2])+273.15
    MFl=np.array([409,407,404,410,416,420, 480])*lbshr
    MFr=np.array([470, 466, 455, 471, 448, 484, 469])*lbshr

    thrusts = []
    for i in range(len(h)):
        case = Case(h[i], V[i], TAT[i], MFl[i], MFr[i])
        thrusts.append(case.thrust)

    print(thrusts)