import control
import numpy as np
from numerical_model.Cit_par import *
import matplotlib.pyplot as plt
import scipy.signal as signal

# symm

C1_s = np.zeros((4, 4))
C1_s[0, 0] = -2 * muc * (c / (V0 ** 2))
C1_s[1, 1] = (CZadot - 2 * muc) * (c / V0)
C1_s[2, 2] = -c / V0
C1_s[3, 3] = -2 * muc * KY2 * (c / V0) ** 2
C1_s[3, 1] = Cmadot * c / V0

C2_s = np.zeros((4, 4))
C2_s[0, 0] = CXu / V0
C2_s[0, 1] = CXa
C2_s[0, 2] = CZ0
C2_s[0, 3] = CXq * c / V0
C2_s[1, 0] = CZu / V0
C2_s[1, 1] = CZa
C2_s[1, 2] = -CX0
C2_s[1, 3] = (CZq + 2 * muc) * c / V0
C2_s[2, 3] = c / V0
C2_s[3, 0] = Cmu / V0
C2_s[3, 1] = Cma
C2_s[3, 3] = Cmq * c / V0

C3_s = np.array([[CXde], [CZde], [0], [Cmde]])

A_s = -1 * np.linalg.inv(C1_s) @ C2_s
B_s = -1 * np.linalg.inv(C1_s) @ C3_s
# print(A_s)
# print(np.linalg.eig(A_s)[0])
C_s = np.eye(4)
D_s = np.zeros((4, 1))

# if __name__== '__main__':
# T = np.linspace(0, 10, 1000)
# sys_s = control.ss(A_s, B_s, C_s, D_s)
# T, yout = control.impulse_response(sys_s, T)
#
# print(np.linalg.eig(A_s)[0])
# print(A_s)
#
# plt.plot(T, yout[1])
# plt.xlim(0,10)
# plt.show()


C1_a = np.zeros((4, 4))
C1_a[0, 0] = (CYbdot - 2 * mub) * (b / V0)
C1_a[1, 1] = -.5 * b / V0
C1_a[2, 2] = -2 * mub * KX2 * (b / V0) ** 2
C1_a[3, 3] = -2 * mub * KZ2 * (b / V0) ** 2
C1_a[2, 3] = 2 * mub * KXZ * (b / V0) ** 2
C1_a[3, 2] = 2 * mub * KXZ * (b / V0) ** 2

C2_a = np.zeros((4, 4))
C2_a[0, 0] = CYb
C2_a[0, 1] = CL
C2_a[0, 2] = CYp * b / (2 * V0)
C2_a[0, 3] = (CYr - 4 * mub) * (b / (2 * V0))
C2_a[1, 2] = b / (2 * V0)
C2_a[2, 0] = Clb
C2_a[2, 2] = Clp * b / (2 * V0)
C2_a[2, 3] = Clr * b / (2 * V0)
C2_a[3, 0] = Cnb
C2_a[3, 2] = Cnp * b / (2 * V0)
C2_a[3, 3] = Cnr * b / (2 * V0)

C3_a = np.array([[CYda, CYdr], [0, 0], [Clda, Cldr], [Cnda, Cndr]])

A_a = -1 * np.linalg.inv(C1_a) @ C2_a
B_a = -1 * np.linalg.inv(C1_a) @ C3_a
C_a = np.eye(4)
D_a = np.zeros((4, 2))

if __name__ == '__main__':
    sys_a = control.ss(A_a, B_a, C_a, D_a)

    T = np.linspace(0, 10, 1000)
    T, yout = control.impulse_response(sys_a, T)

    print(np.linalg.eig(A_a)[0])
    print(A_a)

    plt.plot(T, yout[1])
    plt.xlim(0, 10)
    plt.show()