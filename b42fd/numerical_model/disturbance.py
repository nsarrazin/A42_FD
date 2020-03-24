from b42fd.numerical_model.model import A_s_h, B_s_h, C_s_h, D_s_h, A_a_h, B_a_h, C_a_h, D_a_h
import control
import numpy as np
import matplotlib.pyplot as plt




ss = control.StateSpace(A_s_h, B_s_h, C_s_h, D_s_h)



if __name__== '__main__':
    T = np.linspace(0, 100, 1000)
    T, yout = control.initial_response(ss, T, [0, .05, .05, 0])
    for i in range(4):
        plt.plot(T, yout[i])
        # plt.xlim(0,10)
        plt.show()