import json
import control
import numpy as np
import matplotlib.pyplot as plt
from b42fd.numerical_model.model import A_s_h, B_s_h, C_s_h, D_s_h, A_a_h, B_a_h, C_a_h, D_a_h, A_s_l, B_s_l, C_s_l, D_s_l

with open('input.json', 'r') as f:
    raw = f.read()
    data = json.loads(raw)

print(data.keys())
print(data["phugoid"].keys())

ts = data['phugoid']['t']
ta = data['dutchroll']['t']
dele = data['phugoid']['delta_e']

dela = data['dutchroll']['delta_a']
delr = data['dutchroll']['delta_r']

sys_s = control.ss(A_s_h, B_s_h, C_s_h, D_s_h)
sys_a = control.ss(A_s_l, B_s_l, C_s_l, D_s_l)

# print(list(zip(dela, delr)))


T, yout, bla = control.forced_response(sys_s, ts, dele)
for i in range(4):
    plt.plot(T, yout[i])
    # plt.plot(t, dele)
    plt.show()

# T, yout, bla = control.forced_response(sys_a, ts, dele)
# for i in range(4):
#     plt.plot(T, yout[i])
#     # plt.plot(t, dele)
#     plt.show()

# T, yout, bla = control.forced_response(sys_a, ta, np.array(list(zip(dela, delr))).T)
# for i in range(4):
#     plt.plot(T, yout[i])
#     # plt.plot(t, dele)
#     plt.show()


