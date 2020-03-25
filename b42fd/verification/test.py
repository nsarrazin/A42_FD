
from b42fd.numerical_model.model import get_state_space
from b42fd.analytical_model.analytical_full import Analytical_Eigenmotion
from b42fd.validation.fuelmass import data
from b42fd.analytical_model.time import TimeTool
from b42fd.helpers import load_data
import numpy as np


t_phugoid    = 53*60
t_spm        = 58*60
t_dutchroll  = 60*60   #there is another dutch roll in yawning direction. not sure if we had to use it. Time=61*60
t_ape_roll   = 57*60
t_ape_spiral = 62*60


A_s_phugoid, A_a_phugoid = get_state_space(t_phugoid)[0], get_state_space(t_phugoid)[4]
A_s_spm, A_a_spm = get_state_space(t_spm)[0], get_state_space(t_spm)[4]
A_s_dutchroll, A_a_dutchroll = get_state_space(t_dutchroll)[0], get_state_space(t_dutchroll)[0]
A_s_ape_roll, A_a_ape_roll = get_state_space(t_ape_roll)[0], get_state_space(t_ape_roll)[0]
A_s_ape_spiral, A_a_ape_spiral = get_state_space(t_phugoid)[0], get_state_space(t_phugoid)[0]

eigs = np.linalg.eig
print(eigs(A_s_phugoid)[0], eigs(A_a_phugoid)[0])
print(eigs(A_s_spm)[0], eigs(A_a_spm)[0])
print(eigs(A_s_dutchroll)[0], eigs(A_a_dutchroll)[0])
print(eigs(A_s_ape_roll)[0], eigs(A_a_ape_roll)[0])
print(eigs(A_s_ape_spiral)[0], eigs(A_a_ape_spiral)[0])


m_pax = np.array([95, 102, 89, 82, 66, 81, 69, 85, 96])  # passenger weights in kg
M_e = 9165 * 0.453592  # empty aircraft weight in kg
M_u = 2640 * 0.453592  # mass of fuel

# stationary mesurements results
Cmde = -1.491241347862329
Cma = -0.6746091811758155

CLa = 4.371485054942859
CD0 = 0.016
e = 0.6

phugoid= Analytical_Eigenmotion(data, "phugoid oscillation", t=t_phugoid, M_u=M_u,CLa=CLa, CD0=CD0,  e=e,  Cma=Cma)
dutch_roll= Analytical_Eigenmotion(data, "dutch roll", t=t_phugoid, M_u=M_u, CLa=CLa, CD0=CD0,  e=e,  Cma=Cma)
aperiodic_roll= Analytical_Eigenmotion(data, "aperiodic roll", t=t_ape_roll, M_u=M_u, CLa=CLa, CD0=CD0,  e=e,  Cma=Cma)
aperiodic_spiral=Analytical_Eigenmotion(data, "aperiodic spiral", t=t_ape_spiral, M_u=M_u, CLa=CLa, CD0=CD0,  e=e,  Cma=Cma)

print(phugoid.eigvalues)
