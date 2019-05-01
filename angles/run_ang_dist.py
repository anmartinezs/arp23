"""
Script to simulate angular distribution between two unit vectors.

# Author: Antonio Martinez-Sanchez (Max Planck Institute for Biochemistry)
# Date: 24.01.17
"""

__author__ = 'Antonio Martinez-Sanchez'

from angles import *
from matplotlib import pyplot as plt

####### Input parameters

M = 2000 # Number of simulation
N = 1000 # Pair of vectors per simulation
B = 50 # Number of bins for histogram simulation
CI = 99 # Confidence interval (%)

####### Analytic solution

phi = np.linspace(0, np.pi, B)
P_phi = .5 * np.sin(phi)


###### Simulation

sims, bins = gen_mat_sim(M, N, B)
P_l = func_envelope(sims, per=CI)
P_m = func_envelope(sims, per=50)
P_h = func_envelope(sims, per=100-CI)

###### Plotting

plt.title('Angle distribution random unit vector pairs')
plt.xlabel('phi (rad)')
plt.xlabel('P_phi')
plt.plot(phi, P_phi, 'b')
plt.plot(bins, P_l, 'k--')
plt.plot(bins, P_m, 'k')
plt.plot(bins, P_h, 'k--')
plt.show(block=True)
