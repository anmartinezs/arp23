"""
Utilities to study angular distribution between two vectors.

# Author: Antonio Martinez-Sanchez (Max Planck Institute for Biochemistry)
# Date: 24.01.17
"""

__author__ = 'Antonio Martinez-Sanchez'

import numpy as np


# Angle between two array of unit vectors
def ang_3d_unit_vect(v0, v1):
    return np.arccos(v0[:, 0]*v1[:, 0] + v0[:, 1]*v1[:, 1] + v0[:, 2]*v1[:, 2])


# Generates an array of random 3d uni vectors
# n: the number of vectors
def gen_rnd_unit_3d_vectors(n):

    # Random angles
    rho = 2.*np.pi * np.random.random(n)
    phi = np.arccos(2.*np.random.random(n) - 1.)

    # Array of vectors (spherical to cartesian)
    vects = np.zeros(shape=(int(n), 3), dtype=np.float)
    vects[:, 0] = np.cos(rho) * np.sin(phi)
    vects[:, 1] = np.sin(rho) * np.sin(phi)
    vects[:, 2] = np.cos(phi)

    return vects


# Generates a matrix with angles between arbitrary unit vectors
# m: number of simulation
# n: number of samples per simulation
# b: number of bins for the histogram
def gen_mat_sim(m, n, b):

    # Initialization
    sims = np.zeros(shape=(int(m), int(b)), dtype=np.float)

    # Simulation
    for i in range(m):
        hold = ang_3d_unit_vect(gen_rnd_unit_3d_vectors(n), gen_rnd_unit_3d_vectors(n))
        hist, bins = np.histogram(hold, bins=int(b), density=True)
        sims[i, :] = hist

    return sims, bins[:-1] + .5*bins[0]


# Computes the envelope a simulation matrix
# sims: simulation matrix
# per: percentile for the envelope, default is 50 (median)
def func_envelope(sims, per=50):
    return np.percentile(sims, per, axis=0)




