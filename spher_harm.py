import torch
from scipy.special import lpmn
import numpy as np


def sum_spherical_coefs(m,n,fi,al,koef):

    y = fi / 360.0 * np.pi
    z = np.sin(y)

    r0 = 0.0
    r1 = 0.0
    for n in range(m,koef.shape[0]):
        pmn = lpmn(m, n, z)
        p1 = pmn[m][0][n]
        r0 += koef[n][0]*p1
        r1 += koef[n][1] * p1
    return r0,r1

def get_spherical_harm(m, n, r, fi, al, mass, koef):
    x = al / 180.0 * np.pi
    z = np.sin(x)

    y = fi / 180.0 * np.pi
    z = np.sin(y)
    pmn = lpmn(m, n, z)
    p = pmn[m][0][n]
    r0,r1 = sum_spherical_coefs(m,n,fi,al,koef)
    sph = s**p*mass/r
    return sph