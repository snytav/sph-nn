import torch
from scipy.special import lpmn
import numpy as np


def sum_spherical_coefs(m,n,fi,lb,koef):
    sum = 0.0
    for n in range(koef.shape[0]):
        sum += (koef[n][0]*np.cos(y)+koef[n][1]*np.sin(y))
    return sum

def get_spherical_harm(m, n, r, fi, lb, mass, koef):
    x = lb / 180.0 * torch.pi * torch.ones(1)
    z = torch.sin(x)

    y = fi / 360.0 * np.pi

    pmn = lpmn(m, n, z)
    p = pmn[m][n]
    s = sum_spherical_coefs(m,n,fi,lb,koef)
    sph = s**p*mass/r
    return sph