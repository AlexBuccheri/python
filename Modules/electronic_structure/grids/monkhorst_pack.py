import numpy as np
from typing import List

def even(i):
    if i % 2 == 0:
        return True
    else:
        return False


# Sequence of grid points per dimension, defined by equation 3) in:
# https://journals.aps.org/prb/pdf/10.1103/PhysRevB.13.5188
def mk_function(ni: int):
    assert ni > 0 
    i = np.arange(1, ni + 1, dtype=int)
    #shift = (0.5, 0.)[even(ni)]
    shift = 0
    return (2 * i - ni - 1) / (2 * ni) + shift


# Sequence of points per dimension, defined by equation 1) in:
# https://doi.org/10.1103/PhysRevB.16.1748
# for hexagonal grids

# This appears to be a) Not centred on gamma
# b) Always positive
# c) Always include 0
# https://www.c2x.org.uk/mp_kpoints.html

def mk_function2(ni: int):
    assert ni > 0
    i = np.arange(1, ni + 1, dtype=int)
    return (i - 1) / ni


# Monkhorst-Pack grid, as defined by equation 4) in:
# https://journals.aps.org/prb/pdf/10.1103/PhysRevB.13.5188
# Even grids avoid Gamma point, odd grids include Gamma point
#
# n = Number of sampling points per dimension (uniform grid)
# b = primitive reciprocal lattice vectors, expected column-wise
def grid(n, b):

    ikappa = mk_function(n[0])
    jkappa = mk_function(n[1])
    kkappa = mk_function(n[2])
    kgrid = np.empty(shape=(3, np.prod(n)))

    k_cnt = 0
    for i in ikappa:
        for j in jkappa:
            for k in kkappa:
                kgrid[:, k_cnt] = np.matmul(b, [i, j, k])
                k_cnt += 1

    return kgrid

# Sequence of points per dimension, defined by equation 1) in:
# https://doi.org/10.1103/PhysRevB.16.1748
def hexagonal_grid(n, b):
    assert len(n) == 3

    u_p = mk_function2(n[0])
    u_r = mk_function2(n[1])
    u_s = mk_function(n[2])
    kgrid = np.empty(shape=(3, np.prod(n)))

    k_cnt = 0
    for p in u_p:
        for r in u_r:
            for s in u_s:
                kgrid[:, k_cnt] = np.matmul(b, [p,r,s])
                k_cnt += 1

    return kgrid
