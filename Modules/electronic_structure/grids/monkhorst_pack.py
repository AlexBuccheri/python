import numpy as np

# Sequence of grid points per dimension, defined by equation 3) in:
# https://journals.aps.org/prb/pdf/10.1103/PhysRevB.13.5188
def mk_function(ni):
    i = np.arange(1, ni + 1, dtype=int)
    return (2 * i - ni - 1) / (2 * ni)

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
                kgrid[:, k_cnt] = np.matmul(b, [i,j,k])
                k_cnt += 1

    return kgrid