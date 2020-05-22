# Check lattice parameters are consistent with entos lattice vectors

import numpy as np

from modules.electronic_structure.structure import bravais
from modules.electronic_structure.structure import lattice as lt

# cell volume and cell angles -> Should definitely add this to the entos check

# Anatase


# Rutile
lattice = bravais.simple_tetragonal(a=4.60677734,c=2.99175662)

V = lt.parallelpiped_volume(lattice)
assert np.isclose(V, 63.49224796, atol=1.e-6)

(alpha, beta, gamma) = lt.cell_angles(lattice, return_unit='degrees')
assert np.allclose((alpha, beta, gamma), 90.0)


# Graphite. Ref: http://aflowlib.org/CrystalDatabase/A_hP4_194_bc.html
lattice = bravais.hexagonal(a=2.464, c=6.711)

V = lt.parallelpiped_volume(lattice)
assert np.isclose(V, 35.285743, atol=1.e-6)

(alpha, beta, gamma) = lt.cell_angles(lattice, return_unit='degrees')
assert np.allclose((alpha, beta), 90.0)
assert np.allclose(gamma, 120.0)


# Boron nitride - hexagonal
# Ref: https://materialsproject.org/materials/mp-984/
lattice = bravais.hexagonal(a=2.51242804, c=7.70726501)

V = lt.parallelpiped_volume(lattice)
assert np.isclose(V, 42.132593, atol=1.e-6)

(alpha, beta, gamma) = lt.cell_angles(lattice, return_unit='degrees')
assert np.allclose((alpha, beta), 90.0)
assert np.allclose(gamma, 120.0)