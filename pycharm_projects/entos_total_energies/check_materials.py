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


