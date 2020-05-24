# Check lattice parameters are consistent with entos lattice vectors

import numpy as np
import spglib

from modules.electronic_structure.structure import bravais
from modules.electronic_structure.structure import lattice as lt
from modules.spglib import io as spglib_io

def find_point_group(cell, print_out = False):
    dataset = spglib.get_symmetry_dataset(cell, symprec=1e-5, angle_tolerance=-1.0)
    if print_out: spglib_io.show_spg_symmetry_info(dataset)
    return dataset['pointgroup']

def show_lattice(lattice):
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        print("%s %10.5f %10.5f %10.5f" % (tuple(axis, ) + tuple(vec)))

def show_cell(lattice, positions, numbers):
    show_lattice(lattice)
    print("Atomic points:")
    for p, s in zip(positions, numbers):
        print("%2d %10.5f %10.5f %10.5f" % ((s,) + tuple(p)))



# cell volume and cell angles -> Should definitely add this to the entos check

# Anatase
rad_to_deg = np.pi / 180
lattice = bravais.simple_triclinic(a=5.55734663, b = 5.55734663, c=5.55734663,
                                   alpha = 139.94858990 * rad_to_deg,
                                   beta  = 139.94858990 * rad_to_deg,
                                   gamma = 57.93136581  * rad_to_deg)
V = lt.parallelpiped_volume(lattice)
assert np.isclose(V, 70.436705, atol=1.e-6)

lattice2 = bravais.body_centred_tetragonal(a=5.55734663, c=5.55734663)
assert not np.allclose(lattice, lattice2)
print("These anatase parameters do not appear to be for body-centred tetragonal")

check_anatase = False
if check_anatase:

    print("Initial lattice vectors, stored columnwise")
    print(lattice2)
    print("Anatase: Use SPG Lib to convert to primitive vectors and positions")

    def sglib_data_types(columnwise_lattice):
        lattice = np.transpose(columnwise_lattice)
        positions = [[0.500000, 0.500000, 0.000000],
                     [0.250000, 0.750000, 0.500000],
                     [0.456413, 0.956413, 0.500000],
                     [0.706413, 0.706413, 0.000000],
                     [0.043587, 0.543587, 0.500000],
                     [0.293587, 0.293587, 0.000000]]
        numbers = [81, 81, 8, 8, 8, 8]
        cell = (lattice, positions, numbers)
        return cell


    cell = sglib_data_types(lattice2)
    lattice, positions, numbers = spglib.find_primitive(cell, symprec=1e-1)
    primitive_cell = (lattice, positions, numbers)
    show_cell(lattice, positions, numbers)
    point_group = find_point_group(primitive_cell, True)

    # Switch back to column-wise lattice vectors
    lattice = np.transpose(lattice)
    a = -2 * lattice[0, 0]
    c = -2 * lattice[2, 2]
    print("Lattice constants: ", a, c)
    (alpha, beta, gamma) = lt.cell_angles(lattice, return_unit='degrees')
    print("Primitive cell angles:", (alpha, beta, gamma))

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


# Calcium titanate: CaTiO3
# Structure's consistent but entos bukl calculation fails
lattice = bravais.simple_orthorhombic(a=5.40444906, b=5.51303112, c=7.69713264)
V = lt.parallelpiped_volume(lattice)
assert np.isclose(V, 229.335256, atol=1.e-6)

(alpha, beta, gamma) = lt.cell_angles(lattice, return_unit='degrees')
assert np.allclose(alpha, 90.0)
assert np.allclose(beta,  90.0)
assert np.allclose(gamma, 90.0)




#Is my base-centred monoclinic consistent with the paper??????
# lattice = bravais.base_centred_monoclinic(a = 6.41420528, b = 6.41420528, c=5.87615337, beta=103.27512313)
#
# (alpha, beta, gamma) = lt.cell_angles(lattice, return_unit='degrees')
# assert np.allclose(alpha, 76.72487687)
# assert np.allclose(beta, 103.27512313)
# assert np.allclose(gamma, 152.18872933)

