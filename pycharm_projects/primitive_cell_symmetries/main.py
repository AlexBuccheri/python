# ------------------------------------------------------------------
# Determine the point groups of the reciprocal primitive unit cells
# ------------------------------------------------------------------

import numpy as np
import spglib

from modules.electronic_structure.structure import bravais
from modules.electronic_structure.structure import lattice
from modules.spglib import io as spglib_io

def sglib_data_types(columnwise_lattice):
    lattice = np.transpose(columnwise_lattice)
    # Dummy basis and atomic number as I only care about the point group
    # Should not affect symmetry by placing it at the origin
    positions = [[0, 0, 0]]
    numbers = [1]
    cell = (lattice, positions, numbers)
    return cell


def find_point_group(lattice, print_out = False):
    cell = sglib_data_types(lattice)
    dataset = spglib.get_symmetry_dataset(cell, symprec=1e-5, angle_tolerance=-1.0)
    if print_out: spglib_io.show_spg_symmetry_info(dataset)
    return dataset['pointgroup']

# Not the smartest implementation
def check_all_equal(lst: list, string: str):
    array = np.asarray(lst)
    indices = np.where(array == string)[0]
    return len(indices) == len(lst)

def cubic_point_groups():
    point_groups = []
    al = 1

    cubic_lattice = bravais.simple_cubic(al)
    pg = find_point_group(cubic_lattice)
    point_groups.append(pg)
    print(' Real lattice of simple cubic:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(cubic_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of simple cubic:', pg)

    bcc_lattice = bravais.body_centred_cubic(al)
    pg = find_point_group(bcc_lattice)
    point_groups.append(pg)
    print(' Real lattice of BCC:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(bcc_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of BCC:', pg)

    fcc_lattice = bravais.face_centred_cubic(al)
    pg = find_point_group(fcc_lattice)
    point_groups.append(pg)
    print(' Real lattice of FCC:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(fcc_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of FCC:', pg, "\n")

    assert check_all_equal(point_groups, 'm-3m')
    return


def tetragonal_point_groups():
    point_groups = []
    a = 1
    c = 0.25

    tet_lattice = bravais.simple_tetragonal(a, c)
    pg = find_point_group(tet_lattice)
    point_groups.append(pg)
    print(' Real lattice of simple tetragonal:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(tet_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of simple tetragonal:', pg)

    bc_tet_lattice = bravais.body_centred_tetragonal(a, c)
    pg = find_point_group(bc_tet_lattice)
    point_groups.append(pg)
    print(' Real lattice of body-centred tetragonal:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(bc_tet_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of body-centred  tetragonal:', pg, "\n")

    assert check_all_equal(point_groups, '4/mmm')

    return


def trigonal_point_groups():
    return

def orthorhombic_point_groups():
    return

def monoclinic_point_groups():
    return

def triclinic_point_group():
    return

# Test for ENTOS
# Can test this by uniformly sampling the real-space primitive unit cell and apply the Wyckoff operations
# Repeat for the reciprocal space cell
# Confirm the point groups agree

# NOTE, NOT of the lattice + basis. Hence the dummy basis position at the origin
# In which case, this should be as given on wikipedia:
# https://en.wikipedia.org/wiki/Crystallographic_point_group#Hermann–Mauguin_notation

# Wyckoff Positions 3D Crystallographic Point Groups
# https://www.cryst.ehu.es/cryst/point_wp.html

print("Point groups of Bravais Crystal Lattices \n")

print("Cubic lattices:")
cubic_point_groups()

print("Tetragonal lattices:")
tetragonal_point_groups()

