# -----------------------------------------------------------------------
# Determine the point groups of the reciprocal primitive crystal lattices
# Shown that in all cases, the point group of the real-space lattice
# is equal to the point group of the reciprocal-space lattice
# -----------------------------------------------------------------------

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
    point_groups = []
    a = 1
    c = 1.2
    alpha = 130 * (np.pi/180.)

    hex_lattice = bravais.hexagonal(a, c)
    pg = find_point_group(hex_lattice)
    point_groups.append(pg)
    print(' Real lattice of hexagonal:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(hex_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of hexagonal:', pg)

    #TODO(Alex) This is questionable: 2/m looks like monoclinic on wiki
    rhom_hex_setting_lattice = bravais.rhombohedral_hex_setting(a, c)
    pg = find_point_group(rhom_hex_setting_lattice)
    point_groups.append(pg)
    print(' Real lattice of rhombohedral (hexagonal setting):', pg)

    #TODO(Alex) This is questionable: 2/m looks like monoclinic on wiki
    recip_lattice = lattice.reciprocal_lattice_vectors(rhom_hex_setting_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of rhombohedral (hexagonal setting):', pg)

    #TODO(Alex) This is questionable: 2/m looks like monoclinic on wiki
    rhom_rhom_setting_lattice = bravais.rhombohedral_hex_setting(a, alpha)
    pg = find_point_group(rhom_rhom_setting_lattice)
    point_groups.append(pg)
    print(' Real lattice of rhombohedral (rhom setting):', pg)

    #TODO(Alex) This is questionable: 2/m looks like monoclinic on wiki
    recip_lattice = lattice.reciprocal_lattice_vectors(rhom_rhom_setting_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of rhombohedral (rhom setting):', pg, '\n')

    return



def orthorhombic_point_groups():
    point_groups = []
    a = 1
    b = 2
    c = 3

    orth_lattice = bravais.simple_orthorhombic(a, b, c)
    pg = find_point_group(orth_lattice)
    point_groups.append(pg)
    print(' Real lattice of simple orthorhombic:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(orth_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of simple orthorhombic:', pg)

    orth_lattice = bravais.base_centred_orthorhombic_A(a, b, c)
    pg = find_point_group(orth_lattice)
    point_groups.append(pg)
    print(' Real lattice of base-centred orthorhombic A:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(orth_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of base-centred orthorhombic A:', pg)

    orth_lattice = bravais.base_centred_orthorhombic_C(a, b, c)
    pg = find_point_group(orth_lattice)
    point_groups.append(pg)
    print(' Real lattice of base-centred orthorhombic C:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(orth_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of base-centred orthorhombic C:', pg)

    orth_lattice = bravais.face_centred_orthorhombic(a, b, c)
    pg = find_point_group(orth_lattice)
    point_groups.append(pg)
    print(' Real lattice of face-centred orthorhombic:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(orth_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of face-centred orthorhombic:', pg, '\n')

    assert check_all_equal(point_groups, 'mmm')

    return

def monoclinic_point_groups():
    point_groups = []
    a = 1
    b = 2
    c = 3

    beta =  101 * (np.pi/180.)
    mono_lattice = bravais.simple_monoclinic(a,b,c,beta)
    pg = find_point_group(mono_lattice)
    point_groups.append(pg)
    print(' Real lattice of simple monoclinic:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(mono_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of simple monoclinic:', pg)

    mono_lattice = bravais.simple_monoclinic(a, b, c, beta)
    pg = find_point_group(mono_lattice)
    point_groups.append(pg)
    print(' Real lattice of simple monoclinic:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(mono_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of simple monoclinic:', pg, '\n')

    assert check_all_equal(point_groups, '2/m')

    return

def triclinic_point_group():
    point_groups = []

    a = 1
    b = 2
    c = 3
    alpha = 85 * (np.pi/180.)
    beta  = 91 * (np.pi/180.)
    gamma = 57 * (np.pi/180.)

    tric_lattice = bravais.simple_triclinic(a,b,c, alpha, beta, gamma)
    pg = find_point_group(tric_lattice)
    point_groups.append(pg)
    print(' Real lattice of triclinic:', pg)

    recip_lattice = lattice.reciprocal_lattice_vectors(tric_lattice)
    pg = find_point_group(recip_lattice)
    point_groups.append(pg)
    print(' Reciprocal lattice of triclinic:', pg)

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

# TODO(Alex) Rhombohedral point group disagrees with wiki, totally
print("Trigonal lattices:")
trigonal_point_groups()

print("Orthorhombic lattices:")
orthorhombic_point_groups()

# Note, 2/m has more symmetry than m
print("Monoclinic lattices:")
monoclinic_point_groups()

print("Triclinic lattices:")
triclinic_point_group()