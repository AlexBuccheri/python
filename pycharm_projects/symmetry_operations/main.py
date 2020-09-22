# Use ASE and  lib to convert between  asymmetric, primtivie and conventional unit cells

# !/usr/bin/env python3

from ase.io import read, write
from ase.atoms import Atoms
import numpy as np
import spglib

# -----------------------------------
# Data Input and Output
# -----------------------------------

# All formats ASE can read/write to: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
# CIF: https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#cif

# store_tags grabs alot of useful info. Should use this in production code
data_conventional = read("BaTiO3-Orth-190K.cif", store_tags=False)
print(vars(data_conventional))
print("")

# If sublattice translations are included in the cif data, this perhaps won't work.
data_primitive = read("BaTiO3-Orth-190K.cif", subtrans_included=False, primitive_cell=True)

# Primitive comes out looking erroneous, perhaps due to sublattice options being included in
# initial cif data
write("ase-prim-BaTiO3-Orth-190K.cif", data_primitive)
write("ase-conv-BaTiO3-Orth-190K.cif", data_conventional)

# To create a supercell
# Both xyz supercells come out looking correct but orientated in different directions
write("primitive-BaTiO3-Orth-190K.xyz", data_primitive.repeat(6))
write("conventional-BaTiO3-Orth-190K.xyz", data_conventional.repeat(6))

# Additional ASE notes
# For the CIF format, STAR extensions as save frames, global blocks,
# nested loops and multi-data values are not supported.  Furthermore,
# ASE currently assumes the ``loop_`` identifier, and the following
# loop variable names to be on separate lines.
# https://wiki.fysik.dtu.dk/ase/ase/io/io.html
# https://gitlab.com/ase/ase/issues/15


# -----------------------------------
# Data extraction
# -----------------------------------
print("Looks like cell lengths are angstrom but it does not specify")

print("Lattice vectors, stored column-wise:")
lattice = np.transpose(np.asarray(data_conventional.cell))
inv_lattice = np.linalg.inv(lattice)
print(lattice)

# print("Atomic basis for conventional cell (assume ang), stored column-wise:")
n_atoms = len(data_conventional.numbers)
conventional_basis = np.zeros(shape=(3, n_atoms))
for i, pos in enumerate(data_conventional.positions):
    conventional_basis[:, i] = np.asarray(pos)

print("Atomic basis for conventional cell (fractional), stored column-wise:")
fractional_positions = np.matmul(inv_lattice, conventional_basis)
for iatom in range(0, n_atoms):
    print(fractional_positions[:, iatom])

# --------------------------------------------
# Convert conventional cell to primitive cell
# --------------------------------------------
print("")


# https://github.com/atztogo/spglib/blob/master/python/examples/example.py
def show_symmetry(symmetry, n_symmetries=None):
    if n_symmetries == None:
        for i in range(symmetry['rotations'].shape[0]):
            print("  --------------- %4d ---------------" % (i + 1))
            rot = symmetry['rotations'][i]
            trans = symmetry['translations'][i]
            print("  rotation:")
            for x in rot:
                print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
            print("  translation:")
            print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))
    else:
        for i in range(0, n_symmetries):
            print("  --------------- %4d ---------------" % (i + 1))
            rot = symmetry['rotations'][i]
            trans = symmetry['translations'][i]
            print("  rotation:")
            for x in rot:
                print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
            print("  translation:")
            print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))


def show_lattice(lattice):
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        print("%s %10.5f %10.5f %10.5f" % (tuple(axis, ) + tuple(vec)))


def show_cell(lattice, positions, numbers):
    show_lattice(lattice)
    print("Atomic points:")
    for p, s in zip(positions, numbers):
        print("%2d %10.5f %10.5f %10.5f" % ((s,) + tuple(p)))


# Read cif data with ASE
data_conventional = read("Si_mp-149_conventional_standard.cif")
print(vars(data_conventional))
print("")

# Get ASE into spglib format. All tuples apart from atomic_numbers
lattice = []
for vector in data_conventional.cell:
    lattice.append(tuple(vector))

# ASE converts basis to cartesian, so convert back
inv_lattice = np.linalg.inv(np.transpose(np.asarray(data_conventional.cell)))

basis = []
for pos in data_conventional.positions:
    frac_pos = np.matmul(inv_lattice, pos)
    basis.append(tuple(frac_pos))

atomic_numbers = []
for an in data_conventional.numbers:
    atomic_numbers.append(an)

print("SGLIB data storage for a crystal:")
molecule = (lattice, basis, atomic_numbers)
print(molecule)
print("  Spacegroup of conventional silicon is %s." % spglib.get_spacegroup(molecule))

symmetry = spglib.get_symmetry(molecule)
# show_symmetry(symmetry)
print("  Number of symmetry operations of conventional silicon is %d." % len(symmetry['rotations']))

# Only needs one rotation operation to determine the point group
# Not in my prefered notation. https://en.wikipedia.org/wiki/Crystallographic_point_group#Hermann–Mauguin_notation
print("  Pointgroup of silicon is %s." %
      spglib.get_pointgroup(symmetry['rotations'])[0])

# From the international tables. Not sure how this function differs to get_symmetry
# BUT appear to need this apporach to get Wyckoff symbols/operations
dataset = spglib.get_symmetry_dataset(molecule)
print("  Spacegroup of silicon is %s (%d)." % (dataset['international'],
                                               dataset['number']))
print("  Pointgroup of silicon is %s." % (dataset['pointgroup']))
print("  Hall symbol of silicon is %s (%d)." % (dataset['hall'],
                                                dataset['hall_number']))
print("  Wyckoff letters of silicon are: ", dataset['wyckoffs'])
print("  Mapping to equivalent atoms of silicon are: ")
for i, x in enumerate(dataset['equivalent_atoms']):
    print("  %d -> %d" % (i + 1, x + 1))

# Same as above
# print("  Symmetry operations of silicon unitcell are:")
# for i, (rot,trans) in enumerate(zip(dataset['rotations'],
#                                     dataset['translations'])):
#     print("  --------------- %4d ---------------" % (i + 1))
#     print("  rotation:")
#     for x in rot:
#         print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
#     print("  translation:")
#     print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))


print(" Fine primitive of conventional silicon structure")
lattice, positions, numbers = spglib.find_primitive(molecule, symprec=1e-1)
show_cell(lattice, positions, numbers)

# https://atztogo.github.io/spglib/python-spglib.html#niggli-reduce
# The detailed control of standardization of unit cell can be done using standardize_cell.

