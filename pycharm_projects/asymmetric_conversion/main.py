#!/usr/bin/env python3

# --------------------------------------------------------------------------------
# Convert conventional cell or primitive cell to asymmetric cell
# Alex Buccheri 2020
#
# --------------------------------------------------------------------------------

from ase.io import read, write
from ase.atoms import Atoms, Atom
from ase import spacegroup as ase_spacegroup
from ase.visualize import view as ase_view
import numpy as np
import spglib


# -------------------------------------------
# Functions
# -------------------------------------------

an_to_symbol = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
                10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K',
                20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
                30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y',
                40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In',
                50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr',
                60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm',
                70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
                80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac',
                90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
                100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs',
                109: 'Mt',
                110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}

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

# -----------------------------
# My functions
# -----------------------------
def asymmetric_cell_atom_indices(dataset):
    # atomic indices in supercell
    atom_indices = []
    for x in dataset['equivalent_atoms']:
        atom_indices.append(x)

    # Reduce to unique indices
    atom_indices = list(set(atom_indices))
    atom_indices.sort()
    return atom_indices



def group_reducible_atomic_indices(dataset):
    sets_of_reducibles = []
    one_set = []
    ref_atom = dataset['equivalent_atoms'][0]

    for reducible_atom, i in enumerate(dataset['equivalent_atoms']):

        if (i != ref_atom):
            sets_of_reducibles.append(one_set)
            one_set = []
            ref_atom = reducible_atom

        one_set.append(reducible_atom)

    return sets_of_reducibles


def nearest_neighbour_asymmetric_cell_atom_indices(dataset, molecule):
    sets_of_reducibles = group_reducible_atomic_indices(dataset)
    irreducible_atoms = [sets_of_reducibles[0][0]]
    #Remove first list having initialised with it
    sets_of_reducibles = sets_of_reducibles[1:]
    position = molecule[1]

    for set in sets_of_reducibles:
        separation = np.zeros(shape=(len(set)))

        for i, iatom in enumerate(set):
            pos_i = np.asarray(position[iatom])

            for jatom in irreducible_atoms:
                pos_j = np.asarray(position[jatom])
                separation[i] += np.linalg.norm(pos_j - pos_i)
        imin = np.argmin(separation)

        # Reducible atom from current set with the lowest average seperation from all
        # atoms currently in irreducible_atoms, is added as the next irreducible atom
        irreducible_atoms.append(set[imin])

    return irreducible_atoms


def ase_to_spglib(ase_data):
    # Get ASE into spglib format. All tuples apart from atomic_numbers
    lattice = []
    for vector in ase_data.cell:
        lattice.append(tuple(vector))

    # ASE converts basis to cartesian, so convert back to fractional
    inv_lattice = np.linalg.inv(np.transpose(np.asarray(ase_data.cell)))

    basis = []
    for pos in ase_data.positions:
        frac_pos = np.matmul(inv_lattice, pos)
        basis.append(tuple(frac_pos))

    atomic_numbers = []
    for an in ase_data.numbers:
        atomic_numbers.append(an)

    # SGLIB data storage for a crystal
    molecule = (lattice, basis, atomic_numbers)
    return molecule


def spglib_to_ase(molecule, indices=None):
    basis = molecule[1]
    atomic_numbers = molecule[2]
    if indices == None:
        indices = range(0, len(atomic_numbers))

    lattice = np.transpose(np.asarray(molecule[0]))
    #print(lattice)

    ase_molecule = []
    for ia in indices:
        atomic_symbol = an_to_symbol[atomic_numbers[ia]]
        # Have to store in Cartesian
        pos = np.matmul(lattice, np.asarray(basis[ia]))
        # Fractional
        #pos = basis[ia]
        #print(ia, basis[ia])
        ase_molecule.append(Atom(atomic_symbol, pos))

    return Atoms(ase_molecule, cell=molecule[0])

# Print symmetry data
# Not in my prefered notation:
# https://en.wikipedia.org/wiki/Crystallographic_point_group#Hermann–Mauguin_notation
def print_spg_symmetry_info(dataset, wyckoff=False, equivalent_atoms=False):
    print("  Spacegroup is %s (%d)." % (dataset['international'], dataset['number']))
    print("  Pointgroup is %s." % (dataset['pointgroup']))
    print("  Hall symbol is %s (%d)." % (dataset['hall'], dataset['hall_number']))
    if wyckoff:
        print("  Wyckoff letters are: ", dataset['wyckoffs'])
    if equivalent_atoms:
        print("  Mapping to equivalent atoms are: ")
        for i, x in enumerate(dataset['equivalent_atoms']):
            print("  %d -> %d" % (i + 1, x + 1))
    return


# Given an spg_molcule and set of atomic indices, return a new spg module
def create_spg_molecule(input_molecule, indices):
    lattice        = input_molecule[0]
    basis          = input_molecule[1]
    atomic_numbers = input_molecule[2]
    new_basis = []
    new_atomic_numbers = []

    for iatom in indices:
        new_basis.append(tuple(basis[iatom]))
        new_atomic_numbers.append(atomic_numbers[iatom])

    return (lattice, new_basis, new_atomic_numbers)


# https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
# TODO(Alex) Add some asserts for when this won't work
def rotation_to_align_a_with_b(a, b):
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if not np.allclose(a, a/norm_a):
        print('Input a vector not unit normal - normalising')
        a = a / norm_a
        print(a)
    if not np.allclose(b, b/norm_b):
        print('Input b vector not unit normal - normalising')
        b = b / norm_b
        print(b)

    v = np.cross(a,b)
    #s = np.linalg.norm(v)
    c = np.dot(a,b)
    f = 1./(1. + c)
    vmat = np.array([[    0, -v[2],  v[1]],
                     [ v[2],     0, -v[0]],
                     [-v[1],  v[0],     0]])
    return np.eye(3,3) + vmat + f *(np.matmul(vmat,vmat))


# -----------------------------------
# Main Routine
# -----------------------------------

#Read CIF with ASE and convert to SPG format
ase_input_data = read("aei.cif", store_tags=False)
#print(vars(ase_input_data))
spg_input = ase_to_spglib(ase_input_data)
print("Number of atoms in input", len(spg_input[2]))


# Reduce to primitive cell
find_primitive = True
if find_primitive:
    print(" Find primitive of conventional structure")
    lattice, positions, numbers = spglib.find_primitive(spg_input, symprec=1e-1)
    show_cell(lattice, positions, numbers)
    print("Number of atoms in primitive: ", len(numbers))
    spg_molecule = (lattice, positions, numbers)
    ase_primitive_cell = spglib_to_ase(spg_molecule)
    ase_primitive_cell.set_pbc((1, 1, 1))
    write('aei_primtive_cell.xyz', ase_primitive_cell)
else:
    spg_molecule = spg_input

# Get and print symmetry data
dataset = spglib.get_symmetry_dataset(spg_molecule)
print_spg_symmetry_info(dataset, equivalent_atoms=False, wyckoff=False)


# Reduce cell/primitive cell to asymmetric cell of irreducible atomic positions
neighbouring_atoms = True
if neighbouring_atoms:
    irreducible_atom_indices = nearest_neighbour_asymmetric_cell_atom_indices(dataset, spg_molecule)
else:
    irreducible_atom_indices = asymmetric_cell_atom_indices(dataset)


# NOTE: Don't appear to be able to find the symmetry of the asymmetric cell
#spg_molecule2 = create_spg_molecule(spg_molecule, irreducible_atom_indices)
#dataset2 = spglib.get_symmetry_dataset(spg_molecule2)
#print_spg_symmetry_info(dataset2, equivalent_atoms=False)


# Convert to ASE and write to xyz
# Also need to write as CIF
ase_asymmetric_cell = spglib_to_ase(spg_molecule, irreducible_atom_indices)
ase_asymmetric_cell.set_pbc((1, 1, 1))
write('aei_asymmetric_cell.xyz', ase_asymmetric_cell)


# TEST
# Use ASE to apply symmetry operations and reproduce original cell from asymmetric cell
# THIS WORKS
#aei_cell = ase_spacegroup.crystal(ase_asymmetric_cell, spacegroup=dataset['number'])
#ase_view(aei_cell)
#
# Test; Rotate whole system, output and see if it looks sensible
# Awesome, this works
# O (atom5)  - si (atom1)
# a = ase_asymmetric_cell.positions[4] - ase_asymmetric_cell.positions[0]
# b = [0,0,1]
# R = rotation_to_align_a_with_b(a, b)
# new_positions = []
# for pos in ase_asymmetric_cell.positions:
#     new_pos = np.matmul(R, pos)
#     new_positions.append(new_pos)
# ase_asymmetric_cell.positions = new_positions
# write('rotated_aei_asymmetric_cell.xyz', ase_asymmetric_cell)


# Scale bond lengths oxygens with one silicon neighbour
def scale_bondlength_of_oxy_with_one_neighbour(pair_indices, ase_cell, bond_length):
    z_unit = [0, 0, 1]

    for pair in pair_indices:
        iOxy = pair['o']
        iSi = pair['si']

        # Rotate system so O-Si displacement vector points in +z direction
        disp_o_si = ase_cell.positions[iOxy] - ase_cell.positions[iSi]
        R = rotation_to_align_a_with_b(disp_o_si, z_unit)

        new_ox_position = np.matmul(R, ase_cell.positions[iSi]) + np.array([0, 0, bond_length])
        # Transform back to original reference, utilising R is unitary
        ase_cell.positions[iOxy] = np.matmul(np.transpose(R),new_ox_position)
    return ase_cell


def scale_bondlength_of_oxy_with_two_neighbours(pair_indices, ase_cell, bond_length):
    z_unit = [0, 0, 1]

    for pair in pair_indices:
        iOxy = pair['o']
        iSi = pair['si']

        #Si above O in z
        disp_si_o = ase_cell.positions[iSi] - ase_cell.positions[iOxy]

        #Rotate whole structure
        #For all atoms above iOxy, shift them
        # SEE HARD NOTES


    return ase_cell


def find_oxy_with_one_neighbour():
    return

def find_oxy_with_two_neighbours():
    return

def convert_element(element_a, element_b):
    return


# Take an asymmetric cell and convert from a silicate to a boron oxide
def convert_to_boron_oxide(ase_cell):

    bond_length = 2.

    #Find all oxy with on si neighbour and store like so:
    pair_indices = [{'o': 4, 'si': 0}]
    scale_bondlength_of_oxy_with_one_neighbour(pair_indices, ase_cell, bond_length)

    # Scale bond lengths oxygens with two silicon neighbours
    # => translating everything else attached to that oxygen

    # Convert silicons to borons
    # This is trivial to do

    return ase_cell

ase_asymmetric_cell = convert_to_boron_oxide(ase_asymmetric_cell)
write('scaled_aei_asymmetric_cell.xyz', ase_asymmetric_cell)


# Apply ASE symmetry operations to convert into supercell

# Run in MM package, GFN0 and DFTB+: See what the structure relaxes to

