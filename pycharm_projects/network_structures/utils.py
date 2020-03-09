import numpy as np

import mendeleev

from modules.electronic_structure.structure import atoms, bravais


# Convert ASE atom data into spglib data.
# All tuples apart from atomic_numbers
#
def ase_atom_to_spg_atom(ase_data):

    lattice = []
    for vector in ase_data.cell:
        lattice.append(tuple(vector))

    # ASE converts basis to cartesian, so convert back
    inv_lattice = np.linalg.inv(np.transpose(np.asarray(ase_data.cell)))

    basis = []
    for pos in ase_data.positions:
        frac_pos = np.matmul(inv_lattice, pos)
        basis.append(tuple(frac_pos))

    atomic_numbers = []
    for an in ase_data.numbers:
        atomic_numbers.append(an)

    molecule = (lattice, basis, atomic_numbers)

    return molecule


# Use Python mendeleev to create a symbol:atomic-number map
# https://mendeleev.readthedocs.io/en/stable/data.html
# Only need to do once, then tabulate below
#
def atomic_number_symbol_dict():
    n_elements = 118
    an_symbol = {}
    for an in range(1, n_elements+1):
        element = mendeleev.element(an)
        an_symbol[an] = element.symbol
    return an_symbol

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
                100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt',
                110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}



# Convert ASE atom data into my atom data class
#
def ase_atom_my_atom(ase_data):

    lattice = []
    for vector in ase_data.cell:
        lattice.append(tuple(vector))

    # ASE converts basis to cartesian, so convert back
    inv_lattice = np.linalg.inv(np.transpose(np.asarray(ase_data.cell)))

    basis = []
    for pos in ase_data.positions:
        frac_pos = np.matmul(inv_lattice, pos)
        basis.append(tuple(frac_pos))

    atomic_numbers = []
    for an in ase_data.numbers:
        atomic_numbers.append(an)

    # an_to_symbol = atomic_number_symbol_dict()
    n_atoms = len(atomic_numbers)
    molecule = []

    for ia in range(0, n_atoms):
        species = an_to_symbol[atomic_numbers[ia]]
        molecule.append(atoms.Atom(species, basis[ia]))

    return molecule


# Distance matrix with my atom class
def distance_matrix(molecule):
    n_atoms = len(molecule)
    d = np.zeros(shape=(n_atoms,n_atoms))

    for ia in range(0,n_atoms):
        pos_i = molecule[ia].position
        for ja in range(ia, n_atoms):
            pos_j = molecule[ja].position
            d[ia,ja] = np.linalg.norm(pos_j - pos_i)
            d[ja,ia] = d[ia,ja]

    return d


# Find nearest neighbours of each atom from my Atom Class
#
def neighbour_list(molecule, neighbour_radius):
    n_atoms = len(molecule)
    neighbours = []
    tol = 1.e-5

    for ia in range(0, n_atoms):
        pos_i = molecule[ia].position
        neighbours_of_atom_i = []

        for ja in range(0, n_atoms):
            pos_j = molecule[ja].position
            separation = np.linalg.norm(pos_j - pos_i)
            if (separation <= neighbour_radius + tol) and (ia != ja):
                neighbours_of_atom_i.append(ja)
        neighbours.append(neighbours_of_atom_i)

    return neighbours


# Find nearest neighbours of each atom from ASE data
#
def neighbour_list_from_ase(ase_data, neighbour_radius):
    n_atoms = len(ase_data)
    neighbours = []
    tol = 1.e-5

    for ia in range(0, n_atoms):
        #atom_i = ase_data[ia]
        pos_i = ase_data[ia].position
        neighbours_of_atom_i = []
        for ja in range(0, n_atoms):
            #atom_j = np.asarray(ase_data[ja])
            pos_j = ase_data[ja].position
            separation = np.linalg.norm(pos_i - pos_j)
            if (separation <= neighbour_radius + tol) and (ia != ja):
                neighbours_of_atom_i.append(ja)
        neighbours.append(neighbours_of_atom_i)

    return neighbours


def find_corner_sharing_oxy(molecule, neighbours):
    corner_sharing_oxy = []
    for ia, atom in enumerate(molecule):
        if atom.species.lower() == 'o':
            silicons = neighbours[ia]
            # Oxygen corner-shares with two tetrahedra
            if len(silicons) == 2:
                corner_sharing_oxy.append(ia)
    return corner_sharing_oxy

