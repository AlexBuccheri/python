
# Functions for working with spglib 'molecule' and 'dataset' objects 

import numpy as np
from ase.atoms import Atoms as ase_Atoms, Atom as ase_Atom

from modules.parameters.elements import an_to_symbol


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
        #print(ia, basis[ia])
        ase_molecule.append(ase_Atom(atomic_symbol, pos))

    return ase_Atoms(ase_molecule, cell=molecule[0])



# Given an spg_molcule and set of atomic indices, return a new spg module
def create_spg_molecule(input_molecule, indices):
    lattice        = input_molecule[0]
    old_positions      = input_molecule[1]
    old_atomic_numbers = input_molecule[2]
    positions = []
    atomic_numbers = []

    origin = [0,0,0]
    #origin = np.asarray(old_positions[0])

    for iatom in indices:
        positions.append(tuple(old_positions[iatom] - origin))
        atomic_numbers.append(old_atomic_numbers[iatom])

    return (lattice, positions, atomic_numbers)


# dataset is returned from calling spglib.get_symmetry_dataset(spg_molecule)
def asymmetric_cell_atom_indices(dataset):
    # atomic indices in supercell
    atom_indices = []
    for x in dataset['equivalent_atoms']:
        atom_indices.append(x)

    # Reduce to unique indices
    atom_indices = list(set(atom_indices))
    atom_indices.sort()
    return atom_indices


# dataset is returned from calling spglib.get_symmetry_dataset(spg_molecule)
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
