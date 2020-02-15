import numpy as np

# GFN1 minimal atomic orbital basis

angular_momentum = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

# Create a basis map. Note, basis in shells
class Basis:
    def __init__(self, atom_species, principal_qn, angular_momentum):
        self.species = atom_species
        self.n = principal_qn
        self.l = angular_momentum

# Currently a dummy routine
def get_principal_qn(atom, l):
    return 1

# Make a minimal atomic orbital basis, resolved at the level of shells
def make_basis(molecule, basis_sizes):
    atomic_basis = []
    for atom in molecule:
        l_max_symbol = basis_sizes[atom.species.lower()]
        l_max = angular_momentum[l_max_symbol]
        for l in range(0, l_max+1):
            n = get_principal_qn(atom, l)
            atomic_basis.append(Basis(atom.species, n, l))
            # Probably some other info one should add
    return atomic_basis


# For a given list of atoms, return the start and end basis
# functions that correspond to each of them. Shell-resolved
# NOTE: Assumes basis is ordered shells per atom i.e.
# atom 1   atom 2   atom 3 ...
# s,p,d ,  s,p,d ,  s,p,d ...
#
def find_basis(atom_indices, molecule, basis_sizes):
    shell_indices = []

    for iatom in atom_indices:
        first_index = 0
        # Find first shell index for atom 'iatom'
        for ia in range(0, iatom):
            species = molecule[ia].species
            max_l = angular_momentum[basis_sizes[species.lower()]]
            first_index += max_l+1
        #All shell indices corresponding to iatom
        max_l = angular_momentum[basis_sizes[molecule[iatom].species.lower()]]
        shell_indices.append([i for i in range(first_index, first_index+max_l+1)])

    return shell_indices