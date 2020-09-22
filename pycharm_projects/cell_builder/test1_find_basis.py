import numpy as np

from modules.electronic_structure.structure import atoms, bravais, crystals, supercell
from modules.electronic_structure.basis import gfn1

from modules.parameters import crystals as params
from modules.fileio import write


# Silicon cubic unit cell in angstrom
al = params.silicon['lattice_constant']['angstrom']
fractional_basis_positions = crystals.silicon('conventional')
lattice = bravais.simple_cubic(al)

unit_cell = []
for atom in fractional_basis_positions:
    pos_angstrom = np.matmul(lattice, atom.position)
    unit_cell.append(atoms.Atom(atom.species, pos_angstrom))


# Silicon supercell
# If odd, can centre on T=(0,0,0), else can't
n = [3, 3, 3]
translations = supercell.translation_vectors(lattice, n, centred_on_zero=True)
super_cell = supercell.build_supercell(unit_cell, translations)


icentre = supercell.flatten_supercell_limits([0,0,0], n, centred_on_zero=True)
cells = supercell.list_global_atom_indices_per_cells(unit_cell, translations)
central_cell_atom_indices = cells[icentre]

print("Central cell contains atoms with this indices:", central_cell_atom_indices)
central_cell = []
for iatom in central_cell_atom_indices:
    central_cell.append(super_cell[iatom])


#################################################################
# Test. For one atomic species, compare function output to
# simple, analytic expression
#################################################################

# ----------------
# Approach 1
# ----------------
# GFN1 basis sizes
basis_sizes = {'si': 'd'}

# Find shell-resolved indices that correspond to atoms in central cell
#print("These atoms have basis functions with shell indices:")
shell_indices = gfn1.find_basis(central_cell_atom_indices, super_cell, basis_sizes)

# for i,atom in enumerate(shell_indices):
#     icell = central_cell_atom_indices[i]
#     print(icell, atom)


# ----------------
# Approach 2
# ----------------
# Works if shells_per_atom is fixed
def simple_function(central_cell_atom_indices, shells_per_atom):
    shell_indices = []
    for iatom in central_cell_atom_indices:
        initial_index = iatom * shells_per_atom
        shell_indices.append([i for i in range(initial_index, initial_index + shells_per_atom)])
    return shell_indices

shells_per_atom = 3
shell_indices2 = simple_function(central_cell_atom_indices, shells_per_atom)

assert(len(shell_indices) == len(shell_indices2))

are_same = []
for i in range(0, len(shell_indices2)):
    are_same.append(shell_indices[i] == shell_indices2[i])

if all(are_same) == True:
    print("Test passed, all shell indices are the same using two different methods")


#######################################
# Test 2. Non-central cell
#######################################
icentre = supercell.flatten_supercell_limits([1,1,0], n, centred_on_zero=True)
cells = supercell.list_global_atom_indices_per_cells(unit_cell, translations)
central_cell_atom_indices = cells[icentre]

shell_indices = gfn1.find_basis(central_cell_atom_indices, super_cell, basis_sizes)
shell_indices2 = simple_function(central_cell_atom_indices, shells_per_atom)

are_same = []
for i in range(0, len(shell_indices2)):
    are_same.append(shell_indices[i] == shell_indices2[i])

if all(are_same) == True:
    print("Test passed, all shell indices are the same using two different methods")