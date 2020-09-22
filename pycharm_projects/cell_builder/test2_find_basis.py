import numpy as np

from modules.electronic_structure.structure import atoms, bravais, crystals, supercell
from modules.electronic_structure.basis import gfn1

from modules.parameters import crystals as params
from modules.fileio import write


# Silicon dioxide cubic primitive unit cell. Could be wrong but doesn't matter, just want two atom types
fractional_basis_positions = [atoms.Atom('Si', [0, 0, 0]), atoms.Atom('O', [0.25, 0.25, 0.25])]
lattice = bravais.face_centred_cubic(1.)

unit_cell = []
for atom in fractional_basis_positions:
    pos_angstrom = np.matmul(lattice, atom.position)
    unit_cell.append(atoms.Atom(atom.species, pos_angstrom))


# Supercell
# If odd, can centre on T=(0,0,0), else can't
n = [3, 3, 3]
translations = supercell.translation_vectors(lattice, n, centred_on_zero=True)
super_cell = supercell.build_supercell(unit_cell, translations)


# Central cell
icentre = supercell.flatten_supercell_limits([0,0,0], n, centred_on_zero=True)
cells = supercell.list_global_atom_indices_per_cells(unit_cell, translations)
central_cell_atom_indices = cells[icentre]


#################################################################
# Test. Can compare to the simple, analytic expression if one
# counts the number of si + number of
#################################################################

# ----------------
# Approach 1
# ----------------
# GFN1 basis sizes - again, made up
basis_sizes = {'si': 'd', 'o': 'p'}

# Find shell-resolved indices that correspond to atoms in central cell
print("These atoms have basis functions with shell indices:")
shell_indices = gfn1.find_basis(central_cell_atom_indices, super_cell, basis_sizes)


# ----------------
# Approach 2
# ----------------
# Pointless - end up writing an equally non-trivial function to test this