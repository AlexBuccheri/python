import numpy as np

from modules.electronic_structure.structure import atoms, bravais, crystals, supercell
from modules.electronic_structure.basis import gfn1

from modules.parameters import crystals as params
from modules.fileio import write
from modules.entos import input_strings

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
translations = supercell.translation_vectors(lattice, n, centred_on_zero=False)
super_cell = supercell.build_supercell(unit_cell, translations)



# Information
print('Number of atoms in supercell:', len(super_cell))
#print("Supercell lattice constants: ", , '(Ang)')  Write this


structure_string = input_strings.structure(super_cell, unit=atoms.CoordinateType.XYZ)
#print(structure_string)

xyz_string = write.xyz(super_cell)
print(xyz_string)

quit()


# Find central cell
icentre = supercell.flatten_supercell_limits([0,0,0], n, centred_on_zero=True)
# could return this from the build_supercell command

cells = supercell.list_global_atom_indices_per_cells(unit_cell, translations)
central_cell_atom_indices = cells[icentre]
print("Central cell contains atoms with this indices:", central_cell_atom_indices)
central_cell = []
for iatom in central_cell_atom_indices:
    central_cell.append(super_cell[iatom])

visualise_central_cell = True

# Show central cell has been found:
if visualise_central_cell:
    # Print out all atoms in xyz
    all_atoms_string = write.xyz(super_cell, 'silicon supercell')
    fid = open("all_atoms.xyz", "w")
    print(all_atoms_string, file=fid)
    fid.close()
    # Print out atoms of central cell only
    central_cell_string = write.xyz(central_cell, 'central cell in supercell')
    fid = open("central_cell.xyz", "w")
    print(central_cell_string, file=fid)
    fid.close()


# GFN1 basis sizes
basis_sizes = {'si': 'd'}
#ao = gfn1.make_basis(super_cell, basis_sizes)

# Find shell-resolved indices that correspond to atoms in central cell
print("These atoms have basis functions with shell indices:")
shell_indices = gfn1.find_basis(central_cell_atom_indices, super_cell, basis_sizes)

for i,atom in enumerate(shell_indices):
    icell = central_cell_atom_indices[i]
    print(icell, atom)



# Couple of tests:
# 1 atom species. spd 
# 2 atom species. one sp, the other spd
# Central cell and last cell, when centred on (0,0,0)







# Decide what to do with this stuff

# def atom_indices(lattice, translations, atomic_basis):
#     inv_lattice = np.linalg.inv(lattice)
#     atom_index = 0
#     for translation in translations:
#         n = np.matmul(inv_lattice, translation)
#         for ibasis, atom in enumerate(atomic_basis):
#             print(n, ibasis, atom_index)
#             atom_index += 1
#     return
#
# # To the atom indices so I can loop over shells or components
# def central_charges(n, atomic_basis, nshells):
#     central_cell_atoms = central_cell_indices(n, atomic_basis)
#

#structure_string = entos.input_strings.structure(super_cell, atoms.CoordinateType.XYZ)
