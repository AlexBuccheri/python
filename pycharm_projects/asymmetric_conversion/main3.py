
# External libraries
import numpy as np
from ase.io import read, write
from ase.atoms import Atoms, Atom
from ase import spacegroup as ase_spacegroup
import spglib

# My modules
from modules.spglib.cell import spglib_to_ase
from modules.ase.spglib import ase_to_spglib
from modules.spglib.io import show_cell as spg_show_cell, write as spg_write
from modules.electronic_structure.structure import supercell
from modules.electronic_structure.structure import atoms
from modules.parameters.elements import an_to_symbol
from modules.fileio.write import xyz as alex_xyz
# -------------------
# Main
# -------------------

structure_name = 'aei'
input_dir = 'inputs'
output_dir = 'aei_outputs'

ase_input_data = read(input_dir+"/"+structure_name+".cif", store_tags=False)
spg_input = ase_to_spglib(ase_input_data)
print("Number of atoms in input", len(spg_input[2]))

# Reduce to primitive
print(" Find primitive of conventional structure")
lattice, positions, numbers = spglib.find_primitive(spg_input, symprec=1e-5)
spg_molecule = (lattice, positions, numbers)
spg_show_cell(lattice, positions, numbers)
spg_write(output_dir + '/' + structure_name + '_primtive_cell.xyz', spg_molecule, pbc=(1,1,1))

# -----------------------------------------------
# Confirm the supercell can recovered
# from applying lattice vectors to primitive cell
# -----------------------------------------------

# Need lattice vectors column-wise in np array
lattice_vectors = np.zeros(shape=(3,3))
for i in range(0,3):
    lattice_vectors[:,i] = lattice[i]

# Need positions in angstrom, not fractional
positions_ang = []
for position in positions:
    positions_ang.append(np.matmul(lattice, position))

# Need unit cell in my molecule format
species = [an_to_symbol[an] for an in numbers]
unit_cell = atoms.Atoms(species=species, positions=positions_ang)
n = [2,2,2]
translations = supercell.translation_vectors(lattice_vectors, n, centred_on_zero = False)
s_cell = supercell.build_supercell(unit_cell, translations)

# Output it in xyz
alex_xyz(output_dir + '/' + "aei_supercell", s_cell)



# Modify rings by doing one B-O-B
# Structurally relax the lattice maintaining bond angles
# OR
# Replace around the entirety of each ring before relaxing

# Once one has the supercell, could also randomly apply the snip (whilst retaining global structural connectivity)