
# External libraries
from ase.io import read, write
from ase.atoms import Atoms, Atom
from ase import spacegroup as ase_spacegroup
import spglib

# My modules
from modules.spglib.cell import spglib_to_ase
from modules.ase.spglib import ase_to_spglib
from modules.spglib.io import show_cell as spg_show_cell, write as spg_write

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

# Confirm the supercell can recovered from applying lattice vectors to primtive

# Modify rings by doing one B-O-B
# Structurally relax the lattice maintaining bond angles
# OR
# Replace around the entirety of each ring before relaxing

# Once one has the supercell, could also randomly apply the snip (whilst retaining global structural connectivity)