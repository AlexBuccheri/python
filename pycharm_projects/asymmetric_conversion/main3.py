""" In order to apply substitutions to the system, one needs a system where all atoms are connected.
    As such, I'm translating all loose atoms into positions in other unit cells BUT are connected as
    nearest neighbours to atoms in the central cell. One can then perform substitutions on this system
    and fold any positions that lie outside the central cell of the final structure, back inside
"""

# External libraries
import numpy as np
from ase.io import read, write
from ase.atoms import Atoms, Atom
from ase import spacegroup as ase_spacegroup
import spglib
from scipy import spatial

# My modules
from modules.spglib.cell import spglib_to_ase
from modules.ase.spglib import ase_to_spglib
from modules.spglib.io import show_cell as spg_show_cell, write as spg_write
from modules.electronic_structure.structure import supercell
from modules.electronic_structure.structure import atoms
from modules.parameters.elements import an_to_symbol
from modules.fileio.write import xyz as alex_xyz
import cell_operations

# Ref: http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
def flatten_list(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

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
#spg_show_cell(lattice, positions, numbers)
spg_write(output_dir + '/' + structure_name + '_primtive_cell.xyz', spg_molecule, pbc=(1,1,1))

# --------------------------------------------------------
# Create supercell such that the central primitive cell
# is fully coordinated
# -------------------------------------------------------

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
# [-1:1] per dimension, centred on zero
n = [3, 3, 3]
translations = supercell.translation_vectors(lattice_vectors, n, centred_on_zero=True)

n_atoms_prim = len(unit_cell)
n_atoms_super = np.product(n) * n_atoms_prim

print("N atoms in primitive cell ", n_atoms_prim)
print("N atoms in supercell: ",n_atoms_super)


# Find atoms neighbouring central cell
#coordinating_atoms = cell_operations.find_atoms_neighbouring_central_cell(unit_cell, translations, visualise=True)

# Replace uncoordinated atoms in central cell
unit_cell = cell_operations.replace_uncoordinated_atoms(unit_cell, translations)

# -------------------------------------------------------
# Test folding back in: Do this with py test
# -------------------------------------------------------
def position_in_cell(pos):
    """ For positions outside of the unit cell,
       fold them back in  """

    indices = np.where((pos < 0) | (pos > 1))[0]
    if len(indices) == 0:
        return pos
    else:
        for i in indices:
            pos[i] = pos[i] - np.floor(pos[i])

    return pos


inv_lattice = np.linalg.inv(lattice)
for ia,atom in enumerate(unit_cell):
    fractional_position = np.matmul(inv_lattice, atom.position)
    fractional_position = position_in_cell(fractional_position)
    unit_cell[ia].position = np.matmul(lattice, fractional_position)
    #print(atom.species, unit_cell[ia].position)

alex_xyz(output_dir + '/' + "aei_primitive_cell_folded_back", unit_cell)
print("folded atomic positions back in: aei_primitive_cell_folded_back",)


# -------------------------------------------------------
# Other ideas/notes
# -------------------------------------------------------
# Modify rings by doing one B-O-B
# Structurally relax the lattice maintaining bond angles
# OR
# Replace around the entirety of each ring before relaxing

# Once one has the supercell, could also randomly apply the snip (whilst retaining global structural connectivity)