
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
s_cell = supercell.build_supercell(unit_cell, translations)

# Output it in xyz
alex_xyz(output_dir + '/' + "aei_supercell", s_cell)

n_atoms_prim = len(unit_cell)
n_atoms_super = len(s_cell)
assert n_atoms_super == np.product(n) * n_atoms_prim

print("N atoms in primitive cell ", n_atoms_prim)
print("N atoms in supercell: ",n_atoms_super)


# -----------------------------------------------------------------------
# Smarter way to do this bit is just translate loose atoms by each translation
# vector and see if any positions put it as a NN
# This is the structure one works with
# Finally, once the boron structure has been constructed, one wraps all
# atoms outside the cell to inside
# -----------------------------------------------------------------------

unit_cell_positions = [atom.position for atom in unit_cell]
d = spatial.distance_matrix(unit_cell_positions, unit_cell_positions)

# Index loose atoms

radius = 1.8
loose_atoms = []
for ia in range(0, n_atoms_prim):
    indices = np.where((d[ia, :] > 0.) & (d[ia, :] <= radius))[0]
    if len(indices) == 0:
        loose_atoms.append(ia)

assert len(loose_atoms) == 4

# Should do this for each loose atom and reduce to a set
# target bond length is ~ 1.529 (I THINK) => Set lower bond too
def find_equivalent_position(ia, unit_cell, translations):
    pos_atom = unit_cell[ia].position
    loose_atom_positions = [pos_atom + translation for translation in translations]
    d_loose = spatial.distance_matrix(loose_atom_positions, unit_cell_positions)
    valid_equivalent = []
    for ja in range(0, len(loose_atom_positions)):
        indices = np.where((d_loose[ja, :] > 1.4) & (d_loose[ja, :] <= 1.8))[0]
        if len(indices) > 0:
            print("distance:", d_loose[ja, indices])
            valid_equivalent.append(loose_atom_positions[ja])
    assert len(valid_equivalent) == 1
    return valid_equivalent

valid_equivalents = []
for ia in loose_atoms:
    valid_equivalents.append(find_equivalent_position(ia, unit_cell, translations))


print(valid_equivalents)

neighbours = []
for position in valid_equivalents:
    print("in loop", position)
    neighbours.append(atoms.Atom(position=position[0].tolist(), species='B'))


unit_cell_and_neighbours = unit_cell + neighbours

for atom in unit_cell_and_neighbours:
    print(atom.species, atom.position)

alex_xyz(output_dir + '/' + "aei_central_NN_cell", unit_cell_and_neighbours)





quit()
# -------------------------------------------------------
# Index atoms in central cell
# -------------------------------------------------------

# Central cell index valid when translation vectors centred on zero
central_cell_index = int(0.5 * len(translations))
central_cell_atom_indices = np.arange(n_atoms_prim * central_cell_index,
                                      n_atoms_prim * (central_cell_index + 1))

check_central_cell = False
if check_central_cell:
    central_cell = []
    for ia in central_cell_atom_indices:
        central_cell.append(s_cell[ia])
    alex_xyz(output_dir + '/' + "aei_central_cell", central_cell)


# -----------------------------------------------------------------------
# Find neighbours to atoms of central cell that exist in adjacent cells
# -----------------------------------------------------------------------

def atoms_neighbouring_central_cell(central_cell_atom_indices, super_cell, radius, flatten=True):

    positions = [atom.position for atom in super_cell]
    d = spatial.distance_matrix(positions, positions)

    neighbours_of_central_cell = []
    for ia in central_cell_atom_indices:
        indices = np.where((d[ia, :] > 0.) & (d[ia, :] <= radius))[0]
        neighbours_of_central_cell.append(indices.tolist())

    if flatten:
        neighbours_of_central_cell = flatten_list(neighbours_of_central_cell)
        # Remove duplicates
        neighbours_of_central_cell = set(neighbours_of_central_cell)
        neighbours_of_central_cell = list(neighbours_of_central_cell)

    return neighbours_of_central_cell

def remove_central_cell_atoms(central_cell_atom_indices, neighbours_of_central_cell):

    neigbours_outside_cell = []
    for iN in neighbours_of_central_cell:
        if iN not in central_cell_atom_indices:
            neigbours_outside_cell.append(iN)

    return neigbours_outside_cell


neighbours_of_central_cell = atoms_neighbouring_central_cell(central_cell_atom_indices, s_cell, radius=1.8, flatten=True)
neighbours_outside_cell = remove_central_cell_atoms(central_cell_atom_indices, neighbours_of_central_cell)

neighbours = []
for iN in neighbours_outside_cell:
    position = s_cell[iN].position
    species = s_cell[iN].species
    neighbours.append(atoms.Atom(position=position, species=species))

unit_cell_and_neighbours = unit_cell + neighbours


# -----------------------------------------------------------------------
# Fold these atoms into central cell and see if they coincide with free-floating oxygens
# -----------------------------------------------------------------------

folded_neighbours = []

# Convert into fractional coordinates

# If so, then unit_cell_and_neighbours defines the set of atomic positions with which to work with

# Print out
#alex_xyz(output_dir + '/' + "aei_central_NN_cell", unit_cell_and_neighbours)









# -------------------------------------------------------
# Other ideas/notes
# -------------------------------------------------------
# Modify rings by doing one B-O-B
# Structurally relax the lattice maintaining bond angles
# OR
# Replace around the entirety of each ring before relaxing

# Once one has the supercell, could also randomly apply the snip (whilst retaining global structural connectivity)