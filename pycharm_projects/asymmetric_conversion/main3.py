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


# --------------------------------------------------------------
# For a given unit cell, list all surrounding atoms too
# such that every atom in the central cell is fully coordinated
# --------------------------------------------------------------

# Move central cell translation vector to start of list
translations.pop(int(0.5*len(translations)))
translations = [np.array([0, 0, 0])] + translations

# Unit cell that is fully coordinated
visualise = False
super_cell = supercell.build_supercell(unit_cell, translations)
assert len(super_cell) == np.product(n) * n_atoms_prim
if visualise:
    alex_xyz(output_dir + '/' + "aei_supercell", super_cell)

positions = [atom.position for atom in super_cell]
d_supercell = spatial.distance_matrix(positions, positions)
assert d_supercell.shape[0] == len(super_cell) and d_supercell.shape[1] == len(super_cell)
upper_bound_length = 1.8

# For atoms (ia) in central cell, check for neighbours outside of central cell
coordinating_atom_indices = []
for ia in range(0, n_atoms_prim):
    # This doesn't give the correct answer and I can't see why
    # indices = np.where((d_supercell[ia, n_atoms_prim:] > 0.) &
    #                    (d_supercell[ia, n_atoms_prim:] <= upper_bound_length))[0]
    indices = np.where((d_supercell[ia, :] > 0.) &
                       (d_supercell[ia, :] <= upper_bound_length))[0]
    # Put these straight into the coordinating cell
    coordinating_atom_indices.extend(indices)

# Remove dups
coordinating_atom_indices = list(set(coordinating_atom_indices))

# Remove indices associated with atoms in central cell
coordinating_atom_indices = np.asarray(coordinating_atom_indices)
coordinating_atom_indices = coordinating_atom_indices[coordinating_atom_indices >= n_atoms_prim]

# Store coordinating atoms
coordinating_atoms = []
for ia in coordinating_atom_indices:
    coordinating_atoms.append(super_cell[ia])

alex_xyz(output_dir + '/' + "aei_central_cell", unit_cell)
alex_xyz(output_dir + '/' + "aei_neighbours", coordinating_atoms)

quit()




# -------------------------------------------------------------------------------
# Translate each loose atom by all translation vectors and see if any
# equivalent position put each loose atom as a NN of an atom in the central cell
#
# Retain these and dump the uncoordinated atoms. This is the structure one works
# with for boron substitutions
#
# Once the boron framework has been constructed, one wraps all atoms outside
# the cell to inside
# -------------------------------------------------------------------------------

def index_loose_atoms(unit_cell_positions):
    """ Index loose or uncoordinated atoms """

    d = spatial.distance_matrix(unit_cell_positions, unit_cell_positions)
    upper_bound_length = 1.8
    loose_atom_indices = []
    for ia in range(0, n_atoms_prim):
        indices = np.where((d[ia, :] > 0.) & (d[ia, :] <= upper_bound_length))[0]
        if len(indices) == 0:
            loose_atom_indices.append(ia)

    assert len(loose_atom_indices) == 4
    return loose_atom_indices


def find_equivalent_position(loose_atom_index, unit_cell, translations, verbose=False):

    # Target bond length is ~ 1.529 - 1.6
    lower_bound = 1.4
    upper_bound = 1.8

    loose_atom_position = unit_cell[loose_atom_index].position
    loose_atom_positions = [loose_atom_position + translation for translation in translations]
    d_loose = spatial.distance_matrix(loose_atom_positions, unit_cell_positions)

    valid_equivalent = []
    for ia in range(0, len(loose_atom_positions)):
        # neighbours of loose atom in central cell
        indices = np.where((d_loose[ia, :] > lower_bound) & (d_loose[ia, :] <= upper_bound))[0]
        if len(indices) > 0:
            if verbose: print("distance:", d_loose[ia, indices])
            valid_equivalent.append(loose_atom_positions[ia].tolist())

    assert len(valid_equivalent) == 1
    return valid_equivalent[0]


def replace_loose_atoms(unit_cell, loose_atom_indices, equivalent_positions, visualise_parts=False):
    """ replace_loose_atoms """

    # Remove loose atoms from unit cell - can't use pop as it changes the indexing each time
    new_unit_cell = []
    removed_atoms = []
    for ia in range(0, len(unit_cell)):
        if ia not in loose_atom_indices:
            new_unit_cell.append(unit_cell[ia])
        else:
            removed_atoms.append(unit_cell[ia])

    # I know they're all oxygens but should treat this correctly if making general
    replacements = [atoms.Atom(position=position, species='O') for position in equivalent_positions]

    if visualise_parts:
        alex_xyz(output_dir + '/' + "aei_central_cell", unit_cell)
        alex_xyz(output_dir + '/' + "aei_replacements_cell", replacements)
        alex_xyz(output_dir + '/' + "aei_removed_atoms", removed_atoms)

    return new_unit_cell + replacements


unit_cell_positions = [atom.position for atom in unit_cell]
loose_atom_indices = index_loose_atoms(unit_cell_positions)

# One equivalent position per loose atom
equivalent_positions = []
for ia in loose_atom_indices:
    equivalent_positions.append(find_equivalent_position(ia, unit_cell, translations))

unit_cell = replace_loose_atoms(unit_cell, loose_atom_indices, equivalent_positions)
assert len(unit_cell) == n_atoms_prim
alex_xyz(output_dir + '/' + "aei_primitive_cell_fullcon", unit_cell)
print("Output a cell with full connectivity. "
      "Atoms that lie outside the primitive cell can be folded back in")

# -------------------------------------------------------
# Test folding back in
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