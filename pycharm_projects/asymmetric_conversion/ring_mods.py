# -----------------------
# Read in ring and expand
# -----------------------

import numpy as np
from scipy import spatial
from enum import Enum

# My modules
from modules.fileio import read, write
from modules.electronic_structure.structure import atoms
from modules.maths import geometry


# -----------------------
# Local functions
# -----------------------

#TODO(Alex) Turn into function in module
# Should write something to control return type
def neighbour_list(positions, cutoff, include_self = False):
    n_atoms = len(positions)
    d = spatial.distance_matrix(positions, positions)
    neighbours = []
    if not include_self:
        for i in range(0, n_atoms):
            neighbours.append(np.where((d[i, :] > 0.) & (d[i, :] <= cutoff ))[0])
    elif include_self:
        for i in range(0, n_atoms):
            neighbours.append(np.where((d[i, :] >= 0.) & (d[i, :] <= cutoff))[0])
    return neighbours


def find_tetrahedral_units(species, nn_list):
    tetrahedra_units = []
    n_atoms = len(nn_list)
    for ia in range(0, n_atoms):
        if species[ia] == 'Si' and len(nn_list[ia]) == 4:
            unit = [ia] + nn_list[ia].tolist()
            assert(len(unit) == 5)
            tetrahedra_units.append(unit)

    return tetrahedra_units


def find_oxygens_in_ring(species, nn_list):
    oxy_ring_atoms = []
    n_atoms = len(nn_list)
    for ia in range(0, n_atoms):
        if species[ia] == 'O' and len(nn_list[ia]) == 2:
            oxy_ring_atoms.append(ia)
    return oxy_ring_atoms


def move_atoms_radially(atom_indices, positions, centre, length):
    new_positions = []
    for ia in atom_indices:
        position = positions[ia]
        # Point outwards from the centre
        radial_vector = np.asarray(position) - centre
        # scaled_radial_vector = length * unit vector
        scaled_radial_vector = length *  radial_vector / np.linalg.norm(radial_vector)
        assert (np.isclose(np.linalg.norm(scaled_radial_vector), length))
        new_positions.append(scaled_radial_vector + position)
    return new_positions


def find_closest_ring_oxygens(species, positions, si_oo_unit):

    class unit(Enum):
        Si = 0
        O1 = 1
        O2 = 2

    d = spatial.distance_matrix(positions, positions)
    si = si_oo_unit[unit.Si.value]
    assert (species[si] == 'Si')

    # After itself (Si) and two bonded oxy, get the closest two atoms
    neighbouring_ring_atom_distances = np.sort(d[si, :])[3:5]

    neighbouring_ring_atoms = []
    for distance in neighbouring_ring_atom_distances:
        neighbouring_ring_atoms.append(np.where(d[si, :] == distance)[0][0])

    assert (len(neighbouring_ring_atoms) == 2)
    assert (species[neighbouring_ring_atoms[0]] == 'O')
    assert (species[neighbouring_ring_atoms[1]] == 'O')
    return neighbouring_ring_atoms



# -------------------
# Main routine
# -------------------
species, positions = read.xyz("inputs/ring.xyz")
ring_centre = geometry.find_centre(positions)
n_atoms = len(species)
# In angstrom. Vesta shows bond-lengths differ slightly (:s)
rounded_bl = 1.7
nn_list = neighbour_list(positions, rounded_bl)
tetrahedra_units = find_tetrahedral_units(species, nn_list)
oxy_ring_atoms = find_oxygens_in_ring(species, nn_list)


# Remove ring oxygens from tetrahedral units. Could also look to do this with numpy
si_oo_units = []
for tetrahedron in tetrahedra_units:
    for oxy in oxy_ring_atoms:
        if oxy in tetrahedron:
            tetrahedron.remove(oxy)
    si_oo_units.append(tetrahedron)


# Move each si_oo_unit
length = 8
for unit in si_oo_units:
    # Find unit centre
    unit_centre = geometry.find_centre([positions[iatom] for iatom in unit])
    # Point outwards from the ring centre
    radial_vector = np.asarray(unit_centre) - ring_centre
    # scaled_radial_vector = length * unit vector
    scaled_radial_vector = length *  radial_vector / np.linalg.norm(radial_vector)
    assert (np.isclose(np.linalg.norm(scaled_radial_vector), length))
    # Update each position
    for atom in unit:
        positions[atom] = scaled_radial_vector + positions[atom]

# Move each oxy atom in the ring
oxy_ring_positions = move_atoms_radially(oxy_ring_atoms, positions, ring_centre, 5)
cnt = 0
for iatom in oxy_ring_atoms:
    positions[iatom] = oxy_ring_positions[cnt]
    cnt+=1

print_intermediate = False
if print_intermediate:
    molecule = atoms.Atoms(species, positions)
    write.xyz("select.xyz", molecule)



# For each Si-O-O, move a unit towards each of the two closest ring oxygens,
# until the Si-O bond_length matches the target, then remove one of the two oxygens from the unit
# Assume Br-O bond length

#TODO(Alex) Currently only for one si-o-o unit
si_oo_unit = si_oo_units[0]

neighbouring_ring_atoms = find_closest_ring_oxygens(species, positions, si_oo_unit)

# Move an si-o-o unit to bond to one of the ring oxys
bond_length_bo = 1.7
ring_oxy = neighbouring_ring_atoms[0]
pos_oxy = positions[ring_oxy]
print('pos_oxy', pos_oxy)
pos_si = positions[si_oo_unit[0]]
displacement = np.array(pos_si - pos_oxy)
scaled_displacement = bond_length_bo * displacement / np.linalg.norm(displacement)
assert(np.isclose(np.linalg.norm(scaled_displacement), bond_length_bo))

d_O1 = np.array(positions[si_oo_unit[1]] - pos_si)
d_O2 = np.array(positions[si_oo_unit[2]] - pos_si)

positions[si_oo_unit[0]] = pos_oxy + scaled_displacement
positions[si_oo_unit[1]] = positions[si_oo_unit[0]] + d_O1
positions[si_oo_unit[2]] = positions[si_oo_unit[0]] + d_O2

molecule = atoms.Atoms(species, positions)
write.xyz("one_unit.xyz", molecule)

quit()



new_units = []
bond_length_bo = 1.7
for iring,oxy in enumerate(neighbouring_ring_atoms):
    unit = []
    pos_oxy = positions[oxy]
    print('pos_oxy', pos_oxy)
    pos_si = positions[si_oo_units[0][0]]
    # Scaled so |displacement| = bond_length
    displacement = np.array(pos_oxy - pos_si)
    print("displacement ", displacement)
    scaled_displacement = bond_length_bo * displacement / np.linalg.norm(displacement)
    assert(np.isclose(np.linalg.norm(scaled_displacement), bond_length_bo))
    # Can't do this as need to keep the reference for the unit and shift again
    print("pos_oxy + displacement ",pos_oxy - scaled_displacement)
    unit.append(pos_oxy + scaled_displacement)
    # Can turn this into a list
    if iring == 0:
        unit_oxy = 1
    if iring == 1:
        unit_oxy = 2
    unit.append(unit[0]+ (positions[si_oo_units[0][unit_oxy]] - pos_si))
    new_units.append(unit)

for unit in new_units:
    print(unit)

# Remove old positions by leaving out
print(si_oo_units[0])

new_species = []
new_positions = []
for iatom in range(0, n_atoms):
    if iatom not in si_oo_units[0]:
        print(iatom)
        new_species.append(species[iatom])
        new_positions.append(positions[iatom])
assert(len(new_species) == len(new_positions))

# Add new ones in
for unit in new_units:
    for iatom,xyz in enumerate(unit):
        # Turn into a list
        if iatom == 0:
            new_species.append('Br')
        else:
            new_species.append('C')
        new_positions.append(xyz)

assert(len(new_species) == len(new_positions))

for i in range(0,len(new_species)):
    print(new_species[i], new_positions[i])

molecule = atoms.Atoms(new_species, new_positions)
write.xyz("one_unit.xyz", molecule)


# TODO(Alex)
# Deal with si-o-o and, potentially, floating oxygens
# Need to write something to identify the rings in the primitive cell




quit()
# Move each atom out along vector with ring centre
l = 5
new_positions = []
for ia in range(0, len(species)):
    position = positions[ia]
    # Point outwards from the centre
    radial_vector = np.asarray(position) - ring_centre
    unit_radial_vector = radial_vector / np.linalg.norm(radial_vector)
    scaled_radial_vector = l * unit_radial_vector
    assert( np.isclose(np.linalg.norm(scaled_radial_vector), l) )
    new_positions.append(scaled_radial_vector + position)

molecule = atoms.Atoms(species, new_positions)
write.xyz("test.xyz", molecule)