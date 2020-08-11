# -----------------------
# Read in ring and expand
# -----------------------

import numpy as np
from scipy import spatial
from enum import Enum
from statistics import mean

# My modules
from modules.fileio import read, write
from modules.electronic_structure.structure import atoms
from modules.maths import geometry

import cell_operations

# -----------------------
# Local functions
# -----------------------

# Global to remove
# Assume Br-O bond length
bond_length_bo = 1.7

def neighbour_list(positions, cutoff, include_self=False)->list:
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

def duplicates_exist(my_list: list):
    return not (len(my_list) == len(set(my_list)))

def remove_sharing_oxygens(si_oo_units):
    """
    In the two tetrahedra that don't belong to the ring, remove the
    oxygens that are shared with the tetrahedra in the ring

    :return:
    """
    new_si_oo_units = []
    all_units = flatten(si_oo_units)
    for unit in si_oo_units:
        if len(unit) == 5:
            for i in unit:
                if all_units.count(i) > 1:
                    unit.remove(i)
            new_si_oo_units.append(unit)
        else:
            new_si_oo_units.append(unit)

    assert not duplicates_exist(flatten(new_si_oo_units)), \
        "Should be no reoccuring atomic indices in si-o-o units"

    return new_si_oo_units


def in_plane(posA, posB, tol, plane_index) -> bool:
    delta_height = (posA[plane_index] - posB[plane_index]) / posA[plane_index]
    return abs(delta_height) <= tol

def find_oxygens_in_ring(species, positions, nn_list:list, height_tol, plane_index):
    assert len(species) == len(positions)
    assert len(nn_list) == len(species)

    oxy_ring_atoms = []
    n_atoms = len(nn_list)
    for ia in range(0, n_atoms):
        speciesNN = [species[i] for i in nn_list[ia]]
        if species[ia] == 'O' and speciesNN.count('Si') == 2:
            oxy_ring_atoms.append(ia)

    # Make sure all found atoms are in the same plane as their Si neighbours
    refined_oxy_ring_atoms = []
    for ia in oxy_ring_atoms:
        si_neighbours = [i for i in nn_list[ia] if species[i] == 'Si']
        if in_plane(positions[ia], positions[si_neighbours[0]], height_tol, plane_index=0) and \
                in_plane(positions[ia], positions[si_neighbours[1]], height_tol, plane_index=0):
            refined_oxy_ring_atoms.append(ia)

    return refined_oxy_ring_atoms


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

# Scales a vector to have a magnitude (length) equal to bond_length
def scaled_vector(bond_length, v):
    # scaled_radial_vector = length * unit vector
    scaled_v = bond_length * v / np.linalg.norm(v)
    assert (np.isclose(np.linalg.norm(scaled_v), bond_length))
    return scaled_v

# Translate si-o-o to each of the two closest ring oxygens, convert Si to Br
# and remove alternating oxygens
def translate_si_o_o_unit(species, positions, si_oo_unit, neighbouring_ring_atoms):

    assert (len(neighbouring_ring_atoms) == 2)
    si_index = 0
    o1_index = 1
    o2_index = 2
    assert (species[si_oo_unit[si_index]] == 'Si')
    assert (species[si_oo_unit[o1_index]] == 'O')
    assert (species[si_oo_unit[o2_index]] == 'O')

    pos_si = positions[si_oo_unit[si_index]]
    pos_o1 = positions[si_oo_unit[o1_index]]
    pos_o2 = positions[si_oo_unit[o2_index]]

    d_O1 = np.array(pos_o1 - pos_si)
    d_O2 = np.array(pos_o2 - pos_si)
    d_O = [d_O1, d_O2]

    translated_positions = []
    translated_species = []

    for i,ring_oxy in enumerate(neighbouring_ring_atoms):
        pos_oxy = positions[ring_oxy]
        scaled_displacement = scaled_vector(bond_length_bo, np.array(pos_si - pos_oxy))

        translated_positions.append(pos_oxy + scaled_displacement)
        translated_positions.append(pos_oxy + scaled_displacement + d_O[i])
        translated_species += ['B', 'O']

    return translated_species, translated_positions


# Ref: http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
def flatten(l, ltypes=(list, tuple)):
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

def moleculr_from_indices(species, positions, atomic_indices):
    new_species = []
    new_positions = []
    for ia in atomic_indices:
        new_species.append(species[ia])
        new_positions.append(positions[ia])
    return atoms.Atoms(new_species, new_positions)



def identify_si_o_o_units(species, positions, plane_index) -> list:
    """
    :brief Identify silicons and neighbouring oxygens that are not in the ring

    :param species: list of species
    :param positions: list of positions
    :return: list contain indices of atoms in an si_oo_unit
    """

    assert len(species) == 40, "Expect 40 atoms in both rings"
    assert len(positions) == 40, "Expect 40 atoms in both rings"

    # In angstrom
    # Coordination really doesn't matter when everything's fully-coordinated
    rounded_bl = 2.2  #  prior value 1.7
    nn_list = neighbour_list(positions, rounded_bl)

    tetrahedra_units = find_tetrahedral_units(species, nn_list)
    assert len(tetrahedra_units) == 10, "10 Si, hence 10 tetradrons in the rings"

    # Oxygens in the ring are ~ going to be in the same plane as the Si
    # Assumes x is perpendicular to the plane
    height_tol = 0.075
    oxy_ring_atoms = find_oxygens_in_ring(species, positions, nn_list, height_tol, plane_index)
    assert len(oxy_ring_atoms) == 8, "Did not find 8 oxygen atoms in the ring"

    # Remove ring oxygens from tetrahedral units.
    si_oo_units = []
    for tetrahedron in tetrahedra_units:
        unit = [atom for atom in tetrahedron if atom not in oxy_ring_atoms]
        si_oo_units.append(unit)

    # Should be 8 units with 3 atoms, and two units with 5
    # (those that are not in the plane)
    lengths = [len(unit) for unit in si_oo_units]
    assert lengths.count(3) == 8, "Expected 8 units to have 3 atoms"
    assert lengths.count(5) == 2, "Expected 8 units to have 3 atoms"

    return si_oo_units, oxy_ring_atoms


#TODO(Alex) Routine is too large - split it
def convert_ring(species, positions):


    n_atoms = len(species)
    ring_centre = geometry.find_centre(positions)
    si_oo_units, oxy_ring_atoms = identify_si_o_o_units(species, positions, plane_index=0)
    # Else same atom will get shifted for each time it reappears in a unit
    si_oo_units = remove_sharing_oxygens(si_oo_units)

    # Move each si_oo_unit, length found from trial and error (specific to AEI)
    length = 8
    for unit in si_oo_units:
        # Find unit centre
        unit_centre = geometry.find_centre([positions[iatom] for iatom in unit])
        # Point outwards from the ring centre
        radial_vector = np.asarray(unit_centre) - ring_centre
        scaled_radial_vector = scaled_vector(length, radial_vector)
        # Update each position
        for atom in unit:
            positions[atom] = scaled_radial_vector + positions[atom]

    # Move each oxy atom in the ring
    oxy_ring_positions = move_atoms_radially(oxy_ring_atoms, positions, ring_centre, 2)  # 5
    cnt = 0
    for iatom in oxy_ring_atoms:
        positions[iatom] = oxy_ring_positions[cnt]
        cnt += 1

    print_intermediate = True
    if print_intermediate:
        molecule = atoms.Atoms(species, positions)
        write.xyz("select.xyz", molecule)

    quit('up to here - looks like things move unexpectedly')


    # Assume Br-O bond length
    bond_length_bo = 1.7

    # For each Si-O-O of a former tetrahedron, move a unit towards each of the two closest ring oxygens
    # cleave off one of the oxy and convert Si -> Br
    translated_species = []
    translated_positions = []

    for si_oo_unit in si_oo_units:
        neighbouring_ring_atoms = find_closest_ring_oxygens(species, positions, si_oo_unit)
        ts, tp = translate_si_o_o_unit(species, positions, si_oo_unit, neighbouring_ring_atoms)
        translated_species += ts
        translated_positions += tp

    # Remove old atoms
    atom_indices = np.delete(np.arange(0, n_atoms), flatten(si_oo_units))
    new_species = []
    new_positions = []
    for iatom in atom_indices:
        new_species.append(species[iatom])
        new_positions.append(positions[iatom])

    # Add new atoms
    assert len(translated_species) % 2 == 0
    for iatom in range(0, len(translated_species)):
        new_species.append(translated_species[iatom])
        new_positions.append(translated_positions[iatom])

    assert (len(new_species) == len(new_positions))

    # Smart way to do this is to create sets of pair indices from the translated silicons,
    # then put oxygens between them
    # I always do operations in pairs, so should go [b,o,b,o,b,o...] => elements 0 and 2 are a pair in the ring.
    n_atoms = len(new_species)
    boron_indices = [i for i in range(0, n_atoms) if new_species[i] == 'B']

    # Create list of pairs
    boron_pairs = []
    for i in range(0, len(boron_indices), 2):
        boron_A = boron_indices[i]
        boron_B = boron_indices[i + 1]
        print(new_species[boron_A], new_species[boron_B])
        boron_pairs.append([boron_A, boron_B])

    # Put oxygens between these pairs
    for pair in boron_pairs:
        pos_A = np.asarray(new_positions[pair[0]])
        pos_B = np.asarray(new_positions[pair[1]])
        print(np.linalg.norm(pos_A - pos_B))
        pos_oxy = 0.5 * (pos_A + pos_B)
        new_positions.append(pos_oxy.tolist())
        new_species.append('O')

    return new_species, new_positions


# ------------------------------------
# For operating on the unit cell
# ------------------------------------

def get_ring(unit_cell, in_correct_region, erroneous_indices) ->list:
    """
    :brief Get upper or lower ring
    (determined by the choice of in_correct_region and erroneous_indices)

    :return ring_unit: atoms in ring
    """
    ring_unit = []
    for ia, atom in enumerate(unit_cell):
        if in_correct_region(atom) and (ia not in erroneous_indices):
            ring_unit.append(atom)
    return ring_unit


def get_connecting_chains(unit_cell) -> dict:
    """
    :brief Get connecting chains. Hard-coded as just easier

    Doesn't matter if the central connecting chains take some duplicate atoms
    already found in the rings: Can either label as ghosts or remove duplicates when
    recombining
    Note, the tabulated ghost atoms are silicons in the rings, but the point of connections
    are actually oxygens => these ghost indices are useless

    :return: atoms in connecting chains
    """
    # Note, Vesta indexing starts at 1, so subtract 1 in return statement
    chain1_indices = {'main': [13, 10, 29, 65, 52, 9, 12, 28, 64], 'ghost': [48, 49]}
    chain2_indices = {'main': [68, 32, 2, 54, 4, 33, 69, 1, 5], 'ghost': [44, 45 ]}

    return {'main': [unit_cell[ia-1] for ia in chain1_indices['main']],
            'ghost': [unit_cell[ia-1] for ia in chain1_indices['ghost']]}, \
           {'main': [unit_cell[ia-1] for ia in chain2_indices['main']],
            'ghost': [unit_cell[ia-1] for ia in chain2_indices['ghost']]}


def ring_substitutions(aei_ring_unit):
    """
    :brief

    :result
    """
    assert isinstance(aei_ring_unit, list), \
        "aei_ring_unit should be list[atoms.Atom]"
    assert isinstance(aei_ring_unit[0], atoms.Atom), \
        "aei_ring_unit should be list[atoms.Atom]"
    species = [atom.species for atom in aei_ring_unit]
    positions = [atom.position for atom in aei_ring_unit]
    new_species, new_positions = convert_ring(species, positions)
    return atoms.Atoms(new_species, new_positions)




# # Simple hack as I know the structure beforehand
# # If in the correct spatial region and is bonded (could remove the 2nd assertion)
# def ring_indices(unit_cell, in_ring):
#
#     positions = [atom.position for atom in unit_cell]
#     d = spatial.distance_matrix(positions, positions)
#     max_bond_length = 1.7
#
#     indices = []
#     for ia in range(0, len(unit_cell)):
#         if in_ring(positions[ia]):
#             neighbours = np.where((d[ia, :] > 0.) & (d[ia, :] <= max_bond_length))[0]
#             if len(neighbours) > 0:
#                 indices.append(ia)
#     return indices
#
#
# def get_ring(ring_indices, unit_cell):
#     species = []
#     positions = []
#     for i in ring_indices:
#         species.append(unit_cell[i].species)
#         positions.append(unit_cell[i].position)
#     return species, positions
#
#
# def swap_ring(ring_indices, unit_cell):
#     species, positions = get_ring(ring_indices, unit_cell)
#     # write.xyz("extracted_ring.xyz", atoms.Atoms(species, positions))
#     new_species, new_positions = convert_ring(species, positions)
#     converted_ring = atoms.Atoms(new_species, new_positions)
#
#     # New structure, with old ring swapped out for the new
#     unit_cell_minus_ring = []
#     for i, atom in enumerate(unit_cell):
#         if i not in ring_indices:
#             unit_cell_minus_ring.append(atom)
#
#     return unit_cell_minus_ring + converted_ring
#
#




# -------------------
# Main routine
# -------------------

# Read in CIF and convert to primitive cell
# directories = cell_operations.Directories(structure='aei',
#                                           input='inputs',
#                                           output='aei_outputs')
#
# unit_cell, lattice_vectors = cell_operations.get_primitive_unit_cell(directories)
# translations = cell_operations.translations_for_fully_coordinated_unit(unit_cell, lattice_vectors)
# coordinating_atoms = cell_operations.find_atoms_neighbouring_central_cell(
#     unit_cell, translations)
#
# # Breaks the ring routine if I add coordinating_atoms here
# #unit_cell = unit_cell + coordinating_atoms
#
# # For AEI primitive cell, the rings are:
# # Every connected atom above x value of 11.52979
# # Every connected atom below an x value of 7
# top_ring_indices = ring_indices(unit_cell, lambda pos: pos[0] > 11.5)
#
# structure = swap_ring(top_ring_indices, unit_cell)
# bottom_ring_indices = ring_indices(structure, lambda pos: pos[0] <= 7)
# structure = swap_ring(bottom_ring_indices, structure)
#
# write.xyz("new_ring.xyz", structure)




# def ring_test():
#     species, positions = read.xyz("inputs/ring.xyz")
#     new_species, new_positions = convert_ring(species, positions)
#     molecule = atoms.Atoms(new_species, new_positions)
#     write.xyz("one_unit.xyz", molecule)
#


# quit()
#
# first_boron_index = next(i for i,symbol in enumerate(new_species) if symbol == 'B')
#
# #TODO(Alex) This doesn't work properly
# # Find 2nd-shortest B-B bond length.
# n_atoms = len(new_species)
# boron_indices = [i for i in range(0, n_atoms) if new_species[i] == 'B']
# new_positions = np.asarray(new_positions)
# new_species = np.asarray(new_species)
# boron_positions = new_positions[boron_indices]
# d_boron = spatial.distance_matrix(boron_positions, boron_positions)
# molecule = atoms.Atoms(new_species[boron_indices], boron_positions)
# write.xyz("borons.xyz", molecule)
#
# # This gets lucky
# second_shortest_distance = np.sort(d_boron[0,:])[2]
#
# # Create list of pairs
# boron_pairs = []
# for i in range(0, len(boron_indices)):
#     j = np.where(d_boron[i,:] == second_shortest_distance)
#     print(j)
#     boron_pairs.append([i,j])
#
# # Put oxygens between these boron pairs
#
#
# molecule = atoms.Atoms(new_species, new_positions)
# write.xyz("one_unit.xyz", molecule)
#
#
#
# # TODO(Alex)
# # Deal with si-o-o and, potentially
# # Need to write something to identify the rings in the primitive cell: Just list manually.
#
#
#
#
# quit()
# # Move each atom out along vector with ring centre
# l = 5
# new_positions = []
# for ia in range(0, len(species)):
#     position = positions[ia]
#     # Point outwards from the centre
#     radial_vector = np.asarray(position) - ring_centre
#     unit_radial_vector = radial_vector / np.linalg.norm(radial_vector)
#     scaled_radial_vector = l * unit_radial_vector
#     assert( np.isclose(np.linalg.norm(scaled_radial_vector), l) )
#     new_positions.append(scaled_radial_vector + position)
#
# molecule = atoms.Atoms(species, new_positions)
# write.xyz("test.xyz", molecule)