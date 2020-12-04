"""
One basis position per cell, at (0,0,0)
Find nearest neighbour distance and number of NN to central atom
"""
import numpy as np
from scipy import spatial

from reference_systems import System

from modules.electronic_structure.structure import lattice
from modules.electronic_structure.structure import bravais
from modules.electronic_structure.structure import atoms
from modules.electronic_structure.structure import supercell
from modules.fileio import write


def get_positions(unit_cell):
    n_atoms = len(unit_cell)
    # Require this shape for default distance matrix function input
    positions = np.empty(shape=(n_atoms, 3))
    for i, atom in enumerate(unit_cell):
        positions[i, :] = atom.position[:]
    return positions


def get_neighbour_details(system: System):
    unit_cell = system.unit_cell
    lattice_vectors = system.lattice

    coordinating_translations = supercell.translation_vectors(
        lattice_vectors, [3, 3, 3], centred_on_zero=True)

    central_index = int(0.5 * len(coordinating_translations))
    assert np.all(coordinating_translations[central_index] == 0.), "translation = [0, 0, 0]"

    super_cell = supercell.build_supercell(unit_cell, coordinating_translations)
    positions = get_positions(super_cell)
    d_matrix = spatial.distance_matrix(positions, positions)
    #write.xyz("cell", super_cell)

    NN_dist = np.sort(d_matrix[central_index, :])[1]
    tol = 1.e-5
    NN_indices = np.where((d_matrix[central_index, :] > 0.) &
                          (d_matrix[central_index, :] <= NN_dist + tol))[0]

    #write.xyz("reduced_cell", [super_cell[ia] for ia in NN_indices])
    return NN_dist, len(NN_indices)


def get_surface_atoms(super_cell, radius, n_neighbours_bulk):
    positions = get_positions(super_cell)
    d_matrix = spatial.distance_matrix(positions, positions)

    surface_atom_indices = []
    for ia in range(0, d_matrix.shape[0]):
        indices = np.where((d_matrix[ia, :] > 0.) &
                           (d_matrix[ia, :] <= radius))[0]
        if len(indices) < n_neighbours_bulk:
            surface_atom_indices.append(ia)

    return surface_atom_indices


def get_radial_distance_mean_and_sd(full_cell, surface_atom_indices):
    """
    Get radial distances from the origin and use to
    compute mean and SD
    """
    distances = np.empty(shape=len(surface_atom_indices))
    for i, ia in enumerate(surface_atom_indices):
        distances[i] = np.linalg.norm(full_cell[ia].position)
    return distances


# Dummy FCC
def dummy_fcc():
    system_label = 'dummy_fcc'
    al = 10
    unit_cell = atoms.Atoms(['Si'], [[0, 0, 0]])
    lattice_vectors = bravais.face_centred_cubic(al)
    return System(system_label, unit_cell, lattice_vectors, al)


def test_fcc():
    print("Testing FCC")

    system = dummy_fcc()
    lattice_vectors = system.lattice
    unit_cell = system.unit_cell
    al = system.al

    cutoff = 20
    max_integers_cubic = lattice.simple_cubic_cell_translation_integers(lattice_vectors, cutoff)
    max_integers_gen = lattice.translation_integers_for_radial_cutoff(lattice_vectors, cutoff)

    # "For FCC, the simple cubic max integers should differ from those found "
    # "using the generic expression"
    assert max_integers_cubic[0] == 3
    assert max_integers_cubic[1] == 3
    assert max_integers_cubic[2] == 3

    assert max_integers_gen[0] == 4
    assert max_integers_gen[1] == 4
    assert max_integers_gen[2] == 4

    NN_dist, num_NN = get_neighbour_details(dummy_fcc())

    translations_cubic = lattice.lattice_sum(lattice_vectors, max_integers_cubic, cutoff)
    super_cell_cubic = supercell.build_supercell(unit_cell, translations_cubic)

    surface_atom_indices = get_surface_atoms(super_cell_cubic, NN_dist, num_NN)
    distances = get_radial_distance_mean_and_sd(super_cell_cubic, surface_atom_indices)
    mean_cubic = np.mean(distances)
    sd_cubic   = np.std(distances)
    max_cubic  = np.max(distances)

    translations_gen = lattice.lattice_sum(lattice_vectors, max_integers_gen, cutoff)
    super_cell_gen = supercell.build_supercell(unit_cell, translations_gen)

    print(len(super_cell_cubic), len(super_cell_gen))


    surface_atom_indices = get_surface_atoms(super_cell_gen, NN_dist, num_NN)
    distances = get_radial_distance_mean_and_sd(super_cell_gen, surface_atom_indices)
    mean_genr = np.mean(distances)
    sd_genr = np.std(distances)
    max_genr= np.max(distances)

    # I think this is the best one can do: Evaluate the standard deviatiosn
    print(cutoff)
    print(mean_cubic, sd_cubic, max_cubic)
    print(mean_genr, sd_genr, max_genr)

    assert mean_cubic < mean_genr
    assert sd_cubic > sd_genr



test_fcc()