import numpy as np
from typing import List, Callable, Optional
from scipy.spatial import distance_matrix as scipy_distance_matrix

from ase.geometry import complete_cell, find_mic, wrap_positions
from ase.neighborlist import mic


def fixed_precision_rgkmax(atomic_number: int):
    """ Get rgkmax

    rgkmax which give a consistent total_energy precision, per elemental crystal.
    Found empirically by Sven working on the Delta-DFT project.
    :param atomic_number:
    :return:
    """
    fixed_rgkmax = np.array([
        5.835430, 8.226318, 8.450962, 8.307929,
        8.965808, 9.376204, 9.553568, 10.239864, 10.790975, 10.444355, 10.636286, 10.579793,
        10.214125, 10.605334, 10.356352, 9.932381, 10.218153, 10.466519, 10.877475, 10.774763,
        11.580691, 11.800971, 11.919804, 12.261896, 12.424606, 12.571031, 12.693836, 12.781331,
        12.619806, 12.749802, 12.681350, 12.802838, 12.785680, 12.898916, 12.400000, 10.596757,
        11.346060, 10.857573, 11.324413, 11.664200, 11.859519, 11.892673, 12.308470, 12.551024,
        12.740728, 12.879424, 13.027090, 13.080576, 13.230621, 13.450665, 13.495632, 13.261039,
        13.432654, 11.329591, 13.343047, 13.011835, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 15.134859, 14.955721,
        14.607311, 13.930505, 13.645267, 13.629439, 13.450805, 13.069046, 13.226699, 13.261342,
        13.365992, 13.557571, 13.565048, 13.579543, np.nan, 12.273924])
    return fixed_rgkmax[atomic_number + 1]


# def mufftin_radius_to_conserve_precision():
#     """ Given the minimum MT radius in the system, and the input RGKMAX,
#     return the MT radius to assure a consistent RGKMAX for element x.
#
#     :return:
#     """

def consistent_muffin_tin_radius(rgkmax_input: int, atomic_number_x: int, minimum_mt_radius: float):
    """ Given the smallest MT radius of species in a system, return the MT radius of element X

    Different `rgkmax` values are required for different species to obtain the same accuracy (in same quantity - I
    assume total energy). Because `rgkmax` is set once according to the smallest muffin tin radius, this determines
    the maximum $G$ vector for all species. As such, to obtain the desired `rgkmax` values for all other species in the
    calculation, their muffin tin radii must be set proportionally to min(MT), to obtain the target `rgkmax` values
    which ensure a consistent accuracy.

     \textrm{rgkmax}_{input} = \min(\textrm{MT}) * \textrm{G}_{max}

     therefore:

     \begin{equation}
     \textrm{G}_{max} = \textrm{rgkmax}_{input} / \min(MT).
     \end{equation}

     So for a given species, $X$, in the system:

     \begin{align}
     \textrm{rgkmax}_X &= MT_X * \textrm{G}_{max},  \\
     \textrm{rgkmax}_X &= MT_X * \textrm{rgkmax}_{input} / \min(MT).
     \end{align}

     Therefore, to obtain the desired $\textrm{rgkmax}$, $MT_X$ is the free parameter:

     \begin{equation}
     MT_X = \frac{\textrm{rgkmax}_X}{\textrm{rgkmax}_{input}} * min(MT)
     \end{equation}

    :return:
    """
    return (fixed_precision_rgkmax(atomic_number_x) / rgkmax_input) * minimum_mt_radius


class MoS2WS2Bilayer:
    """ Tabulated E1 System.

    Ref: https://github.com/nomad-coe/greenX-wp2/blob/main/Benchmarks/Heterobilayers/MoS2-WS2.xyz
    :return:
    """
    # lattice = np.array([[3.168394160510246, -3.3524146226426544e-10,  0.0],
    #                     [-1.5841970805453853, 2.7439098312114987, 0.0],
    #                     [6.365167633892244e-18 -3.748609667104484e-30, 39.58711265]])

    # Lattice with ~ zero terms rounded
    lattice = np.array([[ 3.168394160510246,  0.0,                0.0],
                        [-1.5841970805453853, 2.7439098312114987, 0.0],
                        [ 0.0,                0.0,                39.58711265]])

    positions = np.array([[0.00000000,  0.00000000,  16.68421565],
                          [1.58419708,  0.91463661,  18.25982194],
                          [1.58419708,  0.91463661,  15.10652203],
                          [1.58419708,  0.91463661,  22.90251866],
                          [0.00000000,  0.00000000,  24.46831689],
                          [0.00000000,  0.00000000,  21.33906353]])
    elements = ['W', 'S', 'S', 'Mo', 'S', 'S']
    atomic_numbers = [74, 16, 16, 42, 16, 16]

# TODO(Alex)
# This is more work than I would like:
# Need to ensure the cell origin is at the corner
# Would prefer fractional coordinates
# Need to test that the expression works on a simple 2D example I can work out the result to
# dr - np.rint(dr / cell_lengths) * cell_lengths

# def minimum_image_distance_matrix(a: np.ndarray, lattice: np.ndarray):
#
#     assert a.shape[1] == 3, "Expect vectors to be Euclidean"
#     n_vectors = a.shape[0]
#     cell_lengths = np.asarray([np.linalg.norm(lattice[i, :]) for i in range(0, 3)])
#
#     d = np.empty(shape=(n_vectors, n_vectors))
#     np.fill_diagonal(d, 0.)
#
#     # Upper triangle
#     for i in range(0, n_vectors):
#         for j in range(i + 1, n_vectors):
#             dr = a[j, :] - a[i, :]
#             print(dr - np.rint(dr / cell_lengths) * cell_lengths)
#             d[i, j] = np.linalg.norm(dr - np.rint(dr / cell_lengths) * cell_lengths)
#             d[j, i] = d[i, j]
#
#     return d
#
#
# def distance_matrix(a: np.ndarray):
#
#     assert a.shape[1] == 3, "Expect vectors to be Euclidean"
#     n_vectors = a.shape[0]
#     d = np.empty(shape=(n_vectors, n_vectors))
#     np.fill_diagonal(d, 0.)
#
#     # Upper triangle
#     for i in range(0, n_vectors):
#         for j in range(i + 1, n_vectors):
#             d[i, j] = np.linalg.norm(a[j, :] - a[i, :])
#             d[j, i] = d[i, j]
#
#     return d
#
# def optimal_muffin_tin_radii(positions: np.ndarray, lattice: np.ndarray, atomic_numbers: List[int]):
#
#     # Use rgkmax as proxy for MT radius (seem to be correlated)
#     species_min_mt = np.amin([fixed_precision_rgkmax(x) for x in atomic_numbers])
#
#     dm = scipy_distance_matrix(positions, positions)
#     dm2 = distance_matrix(positions)
#     dm3 = minimum_image_distance_matrix(positions - positions[0, :], lattice)
#     print(dm - dm2)
#     print(dm2)
#     print(dm3)
#




def minimum_muffin_tin_radius(positions: np.ndarray, lattice: np.ndarray, atomic_numbers: List[int]):
    """

    :param positions:
    :param lattice:
    :param atomic_numbers:
    :return:
    """
    bonds = find_minimum_bond_lengths(positions, lattice, atomic_numbers)
    an_min = atomic_number_species_with_mt_min(atomic_numbers)
    rgkmax_min = fixed_precision_rgkmax(an_min)

    for an, bond_length in bonds.items():
        denominator = 1. + (fixed_precision_rgkmax(an) / rgkmax_min)
        # Percentage of the bond that should be taken by the radius of an_min and an, respectively
        percentage_min = 1. / denominator
        percentage_x = (fixed_precision_rgkmax(an) / rgkmax_min) / denominator
        #print(percentage_min, percentage_x)
        print((an_min, an), percentage_min * bond_length)


def atomic_number_species_with_mt_min(atomic_numbers: List[int]):
    """ Species expected to have the smallest MT radus in the system.

    Use rgkmax as proxy for MT radius (seem to be correlated)

    :param atomic_numbers:
    :return:
    """
    index_rgkmax_min = np.argmin([fixed_precision_rgkmax(x) for x in atomic_numbers])
    an_min_mt = atomic_numbers[index_rgkmax_min]
    return an_min_mt


def find_minimum_bond_lengths(positions: np.ndarray, lattice: np.ndarray, atomic_numbers: List[int]) -> dict:
    """Find the minimum bond length between the species with (the anticipated) smallest
    MT radius, and each other species in a periodic cell.

    :return:
    """
    # Species with smallest MT radius
    an_min_mt = atomic_number_species_with_mt_min(atomic_numbers)

    # Get the atomic indices for each species
    unique_ans = set(atomic_numbers)
    indices = {}
    for an in unique_ans:
        ans = np.asarray(atomic_numbers)
        indices[an] = np.where(ans == an)[0]

    # Find all pair vectors in the unit cell, sorted according to species X_min_MT and Y
    minimum_bond_lengths = {}

    r_x = positions[indices[an_min_mt]]
    for an in unique_ans:
        r_y = positions[indices[an]]
        d_vectors = wrapped_displacement_vectors(r_x, r_y, lattice)
        minimum_bond_lengths[an] = minimum_bond_length(d_vectors)

    return minimum_bond_lengths


def wrapped_displacement_vectors(r_x: np.ndarray, r_y: np.ndarray, lattice: np.ndarray, remove_self_interaction=True) \
        -> np.ndarray:
    """ Find all displacement vectors between atom/s at position/s r_x
    and atom/s at position/s r_y, in the minimum image convention.

    :return:
    """
    displacement_vectors = []
    for i in range(r_x.shape[0]):
        for j in range(r_y.shape[0]):
            displacement_vectors.append(r_y[j, :] - r_x[i, :])

    if remove_self_interaction:
        zeros = (displacement_vectors == np.array([0., 0., 0.])).all(-1)
        non_zeros = [not x for x in zeros]
        displacement_vectors = np.asarray(displacement_vectors)
        displacement_vectors = displacement_vectors[non_zeros, :]

    # Apply minimum image convention to these vectors
    wrapped_vectors = mic(displacement_vectors, lattice, pbc=True)
    return wrapped_vectors


def minimum_bond_length(vectors: np.ndarray) -> float:
    """ Find the minimum bond length for a set of displacement vectors

    :param vectors:
    :return:
    """
    if vectors.shape[1] != 3:
        raise ValueError('Vectors should be stored (n_vectors, 3)')

    n_vectors = vectors.shape[0]
    norms = np.empty(shape=n_vectors)

    for i in range(n_vectors):
        norms[i] = np.linalg.norm(vectors[i, :])

    return np.min(norms)




# Method for finding the optimal MT radii
# Optimal is almost-touching MT spheres
# Want to find a pair list for each diatomic permutation in the system, then determine the optimal MT radii of
# non-touching spheres from this, using the ratio above
# This can be a starting point, however one may then need to reduce them
# if the charge overlaps too much

system = MoS2WS2Bilayer()
#print(find_minimum_bond_lengths(system.positions, system.lattice, system.atomic_numbers))
minimum_muffin_tin_radius(system.positions, system.lattice, system.atomic_numbers)

for x in [74, 42]:
    print(consistent_muffin_tin_radius(8, x, 1.05))
