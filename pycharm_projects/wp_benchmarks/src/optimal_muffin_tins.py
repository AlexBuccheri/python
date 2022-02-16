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
    # return (fixed_precision_rgkmax(atomic_number_x) / rgkmax_input) * minimum_mt_radius


class MoS2WS2Bilayer:
    """ Tabulated E1 System.

    Ref: https://github.com/nomad-coe/greenX-wp2/blob/main/Benchmarks/Heterobilayers/MoS2-WS2.xyz
    :return:
    """
    # lattice = np.array([[3.168394160510246, -3.3524146226426544e-10,  0.0],
    #                     [-1.5841970805453853, 2.7439098312114987, 0.0],
    #                     [6.365167633892244e-18 -3.748609667104484e-30, 39.58711265]])

    # Lattice with ~ zero terms rounded
    lattice = np.array([[3.168394160510246, 0.0, 0.0],
                        [-1.5841970805453853, 2.7439098312114987, 0.0],
                        [0.0, 0.0, 39.58711265]])

    positions = np.array([[0.00000000, 0.00000000, 16.68421565],
                          [1.58419708, 0.91463661, 18.25982194],
                          [1.58419708, 0.91463661, 15.10652203],
                          [1.58419708, 0.91463661, 22.90251866],
                          [0.00000000, 0.00000000, 24.46831689],
                          [0.00000000, 0.00000000, 21.33906353]])
    elements = ['W', 'S', 'S', 'Mo', 'S', 'S']
    atomic_numbers = [74, 16, 16, 42, 16, 16]


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
        # print(percentage_min, percentage_x)
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


def optimal_smallest_muffin_tin_radius(an_x, an_y, bond_length):
    """For each bond length, compute the smaller MT radius associated with the two
    elements comprising the bond, according to the ratio of rgkmaxs.

    MT_y = (rgkmax_y / rgkmax_x) * MT_x
    where MT_x corresponds to the smallest MT in the system, and rgkmax_i are a consistent set of tabulated
    rgkmax values (leading to the same precision in total energy), the smallest MT radius can be used
    to determine consistent MT radii for all other elements in the system.

    The sum of two muffin tin radii cannot exceed the bond length:
    MT_y + MY_x = bond_length

    therefore substituting for MT_y and rearranging for MT_x:
    MT_x = bond_length / (1 + (rgkmax_y / rgkmax_x))

    :param an_x:
    :param an_y:
    :param bond_length:
    :return:
    """
    mt_x = bond_length / (1 + fixed_precision_rgkmax(an_y) / fixed_precision_rgkmax(an_x))
    return mt_x


if __name__ == "__main__":
    # Get touching spheres, using optimal radius ratios.
    # In the case of the bond length = min(X_min_mt - Y), the MT spheres should touch.
    # In all other cases, the Y_mt will be determined by X_min_mt
    # Basically want to get the touching sphere recommendations, select the smallest MT radius for the element
    # expected to have the smallest MT radius, else I cannot maintain valid ratios with the other species without overlap
    # Then I just reduce all MT radii by 10-30% (printing out the radii) and do the charge plotting to view an optimal

    system = MoS2WS2Bilayer()
    an_mt_min = atomic_number_species_with_mt_min(system.atomic_numbers)
    min_bond_lengths = find_minimum_bond_lengths(system.positions, system.lattice, system.atomic_numbers)

    print(f"Minimum bond lengths between atomic number {an_mt_min} and:")
    for an, min_b_len in min_bond_lengths.items():
        print(f' Atomic number {an}, min bond length {min_b_len}')

    print("For each bond, compute the MT radius associated with the element"
          "with the smallest MT radius. \n Do this for each bond, and choose the "
          "minimum.")
    mt_mins = []
    for an_y, bond_length in min_bond_lengths.items():
        radius = optimal_smallest_muffin_tin_radius(an_mt_min, an_y, bond_length)
        mt_mins.append(round(radius, ndigits=4))

    print(mt_mins)
    # TODO Run this for a few scaling factors, put into the code, check for when the density overlap becomes small
    scaling_factor = 0.9
    mt_min = scaling_factor * min(mt_mins)
    print(f"Smallest muffin tin radius for element {an_mt_min} is {mt_min}")

    print("Use this to determine the MT radii for all other elements in the system")
    for an_y in set(system.atomic_numbers):
        mt_y = fixed_precision_rgkmax(an_y) / fixed_precision_rgkmax(an_mt_min) * mt_min
        print(f"Atomic Number {an_y}, MT radius {mt_y}, sum of MT radii {mt_min + mt_y}, min_bond_length {min_bond_lengths[an_y]}")






