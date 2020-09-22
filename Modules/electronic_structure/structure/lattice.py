# Operations done on lattice vectors

import numpy as np
import math

from modules.maths.math import triple_product, angle_between_vectors

#Volume of parallepiped cell
#Expect lattice vectors stored column-wise
def parallelpiped_volume(lattice_vectors):
    a = lattice_vectors[:, 0]
    b = lattice_vectors[:, 1]
    c = lattice_vectors[:, 2]
    return triple_product(a, b, c)

# Reciprocal lattice vectors, stored column-wise
def reciprocal_lattice_vectors(a):
    b = np.zeros(shape=(3,3))
    b[:,0] = 2 * np.pi * np.cross(a[:,1], a[:,2]) / triple_product(a[:,0], a[:,1], a[:,2])
    b[:,1] = 2 * np.pi * np.cross(a[:,2], a[:,0]) / triple_product(a[:,1], a[:,2], a[:,0])
    b[:,2] = 2 * np.pi * np.cross(a[:,0], a[:,1]) / triple_product(a[:,2], a[:,0], a[:,1])
    return b


def cell_angles(lattice, return_unit='radian'):
    """ Angles of a parallelpiped unit cell
        See https://en.wikipedia.org/wiki/Crystal_structure#Lattice_systems for convention
    :param lattice: 3x3 np.array, lattice vectors stored columnwise
    :return: tuple containing cell angles
    """
    assert lattice.shape == (3,3)
    a = lattice[:, 0]
    b = lattice[:, 1]
    c = lattice[:, 2]
    alpha = angle_between_vectors(b, c)
    beta  = angle_between_vectors(c, a)
    gamma = angle_between_vectors(a, b)

    if return_unit == 'radian':
     return (alpha, beta, gamma)
    elif return_unit == 'degrees':
     return [angle * (180. / np.pi) for angle in (alpha, beta, gamma)]
    else:
        quit("Return unit not valid: ", return_unit)


def original_translation_integers_from_radial_cutoff(lattice, cutoff):

    translation_integers = np.array([0, 0, 0])
    for i in range(0, 2):
        for j in range(i+1, 3):
            angle = angle_between_vectors(lattice[:,i], lattice[:,j])
            t = cutoff / np.sin(angle)
            n_i = math.ceil(t/np.linalg.norm(lattice[:,i]))
            n_j = math.ceil(t / np.linalg.norm(lattice[:, j]))
            translation_integers[i] = max(n_i, translation_integers[i])
            translation_integers[j] = max(n_j, translation_integers[j])

    return translation_integers


# def translation_integers_from_radial_cutoff(lattice, cutoff) :
#     """
#     Should result in translation integers that correspond to summing
#     over a spherical volume (in the limit) or an isotropic sampling of the space,
#     for arbitrarily-shaped unit cells
#
#     :param lattice: lattice vectors, stored column-wise
#     :param cutoff: radial cut-off
#     :return:
#     """
#     # lattice vector magntiudes
#     a = np.linalg.norm(lattice[:, 0])
#     b = np.linalg.norm(lattice[:, 1])
#     c = np.linalg.norm(lattice[:, 2])
#     alpha, beta, gamma = cell_angles(lattice)
#
#     t_alpha = cutoff / np.sin(alpha)
#     t_beta  = cutoff / np.sin(beta)
#     t_gamma = cutoff / np.sin(gamma)
#
#     t_a = t_gamma
#     t_b = t_alpha
#     t_c = t_beta
#
#     return np.array([math.ceil(t_a / a), math.ceil(t_b / b), math.ceil(t_c / c)])


def simple_cubic_cell_translation_integers(lattice, cutoff: float):
    n_vectors = lattice.shape[1]
    return [int(np.ceil(cutoff / np.linalg.norm(lattice[:, i])))
            for i in range(0, n_vectors)]


def translation_integers_for_radial_cutoff(lattice, cutoff) :
    """
    Should result in translation integers that correspond to summing
    over a spherical volume (in the limit) or an isotropic sampling of the space,
    for arbitrarily-shaped unit cells

    :param lattice: lattice vectors, stored column-wise
    :param cutoff: radial cut-off
    :return:
    """
    # lattice vector magntiudes
    a = np.linalg.norm(lattice[:, 0])
    b = np.linalg.norm(lattice[:, 1])
    c = np.linalg.norm(lattice[:, 2])
    alpha, beta, gamma = cell_angles(lattice)

    t_a = max(cutoff / (a * np.sin(gamma)), cutoff / (a * np.sin(beta)))
    t_b = max(cutoff / (b * np.sin(alpha)), cutoff / (b * np.sin(gamma)))
    t_c = max(cutoff / (c * np.sin(beta)),  cutoff / (c * np.sin(alpha)))

    return np.array([math.ceil(t_a), math.ceil(t_b), math.ceil(t_c)])


def another_translation_integers_for_radial_cutoff(lattice, cutoff) :

    # lattice vectors
    a = lattice[:, 0]
    b = lattice[:, 1]
    c = lattice[:, 2]
    vol = triple_product(a, b, c)

    # n0 = cutoff * np.linalg.norm(np.cross(b, c)) / vol
    # n1 = cutoff * np.linalg.norm(np.cross(c, a)) / vol
    # n2 = cutoff * np.linalg.norm(np.cross(a, b)) / vol

    n0 = max(cutoff * np.linalg.norm(np.cross(b, c)) / vol, cutoff / np.linalg.norm(a))
    n1 = max(cutoff * np.linalg.norm(np.cross(c, a)) / vol, cutoff / np.linalg.norm(b))
    n2 = max(cutoff * np.linalg.norm(np.cross(a, b)) / vol, cutoff / np.linalg.norm(c))

    return np.array([math.ceil(n0), math.ceil(n1), math.ceil(n2)])


def get_extremal_integers(n, m, l):
    """
    Given 3 integers, return all permutations,
    including zeros
    Check this out for a potential generalisation:
    https://uk.mathworks.com/matlabcentral/answers/260300-how-can-i-generate-a-binary-matrix-of-permutations
    """
    extremal_integers = np.array([[ n,  0,  0],
                                  [-n,  0,  0],
                                  [ 0,  m,  0],
                                  [ 0, -m,  0],
                                  [ 0,  0,  l],
                                  [ 0,  0, -l],
                                  [ n,  m,  0],
                                  [-n,  m,  0],
                                  [ n, -m,  0],
                                  [-n, -m,  0],
                                  [ 0,  m,  l],
                                  [ 0, -m,  l],
                                  [ 0,  m, -l],
                                  [ 0, -m, -l],
                                  [ n,  0,  l],
                                  [-n,  0,  l],
                                  [ n,  0, -l],
                                  [-n,  0, -l],
                                  [ n,  m,  l],
                                  [-n, -m, -l],
                                  [-n,  m,  l],
                                  [ n, -m,  l],
                                  [ n,  m, -l],
                                  [ n, -m, -l],
                                  [-n,  m, -l],
                                  [-n, -m,  l]])
    return extremal_integers.transpose()

# Doesn't work. Norm of cross products are way too large
# def simplified_translation_integers_for_radial_cutoff(lattice, cutoff):
#
#     a = lattice[:, 0]
#     b = lattice[:, 1]
#     c = lattice[:, 2]
#
#     n0 = max(cutoff / np.sqrt(np.linalg.norm(np.cross(b, c))), cutoff / np.linalg.norm(a))
#     n1 = max(cutoff / np.sqrt(np.linalg.norm(np.cross(a, c))), cutoff / np.linalg.norm(b))
#     n2 = max(cutoff / np.sqrt(np.linalg.norm(np.cross(a, b))), cutoff / np.linalg.norm(c))
#
#     print(cutoff / np.sqrt(np.linalg.norm(np.cross(b, c))), cutoff / np.linalg.norm(a))
#
#     return [math.ceil(n) for n in [n0, n1, n2]]
#





