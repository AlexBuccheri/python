# Operations done on lattice vectors

import numpy as np
import math

from modules.maths.math import triple_product, angle_between_vectors


def parallelpiped_volume(lattice_vectors):
    """
    Volume of parallepiped cell
    :param lattice_vectors: lattice vectors, stored column-wise
    :return: volume
    """
    a = lattice_vectors[:, 0]
    b = lattice_vectors[:, 1]
    c = lattice_vectors[:, 2]
    return triple_product(a, b, c)


def reciprocal_lattice_vectors(a):
    """
    Get the reciprocal lattice vectors of real-space lattice vectors {a}
    :param a: lattice vectors, stored column-wise
    :return: Reciprocal lattice vectors, stored column-wise
    """
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


def simple_cubic_cell_translation_integers(lattice, cutoff: float):
    """
    Get the max translation integers required to cut a sphere of radius 'cutoff'
    from a cubic/cuboid supercell.
    :param lattice: lattice vectors, stored column-wise
    :param cutoff: radial cut-off
    :return: max integers
    """
    n_vectors = lattice.shape[1]
    return [int(np.ceil(cutoff / np.linalg.norm(lattice[:, i])))
            for i in range(0, n_vectors)]


def translation_integers_for_radial_cutoff(lattice, cutoff: float) :
    """
    Get the max translation integers required to cut a sphere of radius 'cutoff'
    from a non-cubic supercell.

    TODO(Alex) Give a mathematical description of how this works
    For cubic lattice vectors that form an orthogonal set, the distance from the
    origin to any face (equivdistant) defines the radius of the largest sphere that can
    fit inside. Hence, for a radial threshold centred at the origin of the box, one can define the
    required lengths of the box as

    For a nonorthogonal set of lattice vectors, one should also consider
    a x b =
    such that the integers used to define the super cell result in a cell
    The smallest integer that then faciliates the radial cutoff is defined as
    n_2 = max( cutoff * a x b / (a x b \cdot c), |c|), which is equivalent to
    n_2 = max( cutoff * a x b / V, |c|) where the triple product defines the volume V
    of a parallelpiped cell. One is therefore taking a ratio of volumes, which is dimensionally
    consistent.

    This is rather difficult to test. (See tests). One could verify visually

    This assumes the the sphere is centred at the origin of the cell and
    that the integers returned will be used like:

    for iz in range(-n2, n2+1):
      for iy in range(-n1, n1+1):
        for ix in range(-n0, n0+1):
          translation = np.matmul(lattice, [ix, iy, iz])

    i.e. that the max integers only correspond to the positive quadrant of the cell

    :param lattice: lattice vectors, stored column-wise
    :param cutoff: radial cut-off
    :return: max integers 
    """
    # lattice vectors
    a = lattice[:, 0]
    b = lattice[:, 1]
    c = lattice[:, 2]
    vol = triple_product(a, b, c)

    # May be no need to evaluate the max here but can't hurt
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
    Really, only corresponds to cubic cells.
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


def lattice_sum(lattice: np.ndarray, n, cutoff):
    """
    Sum over translation vectors, using a radial criterion
    :param lattice:
    :param n:
    :param cutoff:
    :return:
    """
    assert lattice.shape == (3,3), "lattice.shape /= (3,3)"
    assert len(n) == 3, "len(n) /= 3"
    cutoff_squared = cutoff * cutoff
    translations = []
    for k in range(-n[2], n[2]):
        for j in range(-n[1], n[1]):
            for i in range(-n[0], n[0]):
                r = np.matmul(lattice, np.array([i, j, k]))
                if np.dot(r, r) <= cutoff_squared:
                    translations.append(r)
    return translations
