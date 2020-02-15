# Operations done on lattice vectors

from maths.math import triple_product
from maths.matrix_utils import off_diagonal_indices

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


def test_reciprocal_lattice_vectors(a,b):
    dp = np.empty(shape=(3,3))

    # a_i . b_j = 2pi * delta_ij
    for i in range(0,3):
        for j in range(0, 3):
            dp[i,j] = np.dot(a[:,i], b[:,j])

    off_diagonal_index = off_diagonal_indices(np.size(dp,0))
    off_diagonals = dp[off_diagonal_index]

    if np.any(dp.diagonal() != 2 * np.pi) or np.any(off_diagonals !=0):
        return False
    else:
        return True
