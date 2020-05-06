import numpy as np 

def off_diagonal_indices(n: int):
    """
    Get indices of all off-diagonal elements of a square matrix, with dimensions n. 

    Parameters
    ----------
    n : int
        Dimensions of square matrix

    Returns
    -------
    c : int array 
        Rototation matrix .shape(n,n)

    Notes
    -----
    Based on stack exchange example:
    https://stackoverflow.com/questions/35746806/
    how-to-get-indices-of-non-diagonal-elements-of-a-numpy-array
  
    """   
    return np.where(~np.eye(n, dtype=bool))


# TODO(Alex) Add asserts for when this won't work
def rotation_to_align_a_with_b(a, b):
    """
    Generates a rotation matrix to align unit vector 'a' with unit vector 'b', 
    such that 'a' and 'b' point in the same direction. 

    Parameters
    ----------
    a : array_like
        First vector.
    b : array_like
        Second vector.

    Returns
    -------
    c : ndarray 
        Rototation matrix .shape(3,3)

    Notes
    -----
    Based on stack exchange:
    https://math.stackexchange.com/questions/180418/
    calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
  
    """
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if not np.allclose(a, a/norm_a):
        print('Input a vector not unit normal - normalising')
        a = a / norm_a
        print(a)
    if not np.allclose(b, b/norm_b):
        print('Input b vector not unit normal - normalising')
        b = b / norm_b
        print(b)

    v = np.cross(a,b)
    #s = np.linalg.norm(v)
    c = np.dot(a,b)
    f = 1./(1. + c)
    vmat = np.array([[    0, -v[2],  v[1]],
                     [ v[2],     0, -v[0]],
                     [-v[1],  v[0],     0]])
    return np.eye(3,3) + vmat + f *(np.matmul(vmat,vmat))
