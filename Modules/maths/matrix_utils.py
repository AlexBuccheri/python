import numpy as np 
#from typing import List

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


def find_centre(points: list):
    """
    Find the centre of a set of points
    
    Given a list of points, find the (bary) center of the boundary made up from 
    the points. point can correspond to 2D or 3D space 

    Parameters
    ----------
    points : array_like
             List of points, of the form points = [ [x1,y1,z1], [x2,y2,z2], ...]

    Returns
    -------
    s/sL : ndarray 
           Vector containing centre coordinates .shape(3)

    Notes
    -----
    Based on stack exchange:
    https://stackoverflow.com/questions/18305712/
    how-to-compute-the-center-of-a-polygon-in-2d-and-3d-space
  
    """
    n = points[0].shape[0] 
    if n < 2 or n > 3:
        raise Exception("Each point must be for a 2D or 3D coordinate")

    shapes = np.zeros(shape=(len(points)))
    for i,point in enumerate(points):
        shapes[i] = point.shape[0]
    if not np.all(shapes):
        raise Exception("All points must contain the same number of coordinates")
        
    s = np.zeros(shape=(3))
    sL = 0
    for i in range(len(points)):
        r0 = np.asarray(points[i - 1])
        r1 = np.asarray(points[i])
        r01 = np.linalg.norm(r1-r0)
        s += 0.5 * r01 * (r0 + r1)
        sL += r01
    return s/sL
