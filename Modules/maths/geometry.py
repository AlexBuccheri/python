import numpy as np

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
    n = len(points[0])
    if n < 2 or n > 3:
        raise Exception("Each point must be for a 2D or 3D coordinate")

    shapes = np.zeros(shape=(len(points)))
    for i,point in enumerate(points):
        shapes[i] = len(point)
    if not np.all(shapes):
        raise Exception("All points must contain the same number of coordinates")
        
    s = np.zeros(shape=(n))
    sum_r01 = 0
    for i in range(len(points)):
        r0 = np.asarray(points[i - 1])
        r1 = np.asarray(points[i])
        r01 = np.linalg.norm(r1-r0)
        s += 0.5 * r01 * (r0 + r1)
        sum_r01 += r01
    return s/sum_r01
