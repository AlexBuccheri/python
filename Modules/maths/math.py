import numpy as np

# a. (b ^ c) = scalar
def triple_product(a, b, c):
    return np.dot(a, np.cross(b,c))

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between_vectors(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2',
        where clip ensures no failure for when  'v1' and 'v2' have
        a)the same direction and b) opposite directions

        Alternatively could just check this.
    """
    v1 = unit_vector(v1)
    v2 = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))

