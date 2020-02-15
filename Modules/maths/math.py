import numpy as np

# a. (b ^ c) = scalar
def triple_product(a, b, c):
    return np.dot(a, np.cross(b,c))
