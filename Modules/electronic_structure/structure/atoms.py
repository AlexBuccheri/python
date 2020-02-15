# Atom class and functions acting on a list of Atom

import numpy as np
from enum import Enum

class Atom:
    def __init__(self, species, position):
        self.species = species
        #self.atomic_number = atomic_number
        self.position = np.asarray(position)

class CoordinateType(Enum):
    XYZ = 1
    FRACTIONAL = 2
