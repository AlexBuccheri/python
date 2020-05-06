# Atom class and functions acting on a list of Atom

import numpy as np
from enum import Enum
from typing import List

from modules.parameters import elements

class Atom:
    def __init__(self, species: str, position: List[float]):
        self.species = species
        self.position = np.asarray(position)
    def atomic_number(self) -> int:
        return elements.symbol_to_an[self.species]

# Should this be a class or a function?
class Atoms:
    def __init__(self):
        return
    def __new__(self, species: List[str], positions: List[float]) -> List[Atom]:
        molecule = []
        for ia in range(0, len(species)):
            molecule.append(Atom(species=species[ia], position=positions[ia]))
        return molecule

class CoordinateType(Enum):
    XYZ = 1
    FRACTIONAL = 2
