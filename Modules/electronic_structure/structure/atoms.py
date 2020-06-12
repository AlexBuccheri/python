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

class Atoms:
    def __init__(self):
        return

    def __new__(self, species: List[str], positions: List[float]) -> List[Atom]:
        molecule = []
        for ia in range(0, len(species)):
            molecule.append(Atom(species=species[ia], position=positions[ia]))
        return molecule

    def __add__(self, other: List[Atom]) -> List[Atom]:
        return self.molecule + other

    def __getitem__(self, atom_index: int):
        return self.molecule[atom_index]

    def __setitem__(self, atom_index: int, atom: Atom):
        self.molecule[atom_index] = atom

    def __delitem__(self, atom_index: int):
        del self.molecule[atom_index]



class CoordinateType(Enum):
    XYZ = 1
    FRACTIONAL = 2
