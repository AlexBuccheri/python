import numpy as np

from modules.fileio import read
from modules.electronic_structure.structure import atoms

names, positions = read.xyz("molecules/pyridine.xyz")
molecule = []
for iatom in range(0, len(names)):
    position = positions[iatom]
    molecule.append(atoms.Atom(names[iatom],position))
