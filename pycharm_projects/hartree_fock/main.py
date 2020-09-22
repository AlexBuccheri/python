import numpy as np

from modules.fileio import read
from modules.electronic_structure.structure import atoms

import read_basis

# Molecule: Positions in angstrom at this point
names, positions = read.xyz("molecules/pyridine.xyz")
molecule = []
for iatom in range(0, len(names)):
    position = positions[iatom]
    molecule.append(atoms.Atom(names[iatom],position))


# Basis sets: https://www.basissetexchange.org
read_basis.sto_3g("basis_sets/STO-3G/sto-3g.1-2.json")


# Integrals


# Symmetric orthogonalisation of the basis


# SCF