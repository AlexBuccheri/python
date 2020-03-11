# Build entos input files for bulk band structure crystal calculations

import numpy as np

from modules.electronic_structure.structure import atoms


def structure_string():
    return string

def lattice_string():
    return string



header = '! Cubic Boron Nitride \n' + \
         '! https://journals.aps.org/prb/pdf/10.1103/PhysRevB.51.14705'

molecule = [atoms.Atom('N', [0,0,0]), atoms.Atom('B', [0.25, 0.25, 0.25])]

lattice_parameters =
monkhorst_pack = [4,4,4]