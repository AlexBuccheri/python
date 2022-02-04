"""
Hard-coded silicon inputs:
    Converged charges
    Band path
    - Check the format of band.out
https://dftbplus-recipes.readthedocs.io/en/latest/basics/bandstruct.html
"""
import os.path
import subprocess


import seekpath
import ase

from tb_lite.src.runner import SubprocessRunResults


def k_path_info_with_seekpath():
    """
    https://seekpath.readthedocs.io/en/latest/maindoc.html#how-to-use
    Returns
     {'GAMMA': [0.0, 0.0, 0.0], 'X': [0.5, 0.0, 0.5], 'L': [0.5, 0.5, 0.5], 'W': [0.5, 0.25, 0.75],
     'W_2': [0.75, 0.25, 0.5], 'K': [0.375, 0.375, 0.75], 'U': [0.625, 0.25, 0.625]}
     [('GAMMA', 'X'), ('X', 'U'), ('K', 'GAMMA'), ('GAMMA', 'L'), ('L', 'W'), ('W', 'X')]

     # Path from materials project: https://materialsproject.org/materials/mp-149/
    :return:
    """
    cell = [[0., 2.7145, 2.7145], [2.7145, 0., 2.7145], [2.7145, 2.7145, 0.]]
    positions = [[0.00, 0.00, 0.00], [0.25, 0.25, 0.25]]
    unique_indices = [14, 14]
    structure = (cell, positions, unique_indices)

    path_info = seekpath.get_path(structure, True, 'hpkot')
    points = path_info['point_coords']
    paths = path_info['path']


def k_path_with_ase():
    """
    https://wiki.fysik.dtu.dk/ase/gettingstarted/tut04_bulk/bulk.html
    :return:
    """
    cell = ase.atoms.Cell([[0., 2.7145, 2.7145], [2.7145, 0., 2.7145], [2.7145, 2.7145, 0.]])
    path = cell.bandpath()
    special_points = path.special_points
    k_path = path.path
    print(special_points)
    print(k_path)
    print(path.kpts)


k_path_with_ase()


si_input1 = """
Geometry = GenFormat {
    2  F
 Si
    1 1    0.0000000000        0.0000000000        0.0000000000
    2 1    0.2500000000        0.2500000000        0.2500000000
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.00000000000000    2.71450000000000    2.71450000000000
    2.71450000000000    0.00000000000000    2.71450000000000
    2.71450000000000    2.71450000000000    0.00000000000000
}

Hamiltonian = xTB {
  Method = "GFN1-xTB"
  SCC = Yes
  SccTolerance = 1e-06
  KPointsAndWeights = SupercellFolding {
    8 0 0
    0 8 0
    0 0 8
    0.5 0.5 0.5
  }
}

ParserOptions {
  ParserVersion = 7
}

"""


# a) Choose a k-path in fractional coordinates - try a small path to confirm
# b) See if this works with the xTB Hamiltonian

si_input2 = """

Geometry = GenFormat {
    2  F
 Si
    1 1    0.0000000000        0.0000000000        0.0000000000
    2 1    0.2500000000        0.2500000000        0.2500000000
    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
    0.00000000000000    2.71450000000000    2.71450000000000
    2.71450000000000    0.00000000000000    2.71450000000000
    2.71450000000000    2.71450000000000    0.00000000000000
}

Hamiltonian = xTB {
  Method = "GFN1-xTB"
  SCC = Yes
  ReadInitialCharges = Yes
  MaxSCCIterations = 1
  KPointsAndWeights = Klines {
     1   0.375  0.375  0.75   # K
    20   0.0    0.0    0.0    # Gamma
    20   0.5    0.5    0.5    # L
  }
}

ParserOptions {
  ParserVersion = 7
}

"""

# TODOs
# Need to determine how to automatically scale the band paths
# - Can write path.kpts from ASE
# - Doesn't look like DFTB+ can just take k-points (entos can though?)
# - I can map the path from ASE to DFTB+, then to be consistent with entos, parse the kpoint path from
# - the stdout and use that in the entos calcs.
# Automate input set-up
# Automate running each job (as two separate jobs)
# - copy to my machine, as I need dp_bands
# Create plotting, with labelling





