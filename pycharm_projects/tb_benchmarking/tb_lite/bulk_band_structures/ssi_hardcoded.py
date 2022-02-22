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



# def generate_inputs():
#     """ Generate some xTB inputs
#     """
#     root = '/Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/check'
#     for material in ['silicon', 'germanium', 'diamond']:
#         directory = os.path.join(root, material)
#         Path(directory).mkdir(parents=True, exist_ok=True)
#         write_an_xtb1_input(directory, material)
#
#
# def write_an_xtb1_input(directory: str, material: str):
#     """Given a material, write xTB input files to a directory.
#     """
#     atoms, input = get_material_xtb1(material)
#     write_dftb(directory + "/geometry.gen", atoms)
#     with open(directory + "/dftb_in.hsd", "w") as fid:
#         fid.write(input.generate_dftb_hsd())
#
#
# def get_material_xtb1(material_name: str) -> Tuple[ase.atoms.Atoms, DftbInput]:
#     """ TB lite Inputs for bulk crystals of Interest
#
#     Manually-converged inputs:
#      * Si
#      * Ge
#      * Diamond
#     Inputs Requiring Convergence:
#      * Graphite
#      * Graphene
#      * ZrO2
#      * ZnO
#      * WS2
#      * GaAs
#      * InAs
#      * PbS
#      * BN -cubic, hexagonal and wurzite
#      * MoS2
#      * WS2
#      * TiO2 rutile
#      * TiO2 anatase
#
#     :param material_name: Material key
#     :return: Atoms object and DFTB+ Input object.
#     """
#     # Group IV elemental crystals
#     if material_name == 'silicon':
#         input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[8, 8, 8]))
#         return silicon(), input
#
#     elif material_name == 'germanium':
#         input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[12, 12, 12]))
#         return germanium(), input
#
#     elif material_name == 'diamond':
#         input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[8, 8, 8]))
#         return diamond(), input
#
#     elif material_name == 'zinc_oxide':
#         file_path = hexagonal.hexagonal_cifs.get(material_name).file
#         # Should read with pymatgen and convert to ASE atoms... Doesn't seem particularly robust
#         # atoms = read_cif(file_path, index=0, primitive_cell=True)
#         # print(list(atoms))
#         # input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[]))
#         return [], []
#
#     else:
#         print(f'material_name is not valid: {material_name}')
#
#     #elif material_name == '':
#         #
#     # Graphite
#     # Graphene
#     # ZrO2
#     #
#     # WS2
#     # GaAs
#     # InAs
#     # PbS
#     # BN -cubic, hexagonal and wurzite
#     # MoS2
#     # WS2
#     # TiO2 rutile
#     # TiO2 anatase
#     # Consider high throughput on the set form the Sotti paper ~ 400 crystals.
#
#
#
#
# file = hexagonal.hexagonal_cifs.get('zinc_oxide').file
# material, atoms = cif_to_ase_atoms(file)
#
# print(material)






