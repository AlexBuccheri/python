"""
 Module containing cubic crystal dictionaries with the signature:
  key : str
    crystal name
  value : str
    file path to cif

  boron nitride is tabulated.
"""
import numpy as np
from ase.atoms import Atoms

from tb_lite.src.utils import FileUrl


# Hexagonal crystals by bravais lattice
# Space groups: 168 - 194
# All space groups in this system begin with P

hexagonal_cifs = {'zinc_oxide': FileUrl('data/bulk_crystals/cifs/hexagonal/ZnO-WZ/ZnO_mp-2133_primitive.cif', "https://materialsproject.org/materials/mp-2133/"),
                  'molybdenum_disulfide': FileUrl('data/bulk_crystals/cifs/hexagonal/MoS2/MoS2_mp-2815_primitive.cif', 'https://materialsproject.org/materials/mp-2815/'),
                  'boron_nitride': {'wurzite': FileUrl('data/bulk_crystals/cifs/hexagonal/BN-WZ/BN_mp-2653_primitive.cif', "https://materialsproject.org/materials/mp-2653/")},
                  'unbuckled_graphite': FileUrl('data/bulk_crystals/cifs/hexagonal/Unbluckled_graphite/A_hP4_194_bc.cif', 'UNKNOWN'),
                  'tungsten_disulfide': FileUrl('data/bulk_crystals/cifs/hexagonal/WS2/WS2_mp-224_primitive.cif', 'https://materialsproject.org/materials/mp-224/'),
                  'cadmium_selenide': FileUrl('data/bulk_crystals/cifs/hexagonal/CdSe-hexagonal/CdSe_mp-1070_primitive.cif', 'https://materialsproject.org/materials/mp-1070/')
                  }


# Crystals without cif files
# Avoid additional tabulating and use cif_parser_wrapper to generate
# these dictionaries.

# def boron_nitride():
#
#     """
#     # TODO(Alex) Refactor to ASE
#     hexagonal boron nitride.
#
#     Notes
#       Space group: P63/mmc [194]
#       Primitive lattice vectors and atomic basis
#       Indirect and gap: 4.482 eV
#       https://materialsproject.org/materials/mp-984/
#     """
#
#     positions = [[1/3, 2/3, 1/4],
#                  [2/3, 1/4, 3/4],
#                  [1/3, 2/3, 3/4],
#                  [2/3, 1/3, 1/4]]
#     species = ['B', 'B', 'N', 'N']
#     bravais = 'hexagonal'
#     space_group = 194
#     lattice_parameters = {'a': Set(2.51242804, 'angstrom'), 'c': Set(7.70726501, 'angstrom')}
#     data = {'fractional':positions,
#             'species':species,
#             'lattice_parameters':lattice_parameters,
#             'space_group': ('', space_group),
#             'bravais': bravais,
#             'n_atoms': 4}
#     return data


