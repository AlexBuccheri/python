"""
 Module containing cubic crystal dictionaries with the signature:

 Cubic crystals by bravais lattice
 Space groups: 195-230.

"""
import numpy as np
from ase.atoms import Atoms

from tb_lite.src.utils import FileUrl

# Any space group beginning with I
bcc_cifs = {'potassium': FileUrl("bulk_cifs/cubic/BCC/potassium/K_mp-58_primitive.cif",
                                 "https://materialsproject.org/materials/mp-58/"),
            'sio2': FileUrl("bulk_cifs/cubic/BCC/sio2/SiO2_mp-1188220_primitive.cif",
                            "https://materialsproject.org/materials/mp-1188220/"),
            }

# Any space group beginning with F
fcc_cifs = {
    'copper': FileUrl("bulk_cifs/cubic/FCC/Cu/Cu_mp-30_primitive.cif",
                      "https://materialsproject.org/materials/mp-30/"),
    'pbs': FileUrl("bulk_cifs/cubic/FCC/PbS/PbS_mp-21276_primitive.cif",
                   "https://materialsproject.org/materials/mp-21276/"),
    'palladium': FileUrl("bulk_cifs/cubic/FCC/Pa/Pa_mp-10740_primitive.cif",
                         "https://materialsproject.org/materials/mp-10740/"),
    'boron_nitride': FileUrl("bulk_cifs/cubic/FCC/BN-cubic/BN_mp-1639_primitive.cif",
                             "https://materialsproject.org/materials/mp-1639/"),
    'sodium_chloride': FileUrl("bulk_cifs/cubic/FCC/NaCl/NaCl_mp-22862_primitive.cif",
                               "https://materialsproject.org/materials/mp-22862/"),
    'magnesium_oxide': FileUrl("bulk_cifs/cubic/FCC/MgO/MgO_mp-1265_primitive.cif",
                               "https://materialsproject.org/materials/mp-1265/"),
    'lithium_hydride': FileUrl("bulk_cifs/cubic/FCC/LiH/LiH_mp-23703_primitive.cif",
                               "https://materialsproject.org/materials/mp-23703/"),
    'silicon': FileUrl("bulk_cifs/cubic/FCC/Si/Si_mp-149_primitive.cif",
                       "https://materialsproject.org/materials/mp-149/")
}

# FCC but the lattice vectors are cubic
conventional_fcc_cifs = {
    'copper': FileUrl("bulk_cifs/cubic/FCC/Cu/Cu_mp-30_conventional_standard.cif",
                      "https://materialsproject.org/materials/mp-30/"),
    'sodium_chloride': FileUrl("bulk_cifs/cubic/FCC/NaCl/NaCl_mp-22862_conventional_standard.cif",
                               "https://materialsproject.org/materials/mp-22862/"),
    'magnesium_oxide': FileUrl("bulk_cifs/cubic/FCC/MgO/MgO_mp-1265_conventional_standard.cif",
                               "https://materialsproject.org/materials/mp-1265/"),
    'silicon': FileUrl("bulk_cifs/cubic/FCC/Si/Si_mp-149_conventional_standard.cif",
                       "https://materialsproject.org/materials/mp-149/")
}


# -----------------------------------------------------------------------
# Crystals without cif files
# Avoid additional tabulating and use cif_parser_wrapper to generate
# these dictionaries.
# -----------------------------------------------------------------------


def silicon() -> Atoms:
    """ Silicon Crystal

    Notes
      Ref: https://doi.org/10.1103/PhysRevB.24.6121
      Experimental value in Table I
      bravais = 'fcc'
      space_group = 227
    """
    cell = 5.429 * np.array([[0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]])

    atoms = Atoms(symbols=['Si', 'Si'],
                  scaled_positions=[[0.00, 0.00, 0.00], [0.25, 0.25, 0.25]],
                  cell=cell,
                  pbc=True)

    return atoms


def diamond() -> Atoms:
    """ Diamond Crystal

    Notes
      Ref: https://doi.org/10.1103/PhysRevB.24.6121
      Experimental value in Table I
      bravais = 'fcc'
      space_group = 227
    """
    cell = 3.567 * np.array([[0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]])

    atoms = Atoms(symbols=['C', 'C'],
                  scaled_positions=[[0.00, 0.00, 0.00], [0.25, 0.25, 0.25]],
                  cell=cell,
                  pbc=True)

    return atoms


def germanium() -> Atoms:
    """Germanium Crystal

    Notes
      Ref: https://doi.org/10.1103/PhysRevB.24.6121
      Experimental value in Table I
      bravais = 'fcc'
      space_group = 227
    """
    cell = 5.652 * np.array([[0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]])

    atoms = Atoms(symbols=['Ge', 'Ge'],
                  scaled_positions=[[0.00, 0.00, 0.00], [0.25, 0.25, 0.25]],
                  cell=cell,
                  pbc=True)

    return atoms
