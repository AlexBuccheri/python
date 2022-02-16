"""
Silicon GS settings
"""
from pymatgen.io.cif import CifParser


def cif_to_ase_atoms(file: str):
    """Extract CIF Data.

    :param file:
    :return:
    """
    structure = CifParser(file).get_structures()[0]
    print(structure)
    # structure.atomic_numbers
    # structure.lattice
    # structure.frac_coords

# def muffin_tin_radii():
#
cif_to_ase_atoms()