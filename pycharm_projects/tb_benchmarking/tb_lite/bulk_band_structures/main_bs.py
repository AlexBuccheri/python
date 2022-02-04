"""
Generate Bulk Band Structures using DFTB+ (where possible) and TB Lite (xTB1, xTB2)
for these materials:
* IV: Si, Ge, Diamond - GOT
* TODO  GET Wide-gap: ZnO, TiO2 (both phases), ZrO2, WO3
* TODO  GET BN (cubic and hexagonal)
* TODO  GET III-V: GaN, GaP, GaAs, InN, InP, InAs
* TODO  GET Narrow band-gap II-VI: PbS, PbSe, PbTe
* TODO  GET MoS2, WS2
"""
from tb_lite.src.dftb_input import DftbInput, Hamiltonian


def get_material(material_name: str) -> DftbInput:
    """ TB lite settings, with converged inputs



    :param material:
    :return:
    """
    materials = {}

    materials['silicon'] = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6,
                                                             k_points=[8, 8, 8]))

    try:
        material = materials[material_name]
    except KeyError:
        raise KeyError(f'Material {material_name} does not have converged settings for band structure defined')

    return material


# Define iterator for directories


def generate_input_for_converged_charges():
    return None


def generate_band_structure_input():
    # klines are given in fractional coordinates
    # Hamiltonian = DFTB {
    #   Scc = Yes
    #   ReadInitialCharges = Yes
    #   MaxSCCIterations = 1
    #
    #   # ...
    #
    #   KPointsAndWeights = Klines {
    #     1   0.5   0.5  -0.5    # Z
    #    20   0.0   0.0   0.0    # G
    #    45   0.0   0.0   0.5    # X
    #    10   0.25  0.25  0.25   # P
    #   }
    # }
    return None


def run_calculations():
    """
    Run calculation 1 to get converged charges
    Run calculation 2 to get the band structure
    :return:
    """
    return None


def parse_band_structure():
    """
    Parse:
        a) Direct gap
        b) CBm - VBM
        c) Band structure, and return in a plottable format
    :return:
    """
    return None
