"""
Module for converging the SCC of Bulk Systems of Interest
"""

# Define inputs for each system
# Maybe just do this manually

from tb_lite.src.dftb_input import DftbInput, Hamiltonian


def get_material(material_name: str) -> DftbInput:
    """ TB lite settings, with converged inputs



    :param material:
    :return:
    """
    materials = {}

    materials['silicon'] = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6,
                                                             k_grid=[8, 8, 8]))

    try:
        material = materials[material_name]
    except KeyError:
        raise KeyError(f'Material {material_name} does not have converged settings for band structure defined')

    return material
