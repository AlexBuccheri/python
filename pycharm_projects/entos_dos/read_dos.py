import json
import numpy as np

# Keys associated with DOS output
_keys = ['energy_grid', 'fermi_level', 'total_density_of_states']


def get_dos(file_name: str, named_result=None):
    """
    Get the density of states results from a Qcore JSON output

    If named_result is not provided, it is assumed
    that the named_result for the dos output contains
    the substring 'dos'

    :param file_name : File name string
    :param named_result : Optional, named result string
    :return : dictionary for Qcore DOS output
    """

    with open(file_name) as fid:
        data = json.load(fid)
    fid.close()

    if named_result is None:
        keys = [key for key in data.keys()]
        named_result = [key for key in keys if 'dos' in key][0]

    return data[named_result]


def resolve_energy_grid(dos: dict) -> np.ndarray:
    """
    Resolve energy grid for the density of states

    In Qcore, the grid is constructed with arma::linspace, which
    Generates a vector with N elements; the values of the elements are
    linearly spaced from start to (and including) end

    :param dos: Qcore JSON output for DOS results
    :return: Linear energy grid
    """
    print(dos['energy_grid'].keys())
    grid = np.linspace(dos['energy_grid']['min'],
                       dos['energy_grid']['max'],
                       dos['energy_grid']['n_points'],
                       endpoint=True)
    return grid
