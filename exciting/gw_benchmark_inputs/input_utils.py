"""
Input utilities shared by all GW input generation scripts
"""
from pathlib import Path
from typing import List
import numpy as np

from parse.set_gw_input import GWInput, set_gw_input_string
from process.optimised_basis import DefaultLOs, filter_lo_functions, generate_optimised_basis_string


def write_file(file_name: str, string: str):
    """
    Basic wrapper around opening, writing and closing a file.

    :param str file_name: File name
    :param str string: String to write to file of file_name
    """
    fid = open(file_name, "w")
    fid.write(string)
    fid.close()


def write_input_file_with_gw_settings(root_path: str, gs_input: str, gw_input: GWInput) -> str:
    """
    Generate GW path, create a GW input file and write it

    :param str root_path: Path to the calculation directory
    :param str gs_input: ground state input, used as the base for the GW input
    :param GWInput gw_input: GW calculation input settings

    :return: str gw_root: Directory for GW calculations
    """
    q_str = ''.join(str(q) for q in gw_input.ngridq)
    gw_root = root_path + "/gw_q" + q_str + "_omeg" + str(gw_input.nomeg) + "_nempty"+str(gw_input.nempty)

    Path(gw_root).mkdir(parents=True, exist_ok=True)
    gw_input_string = set_gw_input_string(gs_input,  gw_input)
    write_file(gw_root + "/input.xml", gw_input_string)
    return gw_root


def write_optimised_lo_basis(species: str,
                             l_max: int,
                             energy_cutoff: dict,
                             lorecommendations: List[np.ndarray],
                             default_basis_string: str,
                             default_los: DefaultLOs,
                             job_dir: str):
    """
    Write the optimised LO basis for GW calculations, to file.

    The number of l-channels in the default basis determines the number of l-channels in the optimised basis.

    :param str species: lowercase species string.
    :param int l_max: l-value of the LO channel with max(l). Only passed to ensure the basis the user species is the
    one being used.
    :param dict energy_cutoff: energy cut-off associated with a given species, of the form {l: LO_cutoff energy}.
    :param List[array] lorecommendations: LO energy recommendations for all l-channels of a given species.
    :param str default_basis_string: Default basis strings with .format tags for where to insert optimised LOs,
    per l-channel.
    :param DefaultLOs default_los: Default LO basis details.
    :param str job_dir: Calculation directory.
    """
    assert len(lorecommendations) == 7, "Expect 7 l-channels for lorecommendations, " \
                                        "per species, as it is hard-coded in exciting"

    assert len(lorecommendations[0]) == 22, "Expect 22 entries per l-channel lorecommendations, " \
                                            "as it is hard-coded in exciting"

    assert l_max == len(default_los.linear_energies) - 1, \
        "l_max defined in 'set_up_g0w0' differs from the max l-channel present in the default basis"

    optimised_lo = filter_lo_functions(lorecommendations, default_los, energy_cutoff)
    basis_string = generate_optimised_basis_string(default_basis_string, optimised_lo, max_matching_order=1)

    basis_file = species.capitalize() + '.xml'
    write_file(job_dir + '/' + basis_file, basis_string)



def restructure_energy_cutoffs(n_energies:int, energy_cutoffs: dict) -> list:
    """
    Get in a more useful structure to iterate over

    TODO Automate to find n_energies

    """
    restructured_energies = []

    for inum in range(0, n_energies):
        data = {}
        for species, l_channels in  energy_cutoffs.items():
            data[species] = {l:energies[inum] for l, energies in l_channels.items()}
        restructured_energies.append(data)

    return restructured_energies



# # TODO Write this. Need to unpack in the same way I did in other instances.
# # Or just dump the dict in JSON
# def write_basis_settings(energy_cutoffs:dict):
#
#     species = [element for element in energy_cutoffs.keys()]
#     for element in species:
#     return

