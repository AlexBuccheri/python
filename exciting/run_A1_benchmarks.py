"""
Given some converged ground state calculation set up GW calculations
with optimised bases, for a range of energy parameter cutoffs.
"""
from pathlib import Path
import shutil

from parse.lorecommendations_parser import parse_lorecommendations
from parse.parse_linengy import parse_lo_linear_energies
from parse.parse_basis_xml import parse_basis_as_string
from parse.set_gw_input import GWInput, set_gw_input_string
from process.optimised_basis import DefaultLOs, filter_lo_functions, generate_optimised_basis_string
from benchmark_inputs.A1_groundstate import converged_ground_state_input


def write_file(file_name, string):
    fid = open(file_name, "w")
    fid.write(string)
    fid.close()
    return


# TODO Write some wrappers. Finish documentation
def set_up_g0w0(root_path:str):

    species = ['zr', 'o']
    l_max = {'zr': 4, 'o': 3}

    # GW root
    gw_root = root_path + "/zr_lmax" + str(l_max['zr']) + "_o_lmax" + str(l_max['o']) + "_rgkmax7"
    Path(gw_root).mkdir(parents=True, exist_ok=True)

    # exciting input file
    gw_input = GWInput(taskname="g0w0", nempty=600, ngridq=[2, 2, 2], skipgnd=False)
    gw_input_string = set_gw_input_string(converged_ground_state_input,  gw_input)
    write_file(gw_root + "/input.xml", gw_input_string)

    # Default basis settings
    default_linear_energies = parse_lo_linear_energies(root_path + "/groundstate")
    default_los = {}
    default_los['zr'] = DefaultLOs(default_linear_energies['zr'], energy_tol=0.1)
    default_los['o'] = DefaultLOs(default_linear_energies['o'],  energy_tol=0.1)

    # Default basis strings with tags
    default_basis_string = {}
    default_basis_string['zr'] = parse_basis_as_string(root_path + "/groundstate/Zr.xml")
    default_basis_string['o'] = parse_basis_as_string(root_path + "/groundstate/O.xml")

    # LO recommendation energies
    lorecommendations = parse_lorecommendations(root_path + '/lorecommendations.dat', species)

    for energy_cutoff in [20]:

        job_dir = gw_root + '/max_energy_' + str(energy_cutoff)
        # TODO Rather than make, this needs to copy the ground state
        Path(job_dir).mkdir(parents=True, exist_ok=True)
        shutil.copy(gw_root + "/input.xml", job_dir + "/input.xml")
        # TODO Write run script

        # Same for both species and each l-channel
        optimised_lo_cutoff_zr = [energy_cutoff] * (l_max['zr'] + 1)
        optimised_lo_cutoff_o = [energy_cutoff] * (l_max['o'] + 1)

        optimised_lo_zr = filter_lo_functions(lorecommendations['zr'], default_los['zr'], optimised_lo_cutoff_zr)
        optimised_lo_o = filter_lo_functions(lorecommendations['o'], default_los['o'], optimised_lo_cutoff_o)

        basis_string_zr = generate_optimised_basis_string(default_basis_string['zr'], optimised_lo_zr, max_matching_order=1)
        basis_string_o = generate_optimised_basis_string(default_basis_string['o'], optimised_lo_o, max_matching_order=1)

        write_file(job_dir + '/Zr.xml', basis_string_zr)
        write_file(job_dir + '/O.xml', basis_string_o)

    return

set_up_g0w0("/Users/alexanderbuccheri/gw_calc")