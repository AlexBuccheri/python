"""
Given a converged basis from set9/A1_zr_o.py, systematically reduce the
basis function LOs per channel, such that a minimal basis for converged
calculations can be found.

Converged basis directory:
/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8/gw_q222_omeg32_nempty3000/max_energy_i1

Basis Settings
Zr MT radius = 2 and O MT radius = 1.6
rgkmax = 8
L_max=(6,5) LO cut-off = 100 Ha

Tolerate ~ 1 meV change in QP direct gap got a change of basis in any given
l-channel.
"""
import shutil
from collections import OrderedDict
from distutils.dir_util import copy_tree
import os
from typing import List
import copy

from job_schedulers import slurm
from parse.lorecommendations_parser import parse_lorecommendations
from parse.parse_linengy import parse_lo_linear_energies
from parse.parse_basis_xml import parse_basis_as_string
from parse.set_gw_input import GWInput
from process.optimised_basis import DefaultLOs

# I/O Utilities
from gw_benchmark_inputs.input_utils import write_input_file_in_root, write_optimised_lo_bases, \
    restructure_energy_cutoffs, write_file
# Settings
from gw_benchmark_inputs.set9.basis import converged_ground_state_input as A1_gs_input


# TODO Make a note of the directory structure


def gw_input(root_path: str,
             ground_state_dir: str,
             energy_cutoffs: List[int],
             species=['zr', 'o'],
             l_max={'zr': 6, 'o': 5}):
    """

    :param str root_path: Top level path to calculations
    :param ground_state_dir: Path to groundstate directory
    :param List[int] energy_cutoffs: LO energy cut-offs
    :param species:
    :param l_max:
    """

    # Slurm script settings
    env_vars = OrderedDict([('EXE', '/users/sol/abuccheri/exciting/bin/excitingmpismp'),
                            ('OUT', 'terminal.out')
                            ])
    module_envs = ['intel/2019']
    slurm_directives = slurm.set_slurm_directives(time=[0, 72, 0, 0],
                                                  partition='all',
                                                  exclusive=True,
                                                  nodes=4,
                                                  ntasks_per_node=2,
                                                  cpus_per_task=18,
                                                  hint='nomultithread')

    # Need some excessively large number for nempty => exciting takes upper bound
    write_input_file_in_root(root_path,
                             A1_gs_input,
                             GWInput(taskname="g0w0",
                                     nempty=3000,
                                     ngridq=[2, 2, 2],
                                     skipgnd=False,
                                     n_omega=32,
                                     freqmax=1.0)
                             )
    # Default basis settings
    default_linear_energies = parse_lo_linear_energies(ground_state_dir)
    default_los = {'zr': DefaultLOs(default_linear_energies['zr'], energy_tol=0.8),
                   'o': DefaultLOs(default_linear_energies['o'],  energy_tol=0.8)}

    # Default basis strings with .format tags
    default_basis_string = {'zr': parse_basis_as_string(os.path.join(ground_state_dir, "Zr.xml")),
                            'o': parse_basis_as_string(os.path.join(ground_state_dir, "O.xml"))}

    # LO energies
    lorecommendations = parse_lorecommendations(root_path + '/../../lorecommendations.dat', species)

    n_energies_per_channel = 3
    energy_cutoffs = restructure_energy_cutoffs(n_energies_per_channel, energy_cutoffs)

    species_basis_string = "_".join(s.capitalize() + str(l_max[s]) for s in species)

    for ie, energy_cutoff in enumerate(energy_cutoffs):
        # Copy ground state directory to GW directory
        job_dir = root_path + '/max_energy_' + str(ie)
        print('Creating directory, with input.xml, run.sh and optimised basis:', job_dir)
        copy_tree(ground_state_dir, job_dir)

        # Copy input.xml with GW settings
        shutil.copy(root_path + "/input.xml", job_dir + "/input.xml")

        # New Slurm script
        slurm_directives['job-name'] = "gw_A1_lmax_" + species_basis_string + str(ie)
        write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))

        # Write optimised basis
        write_optimised_lo_bases(species, l_max, energy_cutoff, lorecommendations,
                                 default_basis_string, default_los, job_dir)


def generate_g0w0_inputs(root_path: str, ground_state_path: str):
    """
    Generate an input for a given channel
    :param root_path:
    :return:
    """
    n_changes = 3

    original_cutoffs = {'zr': {0: [100] * n_changes,
                               1: [100] * n_changes,
                               2: [100] * n_changes,
                               3: [100] * n_changes,
                               4: [100] * n_changes,
                               5: [100] * n_changes,
                               6: [100] * n_changes},

                        'o': {0: [100] * n_changes,
                              1: [100] * n_changes,
                              2: [100] * n_changes,
                              3: [100] * n_changes,
                              4: [100] * n_changes,
                              5: [100] * n_changes}
                        }

    # TODO Elaborate on how this format means something different
    # Each cut-off has been chosen by inspection to
    # remove 1, 2 and 3 of the highest LOs
    channel_cutoffs = {'zr': {0: [90, 70, 50],  #
                              1: [90, 70, 50],
                              2: [90, 70, 50],
                              3: [80, 60, 45],
                              4: [75, 55, 40],
                              5: [85, 65, 45],
                              6: [75, 55, 40]},

                       'o': {0: [70, 45, 25],
                             1: [85, 60, 40],
                             2: [75, 50, 30],
                             3: [90, 60, 40],
                             4: [75, 50, 30],
                             5: [85, 60, 40]}
                       }

    # Need to remove this from the BASIS!!!
    #        <lo l="2">
    #         <wf matchingOrder="0" trialEnergy="5.64" searchE="false"/>
    #         <wf matchingOrder="1" trialEnergy="5.64" searchE="false"/>
    #        </lo>

    # TODO Note that slurm needs to be replaced for HAWK

    for l_channel, cutoffs_to_vary in channel_cutoffs['zr'].items():
        energy_cutoffs = copy.deepcopy(original_cutoffs)
        energy_cutoffs['zr'][l_channel] = cutoffs_to_vary
        print('Varying cutoff for channel ' + str(l_channel) + 'of Zr', energy_cutoffs['zr'])

        full_path = os.path.join(root_path, 'zr_channel' + str(l_channel))
        gw_input(full_path, ground_state_path, energy_cutoffs)
        energy_cutoffs.clear()

    for l_channel, cutoffs_to_vary in channel_cutoffs['o'].items():
        energy_cutoffs = copy.deepcopy(original_cutoffs)
        energy_cutoffs['o'][l_channel] = cutoffs_to_vary
        print('Varying cutoff for channel ' + str(l_channel) + 'of O', energy_cutoffs['zr'])

        full_path = os.path.join(root_path, 'o_channel' + str(l_channel))
        gw_input(full_path, ground_state_path, energy_cutoffs)
        energy_cutoffs.clear()


def main():
    root_path = "/users/sol/abuccheri/gw_benchmarks/A1_set9/"
    ground_state_path = os.path.join(root_path, 'zr_lmax6_o_lmax5_rgkmax8/groundstate')
    full_path = os.path.join(root_path, 'zr_lmax6_o_lmax5_rgkmax8/channel_refinement')
    generate_g0w0_inputs(full_path, ground_state_path)


if __name__ == "__main__":
    main()
