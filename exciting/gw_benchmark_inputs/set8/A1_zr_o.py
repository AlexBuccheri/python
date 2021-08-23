"""
Set 7 converging slowly, and can't easily visualise the numbers
Run systematic change in basis w.r.t.
 a) l_max channels
 b) N LOs per channel

l_max values: (4, 3) (5, 4) (6, 5) (7, 6)
Cut-off can probably go to ~ 150 Ha per channel

"""

import shutil
from collections import OrderedDict
from distutils.dir_util import copy_tree
from pathlib import Path
import os

from job_schedulers import slurm
from parse.lorecommendations_parser import parse_lorecommendations
from parse.parse_linengy import parse_lo_linear_energies
from parse.parse_basis_xml import parse_basis_as_string
from parse.set_gw_input import GWInput
from process.optimised_basis import DefaultLOs

# I/O Utilities
from gw_benchmark_inputs.input_utils import write_input_file, write_optimised_lo_bases, write_file
# Settings
from gw_benchmark_inputs.set8.basis import converged_ground_state_input as A1_gs_input, set_lo_channel_cutoffs


def input_for_lmax_pair(root_path: str, species: list, l_max: dict):
    """
    Given an l_max pair, create
    :param root_path:
    :param species:
    :param l_max:
    :return:
    """

    # Slurm script settings
    env_vars = OrderedDict([('EXE', '/users/sol/abuccheri/exciting/bin/excitingmpismp'),
                            ('OUT', 'terminal.out')
                            ])
    module_envs = ['intel/2019']
    slurm_directives = slurm.set_slurm_directives(time=[0, 24, 0, 0],
                                                  partition='all',
                                                  exclusive=True,
                                                  nodes=4,
                                                  ntasks_per_node=2,
                                                  cpus_per_task=18,
                                                  hint='nomultithread')

    # Need some excessively large number for nempty => exciting takes upper bound
    gw_root = write_input_file(root_path,
                               A1_gs_input,
                               GWInput(taskname="g0w0",
                                       nempty=2000,
                                       ngridq=[2, 2, 2],
                                       skipgnd=False,
                                       n_omega=32,
                                       freqmax=1.0)
                               )
    # Default basis settings
    default_linear_energies = parse_lo_linear_energies(root_path + "/groundstate")
    default_los = {'zr': DefaultLOs(default_linear_energies['zr'], energy_tol=0.8),
                   'o': DefaultLOs(default_linear_energies['o'],  energy_tol=0.8)}

    # Default basis strings with .format tags
    default_basis_string = {'zr': parse_basis_as_string(root_path + "/groundstate/Zr.xml"),
                            'o': parse_basis_as_string(root_path + "/groundstate/O.xml")}

    # LO energies
    lorecommendations = parse_lorecommendations(root_path + '/lorecommendations.dat', species)
    energy_cutoffs = set_lo_channel_cutoffs(l_max)

    species_basis_string = "_".join(s.capitalize() + str(l_max[s]) for s in species)

    for ie, energy_cutoff in enumerate(energy_cutoffs):
        # Copy ground state directory to GW directory
        job_dir = gw_root + '/max_energy_i' + str(ie)
        print('Creating directory, with input.xml, run.sh and optimised basis:', job_dir)
        print(root_path)
        copy_tree(root_path + '/groundstate', job_dir)

        # Copy input.xml with GW settings
        shutil.copy(gw_root + "/input.xml", job_dir + "/input.xml")

        # New Slurm script
        slurm_directives['job-name'] = "gw_A1_lmax_" + species_basis_string + str(ie)
        write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))

        # Write optimised basis
        write_optimised_lo_bases(species, l_max, energy_cutoff, lorecommendations,
                                 default_basis_string, default_los, job_dir)
        quit()


def generate_g0w0_inputs(root_path: str):

    species = ['zr', 'o']
    #l_max_pairs = [{'zr': 4, 'o': 3}, {'zr': 5, 'o': 4}, {'zr': 6, 'o': 5}, {'zr': 7, 'o': 6}]
    l_max_pairs = [{'zr': 4, 'o': 3}]

    for l_max in l_max_pairs:
        l_max_path = 'zr_lmax' + l_max['zr'] + '_o_lmax' + l_max['o'] + '_rgkmax8'
        full_path = os.path.join(root_path, l_max_path)
        Path(full_path).mkdir(parents=True, exist_ok=True)
        input_for_lmax_pair(full_path, species, l_max)


generate_g0w0_inputs("/users/sol/abuccheri/gw_benchmarks/A1_set8/")
