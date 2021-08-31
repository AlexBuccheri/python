"""
Set 7 converging slowly, and can't easily visualise the numbers
Run systematic change in basis w.r.t.
 a) l_max channels
 b) N LOs per channel

l_max values: (4, 3) (5, 4) (6, 5) (7, 6)

Same as A1_zr_o.py in set 8 but
 a) Using high LO cut-offs i.e > 150 Ha. To check convergence
 b) Setup is for HAWK, using PBS Pro, and limited to 24 hours run-time

"""

import shutil
from collections import OrderedDict
from distutils.dir_util import copy_tree
import os

from job_schedulers.pbs_pro import set_pbs_pro_directives, set_pbs_pro
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
    Given an l_max pair, create G0W0 inputs for a specifid range of LO cut-offs per channel,
    as defined in set8/basis.py

    :param str root_path: Top level path to calculations
    :param List[str] species: List of species
    :param dict l_max: l_max associated with each species
    :return:
    """

    # Notes:
    # Ignore queue name and have HAWK assign the correct one.
    # omplace caused issues on test queue, so ignore it for now.

    omp = 16
    directives = set_pbs_pro_directives(time=[24,00,0],
                                        nodes=1,
                                        mpi_ranks_per_node=8,
                                        omp_threads_per_process=omp,
                                        cores_per_node=128,
                                        node_type='rome',
                                        job_name='GW_gs')

    env_vars = OrderedDict([('EXE', '/zhome/academic/HLRS/pri/ipralbuc/exciting-oxygen_release/bin/exciting_mpismp'),
                            ('OUT', 'terminal.out')
                            ])
    module_envs = ['intel/19.1.0', 'mkl/19.1.0', 'impi/19.1.0']
    # mpi_options = ['omplace -nt ' + str(omp)]
    mpi_options = []

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

    for ie in range(4, 6):
        energy_cutoff = energy_cutoffs[ie]

        # Copy ground state directory to GW directory
        job_dir = gw_root + '/max_energy_i' + str(ie)
        print('Creating directory, with input.xml, run.sh and optimised basis:', job_dir)
        print(root_path)
        copy_tree(root_path + '/groundstate', job_dir)

        # Copy input.xml with GW settings
        shutil.copy(gw_root + "/input.xml", job_dir + "/input.xml")

        # New Slurm script
        directives['N'] = "gw_A1_lmax_" + species_basis_string + str(ie)
        write_file(job_dir + '/run.sh', set_pbs_pro(directives, env_vars, module_envs, mpi_options))

        # Write optimised basis
        write_optimised_lo_bases(species, l_max, energy_cutoff, lorecommendations,
                                 default_basis_string, default_los, job_dir)


def generate_g0w0_inputs(root_path: str):

    species = ['zr', 'o']
    l_max_pairs = [{'zr': 4, 'o': 3}, {'zr': 5, 'o': 4}, {'zr': 6, 'o': 5}, {'zr': 7, 'o': 6}]

    for l_max in l_max_pairs:
        l_max_path = 'zr_lmax' + str(l_max['zr']) + '_o_lmax' + str(l_max['o']) + '_rgkmax8'
        full_path = os.path.join(root_path, l_max_path)
        input_for_lmax_pair(full_path, species, l_max)


generate_g0w0_inputs("/lustre/hpe/ws10/ws10.3/ws/ipralbuc-gw_scratch/gw_calculations/set8_part2/")
