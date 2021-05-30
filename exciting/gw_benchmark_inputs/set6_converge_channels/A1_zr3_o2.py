"""
Given some converged ground state calculation set up GW calculations
with optimised bases, for a range of energy parameter cutoffs.
"""
import shutil
from collections import OrderedDict
from distutils.dir_util import copy_tree

from parse.lorecommendations_parser import parse_lorecommendations
from parse.parse_linengy import parse_lo_linear_energies
from parse.parse_basis_xml import parse_basis_as_string
from parse.set_gw_input import GWInput
from process.optimised_basis import DefaultLOs
from job_schedulers import slurm

from gw_benchmark_inputs.input_utils import write_file, write_input_file_with_gw_settings, write_optimised_lo_basis, \
    restructure_energy_cutoffs

from gw_benchmark_inputs.set6_converge_channels.A1_groundstate import converged_ground_state_input as A1_gs_input


def set_up_g0w0(root_path: str, energy_cutoffs: dict):

    # Material
    species = ['zr', 'o']
    l_max = {'zr': 3, 'o': 2}

    gw_root = write_input_file_with_gw_settings(root_path,
                                                A1_gs_input,
                                                GWInput(taskname="g0w0", nempty=2000, ngridq=[2, 2, 2], skipgnd=False,
                                                        n_omega=32)
                                                )

    # Default basis settings
    default_linear_energies = parse_lo_linear_energies(root_path + "/groundstate")
    default_los = {'zr': DefaultLOs(default_linear_energies['zr'], energy_tol=1.5),
                   'o': DefaultLOs(default_linear_energies['o'], energy_tol=1.5)}

    # Default basis strings with .format tags
    default_basis_string = {'zr': parse_basis_as_string(root_path + "/groundstate/Zr.xml"),
                            'o': parse_basis_as_string(root_path + "/groundstate/O.xml")}

    # LO recommendation energies
    lorecommendations = parse_lorecommendations(root_path + '/lorecommendations.dat', species)

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

    species_basis_string = "".join(s.capitalize() + str(l_max[s]) + '_' for s in species)

    for ie, energy_cutoff in enumerate(restructure_energy_cutoffs(len(energy_cutoffs['zr'][0]), energy_cutoffs)):

        # Copy ground state directory to GW directory
        # Use an index not max energy, as the max energy does not change in 3/4 runs
        job_dir = gw_root + '/max_energy_i' + str(ie)
        print('Creating directory, with input.xml, run.sh and optimised basis:', job_dir)
        copy_tree(root_path + '/groundstate', job_dir)

        # Copy input.xml with GW settings
        shutil.copy(gw_root + "/input.xml", job_dir + "/input.xml")

        # New Slurm script
        slurm_directives['job-name'] = "gw_A1_lmax_" + species_basis_string + str(ie) + 'loEcutoff'
        write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))

        # Write optimised basis
        write_optimised_lo_basis('zr', l_max['zr'], energy_cutoff['zr'], lorecommendations['zr'],
                                 default_basis_string['zr'], default_los['zr'], job_dir)
        write_optimised_lo_basis('o', l_max['o'], energy_cutoff['o'], lorecommendations['o'],
                                 default_basis_string['o'], default_los['o'], job_dir)

    return


def converge_each_zr_channel(root: str):
    """
    Wrapper for the setting up a G0W0 run, converging the number of LOs for a given channel of Zr

    :param str root: root job directory
    """

    # TODO. These numbers need verifying from prior calcs
    s_converged = 200
    p_converged = 200
    d_converged = 200
    f_converged = 200

    base_energy_cutoffs = {'zr': {0: [s_converged] * 5,
                                  1: [p_converged] * 5,
                                  2: [d_converged] * 5,
                                  3: [f_converged] * 5},

                           'o':  {0: [140] * 5,
                                  1: [140] * 5,
                                  2: [140] * 5}
                           }

    # Add these, of the form
    # Ideal spd = [10, 10, 10]
    # Run s= [4, 6, 8, 10, 11], p = 10, d = 10
    # Run s= 10,  p = [4, 6, 8, 10, 11],  d = 10
    # Run s= 10,  p = 10,  d = [4, 6, 8, 10, 11]
    s_range = []
    p_range = []
    d_range = []
    f_range = []

    for l, lo_range in enumerate([s_range, p_range, d_range, f_range]):
        energy_cutoffs = base_energy_cutoffs
        energy_cutoffs['zr'][l] = lo_range
        directory = root + "/" + str(l) + "channel"
        set_up_g0w0(directory, energy_cutoffs)

    return


# TODO Once the above is done, implement/run these:
#
# With the optimal permutation, add one channel to Zr l=4, with 10 LOs
# With the optimal permutation, add one channel to O l=3, with 10 LOs
# - Confirm it has not effect
#
# Run it for q= 4,4,4 and 6,6,6 => Set up on Hawk
# Run a range of unoccupied states for q = [2,2,2]

if __name__ == "__main__":
    print("Converge the QP gap w.r.t. each LO channel of Zr")
    converge_each_zr_channel("/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8")
