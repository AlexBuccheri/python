"""
Given some converged ground state calculation set up GW calculations
with optimised bases, for a range of energy parameter cutoffs.
"""
from pathlib import Path
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

from gw_benchmark_inputs.set4_spd.A1_groundstate import converged_ground_state_input as A1_gs_input


def set_up_g0w0(root_path: str):

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

    # Increased rgkmax to 8. Converge the s-channel
    energy_cutoffs = {'zr': {0: [42,   80, 170, 250, 350],
                             1: [180, 180, 180, 180, 180],
                             2: [180, 180, 180, 180, 180],
                             3: [180, 180, 180, 180, 180]},

                      'o':  {0: [140, 140, 140, 140, 140],
                             1: [140, 140, 140, 140, 140],
                             2: [140, 140, 140, 140, 140]}
                      }

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


if __name__ == "__main__":
    info = "Converge the QP gap w.r.t. the s-channel of Zr, having taken the highest s core state into the valence"
    print(info)
    set_up_g0w0("/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set4_spd/zr_lmax3_o_lmax2_rgkmax8")
