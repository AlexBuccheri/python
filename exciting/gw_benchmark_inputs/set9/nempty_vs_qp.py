"""
Generate inputs for a range of N empty states and extract the band gaps
"""
import shutil
from collections import OrderedDict
from typing import List
import os

from job_schedulers import slurm
from parse.set_gw_input import GWInput

# I/O Utilities
from gw_benchmark_inputs.input_utils import write_input_file_in_root, write_file
from gw_benchmark_inputs.set9.basis import converged_ground_state_input as A1_gs_input

# --------------------------
# Constants
# --------------------------
l_max = {'zr': 6, 'o': 5}

refined_basis_cutoffs = {'zr': {0: [70],
                                1: [50],
                                2: [50],
                                3: [100],
                                4: [40],
                                5: [65],
                                6: [100]},

                         'o': {0: [70],
                               1: [85],
                               2: [100],
                               3: [60],
                               4: [75],
                               5: [100]}
                         }

# Max n_empty taken from the refined basis calculation:
refined_basis_file_location = "/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8/refined"

n_empty_refined_basis = 1208


# --------------------------
# Routines
# --------------------------
def vary_nempty(root: str, basis_file_path: str, nempty_range: List[int]):
    """
    Generate inputs with varying numbers of empty states

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
                                                  ntasks_per_node=1,
                                                  cpus_per_task=36,
                                                  hint='nomultithread')

    for n_empty in nempty_range:

        # Job directory
        job_dir = os.path.join(root, str(n_empty))

        # GW INPUT
        gw_settings = GWInput(taskname="g0w0",
                              nempty=n_empty,
                              ngridq=[2, 2, 2],
                              skipgnd=False,
                              n_omega=32,
                              freqmax=1.0)
        write_input_file_in_root(job_dir, A1_gs_input, gw_settings)

        # Just copy the basis files and STATE.OUT
        for species in ['Zr', 'O']:
            species_file = species + ".xml"
            shutil.copy(basis_file_path + "/" + species_file, job_dir + "/" + species_file)
        shutil.copy(basis_file_path + "/STATE.OUT", job_dir + "/STATE.OUT")


        # Slurm script
        slurm_directives['job-name'] = "gw_nempty_" + str(n_empty)
        write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))


if __name__ == "__main__":
    root_path = "/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8/nempty_vs_qp"
    nempty_range = [200, 400, 600, 800, 1000, 1100]
    vary_nempty(root_path, refined_basis_file_location, nempty_range)
