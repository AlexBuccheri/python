"""
Converge "converged" optimised GW basis w.r.t. q-points

5s 12p 10d 8f
"""
from collections import OrderedDict
from distutils.dir_util import copy_tree

from gw_benchmark_inputs.input_utils import write_file, write_input_file_with_gw_settings
from parse.set_gw_input import GWInput, set_gw_input_string
from job_schedulers import slurm

from gw_benchmark_inputs.q_convergence.A1_groundstate import converged_ground_state_input as A1_gs_input
from gw_benchmark_inputs.q_convergence.optimised_basis import zr_basis, o_basis


def set_up_g0w0(root_path: str, q_grid: list):

    species = ['zr', 'o']
    l_max = {'zr': 3, 'o': 2}
    q_str = "".join(str(q) for q in q_grid)

    # GW input file
    gw_input_string = set_gw_input_string(A1_gs_input,
                                          GWInput(taskname="g0w0", nempty=2000, ngridq=q_grid, skipgnd=False, n_omega=32)
                                          )

    run_settings = {'222': {'nodes': 4,  'time': [0, 24, 0, 0]},
                    '444': {'nodes': 10, 'time': [0, 64, 0, 0]},
                    '666': {'nodes': 12, 'time': [0, 100, 0, 0]}
                    }

    # Slurm script settings
    env_vars = OrderedDict([('EXE', '/users/sol/abuccheri/exciting/bin/excitingmpismp'),
                            ('OUT', 'terminal.out')
                            ])
    module_envs = ['intel/2019']
    slurm_directives = slurm.set_slurm_directives(time=run_settings[q_str]['time'],
                                                  partition='all',
                                                  exclusive=True,
                                                  nodes=run_settings[q_str]['nodes'],
                                                  ntasks_per_node=2,
                                                  cpus_per_task=18,
                                                  hint='nomultithread')

    # Job directory
    job_dir = root + '/' + q_str

    # Write input, basis files and slurm file
    print('Creating directory, with input.xml, run.sh and optimised basis:', job_dir)
    copy_tree(root_path + '/groundstate', job_dir)

    # Copy input.xml with GW settings
    write_file(job_dir + "/input.xml", gw_input_string)

    # New Slurm script
    slurm_directives['job-name'] = "gw_A1_" + q_str
    write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))

    # Optimised bases
    write_file(job_dir + '/Zr.xml', zr_basis)
    write_file(job_dir + '/O.xml', o_basis)


if __name__ == "__main__":
    print("Converge the QP gap w.r.t. each LO channel of Zr")

    q_grids = [[2, 2, 2], [4, 4, 4], [6, 6, 6]]
    root = "/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/q_convergence"

    for q_grid in q_grids:
        q_str = "".join(str(q) for q in q_grid)
        set_up_g0w0(root + '', q_grid)

