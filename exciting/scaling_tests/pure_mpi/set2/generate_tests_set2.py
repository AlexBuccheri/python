"""
Generate pure MPI tests

Dune 3 has 16 nodes (of the same type) and 36 cores_per_node = 576 cores

```
"""
import numpy as np
import re
from collections import OrderedDict
from pathlib import Path

from job_schedulers import slurm
from scaling_tests.pure_mpi.set2 import inputs_set2


def write_file(file_name, string):
    fid = open(file_name, "w")
    fid.write(string)
    fid.close()
    return


def set_up_pure_mpi_scaling_tests(scaling_root:str):
    """
    Pure MPI scaling tests for GW, from 1 to 14 nodes, using all cores per node
    and no threading.

    GW is only MPI-parallelised over q-points => use:
    q = [8,8,8] and n_frequency points = 1, such that this can scale on up to 512 cores.

    """

    input_xml = inputs_set2.input_xml
    zr_basis_xml = inputs_set2.zr_basis_xml
    o_basis_xml = inputs_set2.o_basis_xml

    # Check GW input script settings
    match = re.search('nempty="(.+?)"', input_xml)
    n_empty = int(re.findall(r'\d+', match.group())[0])
    assert n_empty == 100, "n_empty != 100"

    match = re.search('ngridq="(.+?)"', input_xml)
    q_grid = [int(q) for q in re.findall(r'\d+', match.group())]
    assert q_grid == [8, 8, 8], "q_grid != [8, 8, 8]"

    # Slurm script settings
    env_vars = OrderedDict([('EXE', '/users/sol/abuccheri/exciting/bin/excitingmpismp'),
                            ('OUT', 'terminal.out'),
                            ('export MKL_NUM_THREADS', '1')
                            ])

    module_envs = ['intel/2019']

    # Cores per node
    ntasks_per_node = 36

    # OMP threads per MPI rank
    cpus_per_task = 1

    # Nodes to use in scaling tests
    # Dune 3 only appears to have 10 nodes available from nodes 181 - 196
    nodes = np.arange(1, 10+1)

    # These nodes differ in memory or processor to the rest of Dune 3
    # hence exclude 197 - 208
    exclude_nodes = ['node' + str(id) for id in range(197, 208 + 1)]

    # Timing in days, where key = node_count
    times = { 1: [4, 0, 0, 0],
              2: [4, 0, 0, 0],
              3: [4, 0, 0, 0],

              4: [2, 0, 0, 0],
              5: [2, 0, 0, 0],
              6: [2, 0, 0, 0],

              7: [1, 0, 0, 0],
              8: [1, 0, 0, 0],
              9: [1, 0, 0, 0],
             10: [1, 0, 0, 0],
             11: [1, 0, 0, 0],
             12: [1, 0, 0, 0],
             13: [1, 0, 0, 0],
             14: [1, 0, 0, 0]}

    for node_count in nodes:
        job_dir = scaling_root + '/n_nodes_' + str(node_count)
        print("Writing files to:", job_dir)

        Path(job_dir).mkdir(parents=True, exist_ok=True)

        write_file(job_dir + '/input.xml', input_xml)
        write_file(job_dir + '/Zr.xml', zr_basis_xml)
        write_file(job_dir + '/O.xml', o_basis_xml)

        slurm_directives = slurm.set_slurm_directives(job_name='scaling-pure-mpi-GW',
                                                      time=times[node_count],
                                                      partition='all',
                                                      exclusive=True,
                                                      nodes=node_count,
                                                      ntasks_per_node=ntasks_per_node,
                                                      cpus_per_task=cpus_per_task,
                                                      hint='nomultithread',
                                                      exclude=exclude_nodes)
        write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))

    return


set_up_pure_mpi_scaling_tests("/users/sol/abuccheri/gw_benchmarks/scaling/pure_mpi/set2")
