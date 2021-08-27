from collections import OrderedDict

from pbs_pro import set_pbs_pro_directives, set_pbs_pro


def test_set_pbs_pro_without_optionals():

    ref = """#!/bin/bash 

#PBS -N HAWK Rome Pure MPI
#PBS -q ADD ME
#PBS -l walltime= 24:00:00
#PBS -l select=2:ncpus=128:mpiprocs=128:ompthreads=1:
#PBS -M abuccheri@physik.hu-berlin.de
#PBS -j oe

module purge 
module load intel/2019

EXE=/users/sol/abuccheri/exciting/bin/excitingmpismp
OUT=terminal.out
cd $PBS_O_WORKDIR 
export OMP_NUM_THREADS=1

mpirun -np 256 $EXE > $OUT 2>&1
    """

    directives = set_pbs_pro_directives([24,0,0],
                                        'ADD ME',
                                         nodes=2,
                                         mpi_ranks_per_node=128,
                                         omp_threads_per_process=1,
                                         job_name='HAWK Rome Pure MPI')

    env_vars = OrderedDict([('EXE', '/users/sol/abuccheri/exciting/bin/excitingmpismp'),
                            ('OUT', 'terminal.out')
                            ])
    module_envs = ['intel/2019']

    # TODO Think about how to add this
    mpi_args = OrderedDict([('pinning', 'omplace')])

    script = set_pbs_pro(directives, env_vars, module_envs)
    print(script)
    # assert script == ref

# TODO Write me
def test_pbs_pro_with_all_args():
