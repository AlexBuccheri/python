"""
Generate SLURM script
"""
from typing import Optional
from collections import OrderedDict


def set_slurm_script(slurm_directives:OrderedDict, env_vars:OrderedDict, module_envs:Optional[list]=None) -> str:
    """
    Generate simple slurm submission script, suitable for hybrid MPI+OMP applications

    :param slurm_directives: Ordered dictionary of slurm directives
    :param env_vars: Ordered dictionary of environment variables to set
    :param module_envs: Optional list of modules to load
    :return: slurm script string
    """

    env_keys = [key for key in env_vars.keys()]
    assert 'EXE' in env_keys, "EXE must be specified in env_vars"

    script = "#!/bin/bash \n\n"
    slurm_prefix = "#SBATCH --"

    for directive, setting in slurm_directives.items():
        if setting == '':
            script += slurm_prefix + directive + '\n'
        else:
            script += slurm_prefix + directive + "=" + setting + '\n'

    script += '\n'

    if module_envs is not None:
        for module in module_envs:
            script +=  "module load " + module + '\n'
        script += '\n'

    for key, setting in env_vars.items():
        script += key + '=' + setting + '\n'

    script += """
cd ${SLURM_SUBMIT_DIR}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# $SLURM_NTASKS = n_nodes * ntasks-per-node == total number of MPI processes
mpirun -np $SLURM_NTASKS $EXE """

    if 'OUT' in env_keys:
        script += ' > $OUT'

    return script


def time_to_string(time:list):
    """
    Time list to string
    :param time: List of length 4
    :return: string with format 1-00:20:00
    """
    assert len(time) == 4, "expect [days, hours, mins, secs]"

    time_string = ''
    for i,t in enumerate(reversed(time[1:])):
        t_str = str(t)
        if len(t_str) == 1:
            t_str = '0' + t_str
        time_string = ':' + t_str + time_string

    time_string = str(time[0]) + '-' + time_string[1:]

    return time_string


def set_slurm_directives(job_name:Optional[str]='default_name',
                         time:Optional[list]=None,
                         partition:Optional[str]=None,
                         exclusive:Optional[bool]=True,
                         nodes:Optional[int]=1,
                         ntasks_per_node:Optional[int]=None,
                         cpus_per_task:Optional[int]=None,
                         hint:Optional[str]='nomultithread') -> OrderedDict:
    """
    Set ordered dictionary of slurm directives

    :param job_name:
    :param time:
    :param partition:
    :param exclusive:
    :param nodes:
    :param ntasks_per_node:
    :param cpus_per_task:
    :param hint:
    :return:
    """

    time_str = time_to_string(time)

    slurm_directives = OrderedDict([('job-name', job_name),
                                    ('time', time_str),
                                    ('partition', partition),
                                    ('nodes', str(nodes)),
                                    ('ntasks-per-node', str(ntasks_per_node)),
                                    ('cpus-per-task', str(cpus_per_task)),
                                    ('hint', hint)
                                    ])

    # Add directives that don't have a RHS
    if exclusive:
        slurm_directives['exclusive']  = ''

    return slurm_directives
