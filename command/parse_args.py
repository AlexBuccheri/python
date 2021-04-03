import subprocess
import os
import warnings
from pathlib import Path


def run_executable(args, input_file) -> dict:
    """
    Run an executable. Assumes CMake build system

    See here for a more extensive way to manage the env:
    https://stackoverflow.com/questions/2059482/python-temporarily-modify-the-current-processs-environment

    :param args
    :param input: input file name for executable
    :param build_type_string: build type str, as defined in CMake with ${CMAKE_BUILD_TYPE}
    :param executable_name: Name of executable, as defined in CMake using
           RUNTIME_OUTPUT_NAME property of set_target_properties(...)
    :return: Dictionary containing stdout, stderr and returncode (error if not 0).
             The standard output and errors are lists, where each element is a string
             (Conversion from a byte object occurs)
    """
    try:
        Path(args.exe).resolve(strict=True)
    except FileNotFoundError as error:
        raise error

    os.environ["OMP_NUM_THREADS"] = str(args.omp_num_threads)
    run_command = ['./' + args.exe, input_file]

    if 'mpi' or 'hybrid' in args.build_type:
        run_command = ['mpirun', '-np', str(args.np)] + run_command

    process = subprocess.run(run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if process.returncode != 0:
        warnings.warn('process returned to stderr: ' + " ".join(run_command))

    return {'stdout': process.stdout.decode("utf-8").split('\n'),
            'stderr': process.stderr.decode("utf-8").split('\n'),
            'returncode': process.returncode}
