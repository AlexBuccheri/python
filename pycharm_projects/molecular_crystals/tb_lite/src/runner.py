"""
Binary runner and results class
"""
from typing import List, Optional
from pathlib import Path
import os
import subprocess


class SubprocessRunResults:
    """
    Results returned from subprocess.run()
    """

    def __init__(self, stdout, stderr, return_code: int):
        self.stdout = stdout
        self.stderr = stderr
        self.return_code = return_code
        self.success = return_code == 0


class BinaryRunner:
    """
    Compose a run command, and run a binary
    """

    def __init__(self,
                 binary: str,
                 run_cmd: List[str],
                 omp_num_threads: int,
                 time_out: int,
                 directory='./',
                 args=None
                 ) -> None:
        """
        :param str binary: Binary name prepended by full path
        :param List[str] run_cmd: Run commands sequentially as a list. For example:
          * For serial: ['./']
          * For MPI:   ['mpirun', '-np', '2']
        :param int omp_num_threads: Number of OMP threads
        :param int time_out: Number of seconds before a job is defined to have timed out
        :param List[str] args: Optional binary arguments
        """
        if args is None:
            args = ['']
        self.binary = binary
        self.directory = directory
        self.run_cmd = run_cmd
        self.omp_num_threads = omp_num_threads
        self.time_out = time_out
        self.args = args

        if not os.path.isfile(self.binary):
            raise FileNotFoundError("Binary does not exist: " + self.binary)

        if not Path(directory).is_dir():
            raise OSError("Run directory does not exist: " + directory)

        try:
            i = run_cmd.index('np')
            mpi_processes = eval(run_cmd[i + 1])
            assert type(mpi_processes) == int, "Number of MPI processes should be an int"
            assert mpi_processes > 0, "Number of MPI processes must be > 0"
        except ValueError:
            # .index will return ValueError if 'np' not found (serial and omp calculations)
            pass

        # TODO(Alex) if directory of type py._path.local.LocalPath, does it need explicitly converting
        # before passing to subprocess?

        assert omp_num_threads > 0, "Number of OMP threads must be > 0"

        assert time_out > 0, "time_out must be a positive integer"

    def _compose_execution_list(self) -> list:
        """
        Generate a complete list of strings to pass to subprocess.run, to execute the calculation.

        For example, given:
          ['mpirun', '-np, '2'] + ['binary.exe'] + ['>', 'std.out']

        return ['mpirun', '-np, '2', 'binary.exe', '>', 'std.out']
        """
        return self.run_cmd + [self.binary] + self.args

    def run(self, directory: Optional[str] = None, execution_list: Optional[list] = None) -> SubprocessRunResults:
        """
        Run a binary.

        :param str directory: Optional Directory in which to run the execute command.
        :param Optional[list] execution_list: Optional List of arguments required by subprocess.run. Defaults to None.
        """

        if directory is not None:
            directory = self.directory

        if not Path(directory).is_dir():
            raise OSError("Run directory does not exist: " + directory)

        if execution_list is None:
            execution_list = self._compose_execution_list()

        my_env = {**os.environ, "OMP_NUM_THREADS": str(self.omp_num_threads)}
        result = subprocess.run(execution_list,
                                env=my_env,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                timeout=self.time_out,
                                cwd=directory)

        return SubprocessRunResults(result.stdout, result.stderr, result.returncode)
