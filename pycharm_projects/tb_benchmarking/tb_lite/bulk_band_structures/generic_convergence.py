""" Module to converge an output value w.r.t. an input value, for an code

In principle, should work for TB lite, entos, exciting, and any other code
in which a Calculation sub-class is written for.
"""
import abc
import os.path
from typing import Union, List, Callable, Optional
from pathlib import Path

import ase.atoms

from tb_lite.src.runner import BinaryRunner, SubprocessRunResults
from tb_lite.src.parsers import parse_dftb_output
from tb_lite.src.dftb_input import DftbInput
from ase.io.dftb import write_dftb


# Valid directory types
path_type = Union[str, Path]


class CalculationIO(abc.ABC):
    """Abstract base class for a calculation that is performed
    by writing input file/s and parsing the result from a file.

    An IO calculation is expected to have:
        * A name
        * A (run) directory
        * One or more input files
        * A method to run the calculation
        * A parser for the outputs of interest

    Every explicit method should return a result or an exception.
    The type of result is left to the implementation of the sub-class.
    """
    def __init__(self, name: str,  directory: path_type):
        self.name = name
        # TODO(Alex) Check the directory exists
        self.directory = directory

    @abc.abstractmethod
    def write_input(self):
        """ Write one or more input files for calculation.
        """
        ...

    @abc.abstractmethod
    def run(self) -> SubprocessRunResults:
        """ Run the calculation.
        """

    @abc.abstractmethod
    def parse_output(self, *args):
        """ Parse one or more output files for calculation.
        """
        ...


class TbliteCalculation(CalculationIO):
    """ Subclass for a TBLite calculation.

    TODO(Alex) This should be a bit more specialised for convergence calculations
    and provide the relevant methods/API.

    Requirements of arguments:
     runner is an instance of BinaryRunner. Although in principle, it could be
     anything with a run method, that returns a SubprocessRunResults.
    """
    def __init__(self,
                 name: str,
                 directory: path_type,
                 runner: BinaryRunner,
                 input: DftbInput,
                 atoms: Optional[ase.atoms.Atoms]=None):
        super().__init__(name, directory)
        self.runner = runner
        if not isinstance(runner, BinaryRunner):
            raise ValueError('runner must be an instance of BinaryRunner, such that it'
                             'has the expected attributes and run method.')
        self.input = input
        self.atoms = atoms

    def write_input(self):
        """ Write DFTB+'s input files:
         * dftb_hsd.in
         * geometry.gen (if not specified in input)
        """
        dftb_hsd: str = self.input.generate_dftb_hsd()
        file_name = os.path.join(self.directory, "dftb_hsd.in")
        with open(file_name, "w") as fid:
            fid.write(dftb_hsd)

        # TODO(Alex) Add check if geometry tag in dftb_hsd via regex. If missing, then call"
        file_name = os.path.join(self.directory, "geometry.gen")
        write_dftb(file_name, self.atoms)

        # TODO(Alex) Return an error if one is not able to write either file.

    def run(self) -> Union[SubprocessRunResults, NotImplementedError]:
        """ Wrapper for simple BinaryRunner.

        :return: Subprocess results or NotImplementedError.
        """
        try:
            run_results = self.runner.run()
            return run_results
        except NotImplementedError:
            methods = [method_name for method_name in dir(self.runner)
                       if callable(getattr(self.runner, method_name))]
            return NotImplementedError(f'runner does not have a .run() method: {methods}')

    def parse_output(self) -> Union[dict, FileNotFoundError]:
        """ Parse energies from DFTB+'s "detailed.out" in eV.

        :return: Dictionary of final energies, or FileNotFoundError.
        """
        file_name = os.path.join(self.directory, "detailed.out")
        try:
            with open(file_name, "r") as fid:
                file_str = fid.read()
            return parse_dftb_output(file_str)
        except FileNotFoundError:
            return FileNotFoundError(f"DFTB+ output file not found: {file_name}")


# TODO Alex
# Still need:
# set_value_in_input
# convergence_criterion
# Write a convergence class - Only
# First call to calculation writes all files
# Convergence wrapper really modifies one input parameter, rewrites that files and parses the target output.


def set_kgrid_in_tblite(k_grid: Union[int, List[int]], calculation: TbliteCalculation):
    """ Given a Tblite calculation, with an input of type DftbInput
    modify the k_grid attribute.

    NOTE: Could be a lambda.

    :param kgrid: k-grid sampling, either as 3 int, or a single int.
    :param calculation:
    :return: Mutate calculation.
    """
    _k_grid = [k_grid] * 3 if isinstance(k_grid, int) else k_grid
    calculation.input.hamiltonian.k_grid = _k_grid
    return calculation


class InputParameter:
    def __init__(self, name: any, values: List[any]):
        self.name = name
        if not isinstance(values, list):
            self.values = [values]
        else:
            self.values = values


def converge(calculation: CalculationIO,
             input_parameter: InputParameter,
             output_parameter: any,
             set_value_in_input: Callable[[any, CalculationIO], CalculationIO]) -> List[tuple]:
    """ Given a list of input_parameters converge an output_parameter.

    :param calculation: Calculation instance.
    :param set_value_in_input: Function that sets a value in the input
    defined by the calculation.

    :return: List of results. Each result is a tuple(input value, output)
    """
    results = []
    for value in input_parameter.values:
        set_value_in_input(value, calculation)
        calculation.write_input()
        subprocess_result = calculation.run()
        if subprocess_result.success:
            result = calculation.parse_output()
            output = result.get(output_parameter)
            results.append((value, output))
        else:
            results.append((value, subprocess_result))

    return results
