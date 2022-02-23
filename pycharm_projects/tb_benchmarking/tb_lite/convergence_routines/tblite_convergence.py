"""Example of all the boilerplate one needs to write to use the generic
convergence code.
"""
# Intrinsic packages
import os
from typing import Optional, Union, Tuple, List
import copy

# External packages
import ase.atoms

# Abstract classes and routine for convergence
from tb_lite.convergence_routines.generic_convergence import CalculationIO, path_type, \
    ConvergenceCriteria, converge
from tb_lite.src.runner import BinaryRunner, SubprocessRunResults

# TB Lite - specific
from tb_lite.src.parsers import parse_dftb_output
from tb_lite.src.dftb_input import DftbInput
from ase.io.dftb import write_dftb

# Temporary
from tb_lite.src.dftb_input import DftbInput, Hamiltonian
from tb_lite.crystal_references import cubic


class TbliteCalculation(CalculationIO):
    """ Subclass for a TBLite convergence calculation.
    """
    contains_geometry: bool

    def __init__(self,
                 name: str,
                 directory: path_type,
                 runner: BinaryRunner,
                 input: DftbInput,
                 atoms: Optional[ase.atoms.Atoms] = None):
        super().__init__(name, directory)
        self.runner = runner
        self.input = input
        self.atoms = atoms

        runner_methods = [method_name for method_name in dir(self.runner)
                          if callable(getattr(self.runner, method_name))]

        if "run" not in runner_methods:
            raise AttributeError('runner must contain a run method.')
        if not isinstance(input, DftbInput):
            raise ValueError('Require input to be of type DftbInput')
        if not isinstance(atoms, ase.atoms.Atoms):
            raise ValueError('Require atoms to be of type  ase.atoms.Atoms')

    def write_dftb_hsd(self) -> None:
        """ Write dftb_hsd.in

        TODO(Alex) Return if exception if file cannot be written
        """
        dftb_hsd: str = self.input.generate_dftb_hsd()
        # TODO(Alex) Use regex to check for geometry in input
        self.contains_geometry = False

        file_name = os.path.join(self.directory, "dftb_hsd.in")
        with open(file_name, "w") as fid:
            fid.write(dftb_hsd)

    def write_geometry(self) -> None:
        file_name = os.path.join(self.directory, "geometry.gen")
        write_dftb(file_name, self.atoms)

    def write_inputs(self):
        """ Write all input files required by DFTB+.

        Required file:
          * dftb_hsd.in

        Optional file:
          * geometry.gen (if not specified in input)
        """
        self.write_dftb_hsd()
        if not self.contains_geometry: self.write_geometry()

    def run(self) -> SubprocessRunResults:
        """ Wrapper for simple BinaryRunner.

        :return: Subprocess results or NotImplementedError.
        """
        return self.runner.run()

    def parse_output(self) -> Union[dict, FileNotFoundError]:
        """ Parse energies from DFTB+'s "detailed.out" in eV.

        Note, it does not matter when I parse as long as I am explicit
        w.r.t. what I compare in the convergence class.

        :return: Dictionary of final energies, or FileNotFoundError.
        """
        file_name = os.path.join(self.directory, "detailed.out")
        try:
            with open(file_name, "r") as fid:
                file_str = fid.read()
            return parse_dftb_output(file_str)
        except FileNotFoundError:
            return FileNotFoundError(f"DFTB+ output file not found: {file_name}")


class TBliteConvergence(ConvergenceCriteria):
    """Sub-class defining input parameter to vary for convergence,
    the output value to measure the change in, and the method
    to evaluate the convergence.
    """
    def __init__(self, input: list, criteria: dict):
        super().__init__(input, criteria)

    def evaluate(self, current, prior) -> Tuple[bool, bool]:
        delta = current['Total energy'] - prior['Total energy']
        converged = abs(delta) <= self.criteria['Total energy']
        early_exit = False
        return converged, early_exit


def set_kgrid_in_tblite(k_grid: Union[int, List[int]], calculation: TbliteCalculation):
    """ Given a Tblite calculation, with an input of type DftbInput, modify the k_grid attribute
    and write the new input to file.

    Easiest approach is to create a new calculation instance, update the k_grid, then
    write the input to file (overwriting the existing one). Alternatively, one could
    simply modify the existing input file.

    :param k_grid: k-grid sampling, either as 3 int, or a single int.
    :param calculation: TbliteCalculation
    """
    new_k_grid = [k_grid] * 3 if isinstance(k_grid, int) else k_grid
    odd_grid_points = [not bool(x%2) for x in new_k_grid]
    if any(odd_grid_points):
        raise ValueError(f'Cannot specify odd k-grids: {new_k_grid}. \n'
                         f'This would also require modifying the k-weights with the func `set_kgrid_in_tblite`')
    new_calculation = copy.deepcopy(calculation)
    new_calculation.input.hamiltonian.k_grid = new_k_grid
    new_calculation.write_dftb_hsd()


# Example of How to Perform a Convergence Calculation
#
# Some Issues/Things to Note:
#  directory in runner and calculation need to be consistent
#  Units are defined by the parser

# directory = "ADD ME"
# runner = BinaryRunner(binary='dftb+', run_cmd=['./'], omp_num_threads=1, directory=directory, time_out=600)
# input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[4,4,4]))
# atoms = cubic.silicon()
#
# calculation = TbliteCalculation("TB Lite SCC Convergence", directory, runner, input, atoms)
# k_grids = [[4, 4, 4], [6, 6, 6], [8, 8, 8]]
# convergence = TBliteConvergence(k_grids, {'Total energy': 1.e-5})
# converge(calculation, convergence, set_kgrid_in_tblite)
