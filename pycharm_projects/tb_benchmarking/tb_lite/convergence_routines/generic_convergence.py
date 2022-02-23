""" Module providing abstract classes and function to converge a calculation
w.r.t. some input value.
"""
import abc
from collections.abc import Iterable
from typing import Union, List, Callable, Optional, Tuple
from pathlib import Path

from tb_lite.src.runner import SubprocessRunResults


# Valid directory types
path_type = Union[str, Path]


class CalculationIO(abc.ABC):
    """Abstract base class for a calculation that is performed
    by writing input file/s and parsing the result from a file.

    An IO calculation is expected to have:
        * A name,
        * A working (run) directory,
        * A method to write all input files required to run the calculation,
        * A method to run the calculation,
        * A parser for the outputs of interest.
    """
    def __init__(self, name: str, directory: path_type):
        self.name = name
        self.directory = directory
        if not Path.is_dir(directory):
            raise NotADirectoryError(f'Not a directory: {directory}')

    @abc.abstractmethod
    def write_inputs(self) -> None:
        """ Write all input files required for calculation.
        """
        ...

    @abc.abstractmethod
    def run(self) -> SubprocessRunResults:
        """ Run the calculation.
        :return Subprocess result instance.
        """

    @abc.abstractmethod
    def parse_output(self, *args) -> Union[dict, FileNotFoundError]:
        """ Parse one or more output files for calculation.
        :return Dictionary of results.
        """
        ...


class ConvergenceCriteria(abc.ABC):
    """Abstract base class for performing a set of convergence calculations.

    Attributes correspond to input value to vary, and target value to check convergence against.
    Method should supply a convergence criterion or criteria w.r.t. the target value/s.
    """
    def __init__(self, input, criteria: dict):
        """ Initialise an instance of Convergence.

        :param input: A range of input values. Can be in any format, as long it's iterable.
        :param criteria: Dictionary of convergence criteria. {key:value} = {target: criterion}
        """
        self.input = input
        self.criteria = criteria
        if not isinstance(input, Iterable):
            raise ValueError('input must be iterable.')
        if len(input) <= 1:
            raise ValueError('input must have a length > 1')

    @staticmethod
    def check_target(func: Callable):
        """ Provide argument checking.

        :param func: evaluate method.
        :return: Modified evaluate method.
        """
        def func_with_target_check(self, current: dict, prior: dict):

            # Only expect for a failed run
            if isinstance(current, SubprocessRunResults):
                converged, early_exit = current.success, True
                return converged, early_exit

            set_current = set(current)

            # TODO(Alex) Should look at propagating the ValueErrors out
            if set_current != set(prior):
                return ValueError(f'Keys of current and prior results are inconsistent:'
                                  f'{set_current} != {set(prior)}')

            if set_current != set(self.criteria):
                raise ValueError(f'Keys of current result inconsistent with keys of convergence criteria:'
                                 f'{set(current)} != {set(self.criteria)}')

            return func(self, current, prior)

        return func_with_target_check

    @abc.abstractmethod
    @check_target
    def evaluate(self, current: dict, prior: dict) -> Tuple[bool, bool]:
        """ Evaluate a convergence criterion for each target.

        TODO(Alex) See if I can a) double decorate and b) apply a decorator to an abstract method
        Alternative would be to implement the start of the method, then inherit it in sub-classes.

        :param current: Dictionary containing current result/s
        :param prior: Dictionary containing prior result/s
        :return Tuple of bools indicating (converged, early_exit).
        """
        ...


def convergence_step(value,
                     calculation: CalculationIO,
                     set_value_in_input: Callable) -> Union[dict, SubprocessRunResults]:
    """ Perform a single calculation as part of a convergence test.

    :param value: Input value be varied to achieve convergence.
    :param calculation: calculation instance, with methods defined by CalculationIO
    :param set_value_in_input: Function defining how to change the input value.
    This could be achieved by copying calculation, changing value, then write to file OR
    it could be modifying an existing file.
    :return: Dictionary containing the output value used to evaluate when convergence w.r.t. value.
    Note, this can contain other data. The ConvergenceCriteria class determines how this is evaluated.
    If the calculation run fails, a SubprocessRunResults is returned instead.
    """
    set_value_in_input(value, calculation)
    subprocess_result = calculation.run()
    if subprocess_result.success:
        return calculation.parse_output()
    return subprocess_result


def converge(calculation: CalculationIO,
             convergence: ConvergenceCriteria,
             set_value_in_input: Callable[[any, CalculationIO], Union[any, None]]) -> List[tuple]:
    """ Converge a calculation output with respect to an input parameter.

    calculation defines the methods to:
        Write all inputs required,
        Run the calculation,
        Parse the result/s required to measure convergence.

    convergence defines:
        * The input parameter to vary (and its range of values),
        * The target output/s to check,
        * The criterion/criteria with which to evaluate convergence.

    The function `set_value_in_input`:
        Function defining how to change the input value for per convergence calculation.
        This is not a method of Calculation because it does not manage the state of
        a calculation (indeed, it could modify an input file already written), nor
        is it a method of Convergence.

    Note, this is entirely analogous to writing templated code in C++

    :param calculation: Calculation instance.
    :param convergence: Convergence parameters and criteria.
    :param set_value_in_input: Function that sets a value in the input
    defined by the calculation.

    :return: List of results. Each result is a tuple(input value, output)
    """
    # Write all required files to run the calculation
    calculation.write_inputs()

    # Initialise results by running the first value
    first_value = convergence.input[0]
    result = convergence_step(first_value, calculation, set_value_in_input)
    results = [(first_value, result, False, False)]

    if isinstance(result, SubprocessRunResults):
        return results

    for value in convergence.input[1:]:
        result = convergence_step(value, calculation, set_value_in_input)
        converged, early_exit = convergence.evaluate(result, results[-1][1])
        results.append((value, result, converged, early_exit))
        if converged or early_exit:
            return results

    return results
