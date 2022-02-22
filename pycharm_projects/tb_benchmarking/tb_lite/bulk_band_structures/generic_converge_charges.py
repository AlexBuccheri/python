""" Module to converge an output value w.r.t. an input value, for an code

In principle, should work for TB lite, entos, exciting, and any other code
in which a Calculation sub-class is written for.
"""
import abc
from typing import Union
from pathlib import Path


# Valid directory types
path_type = Union[str, Path]


class Calculation(abc.ABC):
    """Abstract base class for a calculation.

    A calculation is expected to have:
        * A name
        * A (run) directory
        * One or more inputs
        * A parser for one or more outputs

    """
    def __init__(self, name: str,  directory: path_type):
        self.name = name
        self.directory = directory

    @abc.abstractmethod
    def set_input(self):
        """ Set one or more input files for calculation.
        """
        ...

    @abc.abstractmethod
    def write_input(self):
        """ Write one or more input files for calculation.
        """
        ...

    @abc.abstractmethod
    def parse_output(self):
        """ Parse one or more output files for calculation.
        """
        ...


class TbliteCalculation(Calculation):
    """ Subclass for a TBLite calculation.
    """
    def __init__(self, name: str, directory: path_type):
        super().__init__(name, directory)

    def set_input(self):
        """

        :return:
        """
        pass

    def write_input(self):
        """

        :return:
        """
        pass

    def parse_output(self):
        """

        :return:
        """
        pass
