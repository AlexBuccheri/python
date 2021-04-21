"""
Tests for compare_result_with_reference

Note, tests must begin with 'test' to be run

From root, run  `pytest test_dict_tools.py`
"""

import numpy as np

from dict_tools import compare_result_with_reference


def test_results_approx_equal_ref():
    """
    Given some simple, nested dictionaries, assert the values in them are consistent.
    """

    parsed_reference = {'setup': {'n_species': 1, 'max_scf': 3},
                        'lattice_constants':np.array([1.2, 1.3, 1.4])}

    noise = 5.e-9

    parsed_result = {'setup': {'n_species': 1 , 'max_scf': 2},
                     'lattice_constants': np.array([1.2, 1.3, 1.4]) + noise}

    tolerance = {'n_species': 0,             # Must always be constant
                 'max_scf':   1,             # Can vary by ±1
                 'lattice_constants': 1.e-8  # Numerical tolerance
                 }

    assert noise < tolerance['lattice_constants'], 'noise must be less than the tolerance of the lattice constants'

    compare_result_with_reference(parsed_result, parsed_reference, tolerance, debug_mode=True)


def test_with_failing_lattice_constant():
    """
    Given some simple, nested dictionaries, assert the values in them are consistent.
    """

    parsed_reference = {'setup': {'n_species': 1, 'max_scf': 3},
                        'lattice_constants':np.array([1.2, 1.3, 1.4])}

    noise = 1.e-6

    parsed_result = {'setup': {'n_species': 1 , 'max_scf': 2},
                     'lattice_constants': np.array([1.2, 1.3, 1.4]) + noise}

    tolerance = {'n_species': 0,             # Must always be constant
                 'max_scf':   1,             # Can vary by ±1
                 'lattice_constants': 1.e-8  # Numerical tolerance
                 }

    assert noise > tolerance['lattice_constants'], \
        'noise > tolerance of the lattice constants, such that the comparison fails '

    compare_result_with_reference(parsed_result, parsed_reference, tolerance, debug_mode=True)

