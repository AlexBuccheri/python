"""
Tests for compare_result_with_reference

Note, tests must begin with 'test' to be run

From root, run  `pytest test_dict_tools.py`
"""

import pytest
import numpy as np

from dict_tools import compare_result_with_reference, get_hashable_entries


def test_results_approx_equal_ref():
    """
    Given some simple, nested dictionaries, assert the values in them are consistent to specified tolerances
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

@pytest.mark.xfail()
def test_with_failing_lattice_constant():
    """
    Given some simple, nested dictionaries, assert the values in them are consistent to specified tolerances
    lattice constant expected to fail
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


def test_keys_not_overwritten():
    """
    If putting nested dictionary elements into a dictionary with no nesting,
    one will overwrite any keys that appear more than once.

    Avoid this by putting nested dictionary elements into a list of form [{key1:value1}, {key2:value2}, {key1:value3}]
    """
    parsed_result = {'scf1': {'total_energy': -22.23878943},
                     'scf2': {'total_energy': -12.23878945}
                     }

    result = [{key: value} for key, value in get_hashable_entries(parsed_result)]

    assert result == [{'total_energy': -22.23878943}, {'total_energy': -12.23878945}]


def test_equivalent_keys_at_different_nested_levels():


    parsed_reference = {'setup': {'n_species': 1, 'max_scf': 3},
                        'lattice_constants':np.array([1.2, 1.3, 1.4]),
                        'scf1':{'total_energy': -12.23878943},
                        'scf2': {'total_energy': -12.23878943}
                        }

    noise = 5.e-9

    parsed_result = {'setup': {'n_species': 1, 'max_scf': 3},
                        'lattice_constants':np.array([1.2, 1.3, 1.4]),
                        'scf1':{'total_energy': -12.23878943},
                        'scf2': {'total_energy': -12.23878945}
                        }

    tolerance = {'n_species': 0,                # Must always be constant
                 'max_scf':   1,                # Can vary by ±1
                 'lattice_constants': 1.e-8,    # Numerical tolerance
                 'total_energy': 1.e-10          # Diff of 2.e-8 in total_energy of scf2
                 }

    assert noise < tolerance['lattice_constants'], 'noise must be less than the tolerance of the lattice constants'

    compare_result_with_reference(parsed_result, parsed_reference, tolerance, debug_mode=True)
