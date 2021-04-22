import numpy as np
import collections


_default_tolerance = {int: 0, float: 1.e-8}


def assert_equal(data, reference, abs_tol=None, message=''):
    """
    Assertion wrapper

    TODO Test and extend for complex floats + arrays
    """
    assert type(data) == type(reference), "Reference and output data must be of same type"

    if abs_tol is None:
        abs_tol = _default_tolerance[type(reference)]

    if type(reference) in [int, float, np.ndarray]:
        np.testing.assert_allclose(data, reference, atol = abs_tol, err_msg = message)

    elif type(reference) == list:
        # Check reference only contains hashable/immutable data types:
        # int, float, decimal, complex, bool, string, tuple, range, frozenset, bytes
        hashable_mask = [isinstance(entry, collections.Hashable) for entry in reference]
        assert all(hashable_mask) == True, 'ref_value list contains mutable data type'
        np.testing.assert_allclose(np.array(data), np.array(reference), atol = abs_tol, err_msg = message)

    else:
        assert False, "Reference data type not handled by assertions"
