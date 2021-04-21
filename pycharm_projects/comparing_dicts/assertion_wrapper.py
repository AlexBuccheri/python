import numpy as np
import collections


default_tolerance = {int: 0, float: 1.e-8}


def assert_equal(data, reference, atol=None, message=None):
    """
    Assertion wrapper
    """
    assert type(data) == type(reference), "Reference and output data must be of same type"

    if atol is None:
        atol = default_tolerance[type(reference)]

    if message is None:
        message = ''

    if type(reference) in [int, float, np.ndarray]:
        np.testing.assert_allclose(data, reference, atol = atol, err_msg = message)

    elif type(reference) == list:
        # Check reference only contains hashable/immutable data types:
        # int, float, decimal, complex, bool, string, tuple, range, frozenset, bytes
        hashable_mask = [isinstance(entry, collections.Hashable) for entry in reference]
        assert all(hashable_mask) == True, 'ref_value list contains mutable data type'
        np.testing.assert_allclose(np.array(data), np.array(reference), atol = atol, err_msg = message)

    else:
        assert False, "Reference data type not handled by assertions"
