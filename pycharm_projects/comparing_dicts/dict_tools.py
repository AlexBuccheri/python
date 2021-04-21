from collections import abc
import numpy as np
import math

from assertion_wrapper import assert_equal


def get_hashable_entries(nested:dict):
    """
    Create a dictionary of all hashable values, removing nesting.

    Note, this will fail if a cross edge or a back edge is included in the dict
    (i.e. a self-reference of the form ...). However, for the use case of flattening parsed results,
    this should never be an issue.

    # TODO Question of a) how robust is this and b) how fast is it?
    https://stackoverflow.com/questions/10756427/loop-through-all-nested-dictionary-values

    :param nested:
    :return:
    """
    for key, value in nested.items():
        if isinstance(value, abc.Mapping):
            yield from get_hashable_entries(value)
        else:
            yield (key, value)


def check_key_consistency(a:dict, b:dict) -> bool:
    keys_a = [key for key in a.keys()]
    keys_b = [key for key in b.keys()]
    assert keys_a == keys_b, "Keys in reference and output data are not consistent "


def compare_result_with_reference(parsed_result:dict, parsed_reference:dict, tolerance:dict, debug_mode=False):
    """
    # TODO Used OrderedDicts
    # TODO Fill missing tolerances with defaults

    Provides some way of comparing nested data from two dictionaries.
    If an assertion fails, the code crashes.

    Not sure if there's a way of making this work for two differing nested structures.
    i.e. if we switch to NOMAD parser, existing references break, if we switch to structured output (removing
    the need for NOMAD parsers, because parsing becomes trivial), existing references also break.
    Implying we need to generate new reference data each time we change parser OR structure of output data

    :param dict result: results
    :param dict reference: references
    :param dict tolerance: tolerances for result to be considered equal to reference
    :param bool debug_mode: Perform debug mode checks
    """

    result = {key: value for key, value in get_hashable_entries(parsed_result)}
    reference = {key: value for key, value in get_hashable_entries(parsed_reference)}

    if debug_mode:
        check_key_consistency(reference, result)

    for key in result.keys():
        assert_equal(result[key], reference[key], tolerance[key], message=key)

