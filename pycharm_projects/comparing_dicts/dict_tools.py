from collections import abc
import numpy as np
import math

from assertion_wrapper import assert_equal


def get_hashable_entries(nested:dict):
    """
    Create a dictionary of all hashable values, removing nesting.

    Note, this will fail if a cross edge or a back edge is included in the dict
    i.e. a self-reference of the form:

    A = {key1:value, key2:A}

    However, for the use case of flattening parsed results, this should never be an issue.

    # TODO Question of a) how robust is this and b) how fast is it?
    https://stackoverflow.com/questions/10756427/loop-through-all-nested-dictionary-values

    :param nested: Nested dictionary with a mix and hashable and non-hashable values.
    :return: Generator
    """
    for key, value in nested.items():
        if isinstance(value, abc.Mapping):
            yield from get_hashable_entries(value)
        else:
            yield (key, value)


# TODO This should be replaced by a check of key-consistency in the nested dicts
# Maybe easy to do with deepdiff
def check_key_consistency(a:list, b:list) -> bool:
    keys_a = [dictionary.keys() for dictionary in a]
    keys_b = [dictionary.keys() for dictionary in b]
    assert keys_a == keys_b, "Keys in reference and output data are not consistent "


def compare_result_with_reference(parsed_result:dict, parsed_reference:dict, tolerance:dict, debug_mode=False):
    """
    # TODO Fill missing tolerances with defaults

    Provides some way of comparing nested data from two dictionaries.

    It takes all hashable entries:

    dict = { key1:value1, key2:value2, key3:{key4:value4}}

    and puts them sequentially into a list, of the form:

        [{key1:value1}, {key2:value2}, {key4:value4}]

    Note, this is not serialising the data. In flattening the dictionary, we lose some key information
    (key3 in the above example). This is acceptable if the reference and result data have the same structure.
    As such, some pre-check that the keys of the nested dictionaries are consistent before comparing flattened data
    is required.

    TODO Add check that the keys of the nested dictionaries are consistent
    Only likely to occur if structure of output or parser changes, and the reference data has not bee updated.
    Can just be run in debug mode.

    Not sure if there's a way of making this work for two differing nested structures.
    i.e. if we switch to NOMAD parser, existing references break, if we switch to structured output (removing
    the need for NOMAD parsers, because parsing becomes trivial), existing references also break.
    Implying we need to generate new reference data each time we change parser OR structure of output data

    If an assertion fails, the script fails.

    :param dict result: results
    :param dict reference: references
    :param dict tolerance: tolerances for result to be considered equal to reference
    :param bool debug_mode: Perform debug mode checks
    """

    result = [{key: value} for key, value in get_hashable_entries(parsed_result)]
    reference = [{key: value} for key, value in get_hashable_entries(parsed_reference)]

    if debug_mode:
        check_key_consistency(reference, result)

    def key_from_single_entry(a:dict):
        return [x for x in a.keys()][0]

    def value_from_single_entry(a:dict):
        return [x for x in a.values()][0]

    for i in range(0, len(result)):

        result_key = key_from_single_entry(result[i])
        reference_key = key_from_single_entry(reference[i])
        assert result_key == reference_key

        result_value = value_from_single_entry(result[i])
        reference_value = value_from_single_entry(reference[i])

        assert_equal(result_value, reference_value, tolerance[reference_key], message=reference_key)

