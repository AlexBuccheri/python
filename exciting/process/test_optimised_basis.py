import numpy as np

from optimised_basis import filter_default_function, filter_default_functions, filter_high_energy_functions


def test_filter_low_node_function_l0():
    """
    l = 0
    Expect (-956.66885744, -379.21022537, -14.29585798  -1.39085409) to be semi-core
    and valence states (as last lo in yhr default basis has energy -1.39). The first energy
    parameter in the optimised_basis_energies should therefore be 3.41674973
    """

    # energy parameter recommendations for lo's in l=0 channel
    lorecommendations_l0 = np.array([-956.66885744, -379.21022537, -14.29585798,  -1.39085409,
                                         3.41674973,   11.85763359,  23.55534883,  38.22251395,
                                        55.71839349,   75.96314341,  98.90623003, 124.50690254,
                                       152.73632137,  183.57226937, 216.99619737, 252.99363124,
                                        291.5523508,  332.66195298, 376.31355739, 422.49933633,
                                       471.21242405])

    # Representative of what one would extract from LINENGY.OUT, for a given atom
    linear_energies = {0: -1.390000000,
                       1: -0.510000000,
                       2:  0.330000000,
                       3:  1.000000000,
                       4:  1.000000000}

    # Only meaningful if linear_energy does not match any recommendations
    arbitrary_negative_number = -18945
    max_nodes_default_basis = arbitrary_negative_number

    l_value = 0
    tol = 0.1


    optimised_basis = filter_default_function(lorecommendations_l0,
                                             linear_energies[l_value],
                                             max_nodes_default_basis,
                                             l_value,
                                             energy_tolerance=tol)




    # optimised_basis = filter_default_function(lorecommendations_l0,
    #                                          linear_energies[l_value],
    #                                          max_nodes_default_basis,
    #                                          l_value,
    #                                          energy_tolerance=tol)
    #
    # assert optimised_basis.l_value == l_value, "l_value consistent"
    #
    # # Although does this correspond to a valence or conduction state? Not very meaningful
    # assert optimised_basis.first_n_nodes == 4, "first radial function for optimised basis has 4 nodes"
    #
    # np.testing.assert_approx_equal(optimised_basis.energies[0], 3.41674973,
    #                                err_msg="first energy parameter in optimised basis should come "
    #                                        "immediately after -1.39085409 in the reference default basis")
    #
    # n = optimised_basis.first_n_nodes
    # assert lorecommendations_l0[n] == optimised_basis.energies[0], \
    #     "Energy at index defined by n (nodes) in lorecommendations should correspond to the first " \
    #     "energy recommendation in optimised_basis"

    return


def test_filter_high_energy_functions_l0():

    lorecommendations_l0 = np.array([-956.66885744, -379.21022537, -14.29585798,  -1.39085409,
                                         3.41674973,   11.85763359,  23.55534883,  38.22251395,
                                        55.71839349,   75.96314341,  98.90623003, 124.50690254,
                                       152.73632137,  183.57226937, 216.99619737, 252.99363124,
                                        291.5523508,  332.66195298, 376.31355739, 422.49933633,
                                       471.21242405])

    lorecommendations = [lorecommendations_l0]

    l_value = 0
    energy_cutoff = 20.
    optimised_recommendations_index = filter_high_energy_functions(lorecommendations, energy_cutoff)

    assert optimised_recommendations_index[l_value] == 5, \
        "Index of highest-retained energy that is smaller than energy_cutoff"

    i_max = optimised_recommendations_index[l_value] + 1
    optimised_recommendations = lorecommendations_l0[0:i_max]

    np.testing.assert_array_equal(optimised_recommendations,
                                  np.array([-956.66885744, -379.21022537, -14.29585798,
                                            -1.39085409,      3.41674973,  11.85763359]),
                                  err_msg="expect all energies >= energy_cutoff to be filtered")

    return


def test_filter_low_node_function_l2():
    """
    l = 2
    Expect (-5.85379645) to be a valence state
    As no energies in lorecommendations_l2 match the energy parameter 0.33
    of the lo function in the default basis, the lowest recommendation should be 0.90588827
    """

    lorecommendations_l2 = np.array([-5.85379645, 0.90588827, 5.63201191, 13.67427176,
                                     24.54381621, 38.11976977, 54.32810605, 73.13682452,
                                     94.52501993, 118.47340661, 144.97035705, 174.00450227,
                                     205.56671156, 239.64966278, 276.24676259, 315.35282809,
                                     356.96353638, 401.07535931, 447.68537818, 496.79100982,
                                     548.38988507])

    # Representative of what one would extract from LINENGY.OUT, for a given atom
    linear_energies = {0: -1.390000000,
                       1: -0.510000000,
                       2: 0.330000000,
                       3: 1.000000000,
                       4: 1.000000000}

    l_value = 2
    tol = 0.1

    # Only meaningful if linear_energy does not match any recommendations
    max_nodes_default_basis = 0
    assert max_nodes_default_basis == len([linear_energies[l_value]]) - 1, \
        "Max nodes associated with a set of los == number of energy parameters in a given l-channel - 1"

    matched_energy_parameter = np.amin(np.abs(lorecommendations_l2 - linear_energies[l_value])) <= tol
    assert not matched_energy_parameter, "linear energy from LINENGY (0.330) should not be in lorecommendations_l2"

    optimised_basis = filter_default_function(lorecommendations_l2,
                                             linear_energies[l_value],
                                             max_nodes_default_basis,
                                             l_value,
                                             energy_tolerance=tol)

    assert optimised_basis.l_value == l_value, "l_value consistent"

    assert optimised_basis.first_n_nodes == max_nodes_default_basis + 1, \
        "first radial function for optimised basis has 1 node more than that of the max in the default basis"

    np.testing.assert_approx_equal(optimised_basis.energies[0], 0.90588827,
                                   err_msg="first energy parameter in optimised basis")

    n = optimised_basis.first_n_nodes
    assert lorecommendations_l2[n] == optimised_basis.energies[0], \
        "Energy at index defined by n (nodes) in lorecommendations should correspond to the first " \
        "energy recommendation in optimised_basis_energies"

    return


def test_filter_low_node_function_l4():
    """
    l = 4
    All energies exceed 1 Ha (the typical maximum value for an lo in the default basis)
    so expect all states to be conduction states.
    As no energies in lorecommendations_l4 match the energy parameter 1.00 (see linear_energies[4])
    of the lo function in the default basis the lowest recommendation should be determined by
    the function with the largest number of nodes in the default basis. Then go one higher with the recommendation
    This is somewhat meaningless to me
    Maybe 4.77449195 should indeed actually be 1.000, and Andris returns something erroneously in lorecommendations?
    """

    lorecommendations_l4 = np.array([4.77449195, 11.17705627, 20.04560097, 31.43639778,
                                     45.38373733, 61.8787738, 80.91492554, 102.49376441,
                                     126.60673075, 153.24771388, 182.40983259, 214.08650332,
                                     248.27299967, 284.9649611, 324.15906139, 365.85254928,
                                     410.04302461, 456.72834158, 505.90639589, 557.57511102,
                                     611.73245043])

    # Representative of what one would extract from LINENGY.OUT, for a given atom
    #TODO refactor this format (see basis.py)
    linear_energies = {0: -1.390000000,
                       1: -0.510000000,
                       2: 0.330000000,
                       3: 1.000000000,
                       4: 1.000000000}

    l_value = 4
    tol = 0.1

    # Only one energy parameter for the l=4 channel => only one function
    # function number of nodes start counting at 0 =>
    #TODO refactor getting this (see basis.py)
    max_nodes_default_basis = 0

    optimised_basis = filter_default_function(lorecommendations_l4,
                                             linear_energies[l_value],
                                             max_nodes_default_basis,
                                             l_value,
                                             energy_tolerance=tol)

    assert optimised_basis.l_value == l_value, "l_value consistent"
    n = optimised_basis.first_n_nodes

    assert optimised_basis.first_n_nodes == 1, \
        "only one energy parameter defining this l-channel of the default basis, hence " \
        "one basis function with 0 nodes => first recommendation should have 1 node"

    assert n == max_nodes_default_basis + 1, \
        "first radial function for optimised basis has 1 node more than that of the max in the default basis"

    np.testing.assert_approx_equal(optimised_basis.energies[0], 11.17705627,
                                   err_msg="first energy parameter in optimised basis")

    assert lorecommendations_l4[n] == optimised_basis.energies[0], \
        "Energy at index defined by n (nodes) in lorecommendations should correspond to the first " \
        "energy recommendation in optimised_basis_energies"

    return


# TODO Handle deliberately failing test
# def test_filter_low_node_functions():
#     """
#     Different API that accepts multiple l-channels
#     :return:
#     """
#
#     # Could be any numbers
#     lorecommendations_l4 = np.array([4.77449195, 11.17705627, 20.04560097, 31.43639778,
#                                      45.38373733, 61.8787738, 80.91492554, 102.49376441,
#                                      126.60673075, 153.24771388, 182.40983259, 214.08650332,
#                                      248.27299967, 284.9649611, 324.15906139, 365.85254928,
#                                      410.04302461, 456.72834158, 505.90639589, 557.57511102,
#                                      611.73245043])
#
#     linear_energies = {0: -1.390000000,
#                        1: -0.510000000,
#                        2: 0.330000000,
#                        3: 1.000000000,
#                        4: 1.000000000}
#
#     junk_default_lo_nodes = {0:  -1,
#                              1:  -1,
#                              2:  -1,
#                              3:  -1}
#
#     l_value = 5
#     tol = 0.1
#
#     optimised_basis = filter_default_functions(lorecommendations_l4, linear_energies, junk_default_lo_nodes)
#
#     return
