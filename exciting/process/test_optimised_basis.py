import numpy as np

from optimised_basis import filter_default_functions, filter_high_energy_functions



def test_filter_low_node_functions():
    # energy parameter recommendations for lo's. Must have all l-channels present
    lorecommendations = \
            [np.array([-956.66885744, -379.21022537, -14.29585798,  -1.39085409,           # l = 0
                          3.41674973,   11.85763359,  23.55534883,  38.22251395,
                         55.71839349,   75.96314341,  98.90623003, 124.50690254,
                        152.73632137,  183.57226937, 216.99619737, 252.99363124,
                         291.5523508,  332.66195298, 376.31355739, 422.49933633,
                        471.21242405]),
              np.array([-1.15376230e+02, -1.10985082e+01, -5.10197715e-01, 4.25309290e+00, # l = 1
                        1.26683041e+01, 2.41783693e+01, 3.85543139e+01, 5.56730176e+01,
                        7.54659089e+01, 9.78905522e+01, 1.22912605e+02, 1.50509392e+02,
                        1.80663611e+02, 2.13361294e+02, 2.48592067e+02, 2.86347158e+02,
                        3.26619318e+02, 3.69402477e+02, 4.14691420e+02, 4.62481804e+02,
                        5.12769999e+02]),
              np.array([-5.85379645, 0.90588827, 5.63201191, 13.67427176,                  # l = 2
                        24.54381621, 38.11976977, 54.32810605, 73.13682452,
                        94.52501993, 118.47340661, 144.97035705, 174.00450227,
                        205.56671156, 239.64966278, 276.24676259, 315.35282809,
                        356.96353638, 401.07535931, 447.68537818, 496.79100982,
                        548.38988507]),
              np.array([2.35650373, 6.64673662, 13.67649278, 23.41064726,                  # l = 3
                        35.80378511, 50.81286774, 68.42348545, 88.60939206,
                        111.35263121, 136.64066861, 164.46209382, 194.80989773,
                        227.67761387, 263.06002784, 300.9526912, 341.3513524,
                        384.25219036, 429.65163923, 477.54653223, 527.93416228,
                        580.81227924]),
              np.array([4.77449195, 11.17705627, 20.04560097, 31.43639778,                 # l = 4
                        45.38373733, 61.8787738, 80.91492554, 102.49376441,
                        126.60673075, 153.24771388, 182.40983259, 214.08650332,
                        248.27299967, 284.9649611, 324.15906139, 365.85254928,
                        410.04302461, 456.72834158, 505.90639589, 557.57511102,
                        611.73245043]),
              np.array([7.18991969, 15.3531858, 25.8441029, 38.79380026,                   # l = 5
                        54.23784225, 72.19494694, 92.65832815, 115.6309169,
                        141.11461203, 169.10725897, 199.60765626, 232.61252883,
                        268.11869938, 306.12323706, 346.62326487, 389.61642489,
                        435.10068964, 483.0744059, 533.53621293, 586.48490917,
                        641.91935604]),
              np.array([9.80163424, 19.61275006, 31.61368533, 46.03733346,                 # l = 6
                        62.91812641, 82.28555374, 104.14407508, 128.49176281,
                        155.33294529, 184.66830233, 216.49830788, 250.82304796,
                        287.64120295, 326.95142367, 368.75203372, 413.04128971,
                        459.8175576, 509.07931415, 560.82525336, 615.05428752,
                        671.76553404])]

    # Representative of what one would extract from LINENGY.OUT, for a given atom
    linear_energies = {0: -1.390000000,
                       1: -0.510000000,
                       2:  0.330000000,
                       3:  1.000000000,
                       4:  1.000000000}

    default_lo_nodes = {0: 3,
                        1: 2,
                        2: 0,
                        3: 0,
                        4: 0}

    assert len(linear_energies) == len(default_lo_nodes), \
        "Need consistent l-channels for information extracted from the default basis"

    # hartree
    tol = 0.1
    optimised_basis = filter_default_functions(lorecommendations, linear_energies,
                                               default_lo_nodes, energy_tolerance=tol)

    assert len(optimised_basis) == len(linear_energies), \
        "Need consistent l-channels between the optimised basis and the default basis"

    assert optimised_basis[0].l_value == 0, "l_value consistent with optimised_basis list index"
    assert optimised_basis[0].first_n_nodes == 4, "first radial function for optimised basis has 4 nodes"

    np.testing.assert_approx_equal(optimised_basis[0].energies[0], 3.41674973,
                                   err_msg="first energy parameter in optimised basis should come "
                                           "immediately after -1.39085409 in the reference default basis")

    i = optimised_basis[0].first_n_nodes
    assert optimised_basis[0].energies[0] == lorecommendations[0][i], "some message"

    return