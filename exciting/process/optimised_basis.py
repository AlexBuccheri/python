"""

See create_basis.py in zr_lmax6_lmax5 for more functions to adapt

"""
import numpy as np
from typing import List, Dict


class LOEnergies:
    def __init__(self, l_value: int, first_n_nodes: int, energies, last_n_nodes=None):
        """
        :param l_value: Orbital momentum associated with the l-channel
        :param first_n_nodes:  Number of nodes associated with first function in optimised basis
        :param energies: Energy parameter associated with each
        """
        self.l_value = l_value
        self.first_n_nodes = first_n_nodes
        self.last_n_nodes = last_n_nodes
        self.energies = energies


def filter_default_functions(lo_recommendations: List[np.ndarray],
                            linear_energies: dict,
                            default_lo_nodes: dict,
                            energy_tolerance=0.1) -> List[LOEnergies]:
    """

    Every input is w.r.t. one species!

    Given some parsed lorecommendations, systematically exclude lo's that are
    already in the default (ground state) basis.

    The energy parameter of low-l-value functions (found in the default basis)
    should be comparable to that of the final energy parameter in LINENGY.OUT.
    Therefore exclude all functions from lorecommendations <~ energy_parameter
    of the l-channelin LINENGY.OUT.

    For high-l-values that correspond to conduction states, the lowest energy parameter
    in lorecommendations may exceed the energy parameter of the l-channel in the default
    basis. In this case, all lorecommendations with n_nodes <= max_node(lo) in the
    corresponding l-channel of the default basis should be excluded.

    The number of nodes directly correlates with the number of energy parameters for
    a given l-channel. For example if:

      l.o. =  7, l =  3, order =  1 :    1.000000000
      l.o. =  7, l =  3, order =  2 :    1.000000000
      l.o. =  8, l =  3, order =  1 :    1.000000000
      l.o. =  8, l =  3, order =  2 :    1.000000000

    For l-channel = 3, max_node = 0 as there's only one energy parameter for all lo's
    present in the channel. As such, for lorecommendations.dat:

    l= 3
    n=  0  2.35656385982863
    n=  1  6.64679733266455
    n=  2  13.6765531715818
    n=  3  23.4107076815681
    n=  4  35.8038458268313

    the n = 0 term would be excluded, and all others would be retained.

    :param: lo_recommendations lorecommendations for all l-channels of a given species.
                               Elements indexed by l-value: 0 to 6.
                               (l_max = 6 is hard-coded in exciting's lorecommendations).
                               Each value is a np array indexed according to 'n' nodes [0, 21].
    :param: linear_energies    linear energies for all default l-channels of a given species
                               This information would typically come from LINENGY.OUT
                               This ASSUMES one linear energy per l-channel.
                               TODO add more structure to this data to make assumption explicit.
    :param: default_lo_nodes   Largest number of radial nodes associated with the functions
                               of a given l-channel, for a given species.

    :param: energy_tolerance  Tolerance for two energy parameters to be considered the same (Ha).

    :return: List of LOEnergies, with length = n_default_l_channels.
             Each element contains an object with the l-channel, number of nodes for
             the first function in the basis, and the recommended energy parameters for
             functions: indexed from first_n_nodes (n in lorecommendations file) to n=20
             (the hard-coded limit in exciting for the max number of nodes associated
             with a function per l-channel).
    """

    n_l_channels = len(lo_recommendations)
    assert n_l_channels == 7, \
        "l_max hard-coded in exciting to be 6 => 7 l-channels per lorecommendation species"

    # Don't add lo recommendations to l-channels that are not in the default basis
    n_default_l_channels = len(default_lo_nodes)

    assert len(linear_energies) == n_default_l_channels, \
        "linear_energies should have N l-channels == default basis"

    optimised_los = []
    for l_value, lo_energies in enumerate(lo_recommendations[:n_default_l_channels]):
        optimised_los.append(filter_default_function(lo_energies,
                                                    linear_energies[l_value],
                                                    default_lo_nodes[l_value],
                                                    l_value,
                                                    energy_tolerance)
                             )

    return optimised_los


def filter_default_function(lo_energies:np.ndarray,
                           linear_energy:float,
                           max_node: int,
                           l_value:int,
                           energy_tolerance=0.1) -> LOEnergies:
    """

    :param lo_energies:   np.ndarray containing energy parameter recommendations for l-channel l_value,
                          for los with nodes ranging from n=0 to n=20.
    :param linear_energy: Highest linear energy for l-channel 'l_value', in the default basis
    :param l_value:       l-channel
    :return: LOEnergies
    """
    assert lo_energies.size == 21, "n(ode) index should vary from 0 to 20"
    assert l_value >=0, "Angular momentum cannot be less than 0"
    assert energy_tolerance >= 0., "energy tolerance must >= 0 Ha"

    #TODO Would be better to just store indices

    # First criterion. Exclude default basis functions from the optimised lo's
    lo_matches = np.amin(np.abs(lo_energies - linear_energy)) <= energy_tolerance
    if lo_matches:
        i = np.argmin(np.abs(lo_energies - linear_energy))
        return LOEnergies(l_value=l_value,
                          first_n_nodes=i + 1,
                          energies=lo_energies[i + 1:]
                         )

    # Second criterion.
    # If the linear_energy energy of los in a given l-channel doesn't match any energy recommendations
    # OR
    # The lowest lo energy recommend exceeds the energy of the largest default function of a given l-channel,
    # only include radial functions with nodes > max_node(functions) in the l-channel of the default basis.
    else:
        assert max_node >= 0, "Function cannot have a negative number of nodes"
        i = max_node + 1
        assert lo_energies[i] > linear_energy + energy_tolerance, \
            "First recommended energy parameter should exceed any already present in the default basis "
        return LOEnergies(l_value=l_value,
                          first_n_nodes=i,
                          energies=lo_energies[i:]
                          )

def filter_high_energy_functions(lo_recommendations: List[np.ndarray], energy_cutoff: float):
    """

    Nodes are interchangable with trial energy parameters, as basis functions
    are solutions of isolated atoms. So each radial function corresponds
    to an eigenstate of the atom with princial QN, l, n nodes and an energy.

    :return:
    """
    optimised_los = {}
    for l_value, lo_energies in enumerate(lo_recommendations):
        n_index = len(lo_energies[lo_energies < energy_cutoff]) - 1
        optimised_los[l_value] = n_index

    return optimised_los


def filter_high_energy_functions_from_optimised_basis(optimised_basis: List[LOEnergies]):

    for l_value, los in enumerate(optimised_basis):
        lo_energies = optimised_basis[l_value]
        assert lo_energies.l_value == l_value

        first_index = 3
        last_index = 6
        n_index = len([3, 4, 5, 6]) -1 = 3
        last_index = first_index +

        [0, 1, 2,      3, 4, 5, 6]

        lo_energies.


        # Would make more sense to hold all energies and just store the start/end indices
        # n_index = len(lo_energies.energies[lo_energies.energies < energy_cutoff]) - 1
        # lo_energies.last_n_nodes = n_index




def generate_lo_l_string(l:int, energies, max_matching_order:int):
    """

    For a given l-channel, generate XML strings for the additional lo functions that comprise
    the optimised GW basis

    :param l: l-channel
    :param energies:
    :param max_matching_order:
    :return: xml string for exciting species.xml file
    """
    template = """
       <lo l="{l}">
        <wf matchingOrder="{mo1}" trialEnergy="{te}" searchE="false"/> <wf matchingOrder="{mo2}" trialEnergy="{te}" searchE="false"/>
       </lo>
        """
    string = ''
    # Turn this into a routine
    matching_orders = [(i, i + 1) for i in range(0, max_matching_order)]
    for trial_energy in energies:
        for mo1, mo2 in matching_orders:
            string += template.format(l=l, te=round(trial_energy, 2), mo1=mo1, mo2=mo2)

    return string
