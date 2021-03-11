"""

See create_basis.py in zr_lmax6_lmax5 for more functions to adapt

"""
import numpy as np
from typing import List, Dict, Optional
from copy import deepcopy


class LOEnergies:
    def __init__(self, l_value: int, energies: np.ndarray, first_node=None, last_node=None):
        """
        :param l_value: Orbital momentum associated with the l-channel
        :param energies: Energy parameters associated with lo recommendations, for an l-channel.
                         Should range from n = 0 to n = 20 (hard-coded in exciting)

        :param first_node:  Number of nodes associated with first function in optimised basis
                           == Index in energies
        :param last_node:  Number of nodes associated with last function in optimised basis
                           == Index in energies
        """
        assert energies.size == 21
        self.l_value = l_value
        self.energies = energies
        self.first_node = first_node if first_node is not None else 0
        self.last_node = last_node if last_node is not None else energies.size - 1

    def get_optimised_energies(self):
        return self.energies[self.first_node: self.last_node+1]


class DefaultLOs():

    def get_max_nodes(self) -> dict:
        """
        TODO MAYBE ONLY MAKES SENSE COUNTING NODES FOR ENERGIES ABOVE ZERO?
        Would also be useful to distibguish between valence and conduction states.

        Get the number of nodes associated with the lo function with max nodes, per l-channel.

        Nodes appear to be interchangable with the energy parameters, as basis functions
        are solutions of isolated atoms. So each radial function corresponds
        to an eigenstate of the atom with princial QN, l, n nodes and an energy.

        As such, max_nodes = number of unique trial energies - 1.

        Expect linear_energies of the form
            {0: [-5.12000000, -1.390000000],
             1: [-0.51000000, -0.510000000],
             2: [0.330000000,  0.330000000],
             3: [1.000000000,  1.000000000],
             4: [1.000000000,  1.000000000]}

        :return: dictionary of max_nodes
        """
        self.nodes = {}
        for l_value, energies in self.linear_energies.items():
            assert type(energies) == list
            self.nodes[l_value] = len(set(energies)) - 1
        return self.nodes


    def __init__(self, linear_energies:Optional[dict], nodes:Optional[dict]=None, energy_tol=0.1):
        """
        Initialise class
        :param linear_energies: Dictionary of all linear energies for an atom's local orbitals
                                key = l-channel
        :param nodes: Dictionary of the number of nodes associated with the lo with the most nodes,
                      per l-channel of an atom.
        :param energy_tol: energy tolerance for two lo's energy parameters to be considered equivalent.
                           One will see a slight mismatch between lo recommendations and those of the
                           final default basis, output in LINENGY.OUT
        """
        self.linear_energies = linear_energies
        self.nodes = nodes if nodes is not None else self.get_max_nodes()
        self.energy_tol = energy_tol
        assert len(self.linear_energies) == len(self.nodes), \
            "N l-channels in default linear_energies should equal N l-channels in default max nodes"




def filter_default_functions(lo_recommendations: List[np.ndarray],
                             default_los: Optional[DefaultLOs]=None,
                             optimised_lo_cutoff: Optional[list]=None) -> List[LOEnergies]:
    """

    Every input is w.r.t. one species.

    Given some parsed lorecommendations, systematically exclude lo's that are
    already in the default (ground state) basis and/or exceed the energy cutoff
    for the optimised lo basis.

    :param: lo_recommendations lorecommendations for all l-channels of a given species.
                               Elements indexed by l-value: 0 to 6.
                               (l_max = 6 is hard-coded in exciting's lorecommendations).
                               Each value is a np array indexed according to 'n' nodes [0, 21].
    :param: default_los    linear energies for all default l-channels of a given species
                               This information would typically come from LINENGY.OUT
                               This ASSUMES one linear energy per l-channel.
    :param: optimised_lo_cutoff


    :return: List of LOEnergies, with length = n_default_l_channels.
             Each element contains an object with the l-channel, number of nodes for
             the first function in the basis, and the recommended energy parameters for
             functions: indexed from first_n_nodes (n in lorecommendations file) to n=20
             (the hard-coded limit in exciting for the max number of nodes associated
             with a function per l-channel).
    """

    assert len(lo_recommendations) == 7, "Expect 7 l-channels for lorecommendations, " \
                                         "per species, as it is hard-coded in exciting"

    # Don't evaluate lo recommendations for l-channels that are not in the default basis
    n_default_l_channels = len(default_los.linear_energies)

    optimised_los = []
    for l_value, lo_energies in enumerate(lo_recommendations[:n_default_l_channels]):

        optimised_los_lchannel = LOEnergies(l_value, lo_energies)

        if default_los is not None:
            optimised_los_lchannel = filter_default_function(optimised_los_lchannel,
                                                            max(default_los.linear_energies[l_value]),
                                                            max(default_los.nodes[l_value]),
                                                            default_los.energy_tol
                                                            )
        if optimised_lo_cutoff is not None:
            optimised_los_lchannel = filter_high_energy_optimised_functions(optimised_los_lchannel,
                                                                            optimised_lo_cutoff[l_value])

        optimised_los.append(deepcopy(optimised_los_lchannel))
        del optimised_los_lchannel

    return optimised_los


def filter_default_function(optimised_los:LOEnergies,
                            default_linear_energy: float,
                            default_max_node: int,
                            energy_tolerance: float) -> LOEnergies:
    """
    Find an index for lo recommendations that screens all contributions already present in the default basis.


The energy parameter of low-l-value functions (found in the default basis)
    should be comparable to that of the final energy parameter in LINENGY.OUT.
    Therefore exclude all functions from lorecommendations <= energy_parameter
    of the l-channel in LINENGY.OUT.

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





    First criterion. Match an lo recommendation with the max energy default lo in the same l-channel.
                     Exclude that and all lower-energy default basis functions from the optimised lo's.

    Second criterion.
        If the max default_linear_energy of the default los in a given l-channel doesn't match any
        energy recommendations, evaluate according to the number of nodes.

        In which case, only include radial functions with nodes > max_node(functions) in the l-channel
        of the default basis.

        This will always occur if the lowest energy recommendation exceeds ~ 1 Ha, implying it's an unbound solution
        of the radial Schrodinger equation. The largest energy parameter in valence los of the default basis
        typical doesn't exceed 0.15-1 Ha.

        This will sometimes occur if one has recommendations ~ [-5.85379645, 0.90588827, ...] but the max energy
        associated with the default los is inbetween the two i.e. ~0.3.
        In general lo recommendations become more inconsistent with the default basis as the energy >= 0 Ha.
        This is needs resolving.

    :param optimised_los: LOEnergies initialised with l-channel and lo recommendation energies.
    :param default_linear_energy: Highest linear energy for l-channel 'l_value', in the default basis.
    :param default_max_node: lo function with the most nodes for an l-channel of the default basis.
    :return: LOEnergies with first_node set
    """
    assert optimised_los.energies.size == 21, "n(ode) index should vary from 0 to 20"
    assert optimised_los.l_value >=0, "Angular momentum cannot be less than 0"
    assert energy_tolerance >= 0., "energy tolerance must >= 0 Ha"

    lo_matches = np.amin(np.abs(optimised_los.energies - default_linear_energy)) <= energy_tolerance
    if lo_matches:
        i = np.argmin(np.abs(optimised_los.energies - default_linear_energy))
        optimised_los.first_node = i + 1

    else:
        assert default_max_node >= 0, "Function cannot have a negative number of nodes"
        i = default_max_node + 1
        assert optimised_los.energies[i] > default_linear_energy + energy_tolerance, \
            "First recommended energy parameter should exceed any already present in the default basis "
        optimised_los.first_node = i

    return optimised_los


def filter_high_energy_optimised_functions(optimised_basis: LOEnergies,
                                           energy_cutoff: float):
    """
    Find the index the final lo (for a given l-channel) with an energy parameter less than energy_cutoff.
    This allows one to remove local orbitals defined with high energy parameters, from the basis.

    :param optimised_basis: LOEnergies initialised with l-channel and lo recommendation energies
    :return: LOEnergies with last_node set
    """
    assert type(optimised_basis.energies) == np.ndarray, "optimised_basis.energies has wrong type"
    last_node = len(optimised_basis.energies[optimised_basis.energies < energy_cutoff]) - 1
    optimised_basis.last_node = last_node
    return optimised_basis


def generate_matching_orders(max_matching_order: int):
    """
    Document me
    :param max_matching_order:
    :return:
    """
    return [(i, i + 1) for i in range(0, max_matching_order)]


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
    matching_orders = generate_matching_orders(max_matching_order)
    for trial_energy in energies:
        for mo1, mo2 in matching_orders:
            string += template.format(l=l, te=round(trial_energy, 2), mo1=mo1, mo2=mo2)

    return string
