"""

See create_basis.py in zr_lmax6_lmax5 for more functions to adapt

"""
import numpy as np
from typing import List, Dict, Optional
from copy import deepcopy
import re


# TODO Change name to LOOptimisedEnergies
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
        #TODO Change name to energy_recommendations
        self.energies = energies
        self.first_node = first_node if first_node is not None else 0
        self.last_node = last_node if last_node is not None else energies.size - 1

    def get_optimised_energies(self):
        return self.energies[self.first_node: self.last_node+1]


# TODO change name to DefaultLOEnergies
class DefaultLOs():

    def get_max_nodes(self) -> dict:
        """
        TODO MAYBE ONLY MAKES SENSE COUNTING NODES FOR ENERGIES ABOVE ZERO?
        Would also be useful to distinguish between valence and conduction states.

        Get the number of nodes associated with the lo function with max nodes, per l-channel.

        Nodes appear to be interchangeable with the energy parameters, as basis functions
        are solutions of isolated atoms. So each radial function corresponds
        to an eigenstate of the atom with principal QN, l, n nodes and an energy.

        As such, max_nodes = number of unique trial energies - 1.

        Expect linear_energies of the form
            {0: [-5.12000000, -1.390000000],
             1: [-0.51000000, -0.510000000],
             2: [0.330000000,  0.330000000],
             3: [1.000000000,  1.000000000],
             4: [1.000000000,  1.000000000]}

        self.nodes = {0: 1
                      1: 0
                      2: 0
                      3: 0
                      4: 0}

        where set(energies) removes duplicates.

        :return: dictionary of max_nodes
        """
        self.nodes = {}
        for l_value, energies in self.linear_energies.items():
            assert type(energies) == list
            self.nodes[l_value] = len(set(energies)) - 1
        return self.nodes

    def __init__(self, linear_energies: dict, nodes: Optional[dict] = None, energy_tol=0.1):
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


def filter_lo_functions(lo_recommendations: List[np.ndarray],
                        default_los: Optional[DefaultLOs] = None,
                        optimised_lo_cutoff: Optional[dict] = None) -> List[LOEnergies]:
    """

    Every input is w.r.t. one species.

    Given some parsed lorecommendations, systematically exclude lo's that are
    already in the default (ground state) basis and/or exceed the energy cutoff
    for the optimised lo basis.

    :param: lo_recommendations lorecommendations for all l-channels of a given species.
                               Elements indexed by l-value: 0 to 6.
                               (l_max = 6 is hard-coded in exciting's lorecommendations).
                               Each value is a np array indexed according to 'n' nodes [0, 20].
    :param: default_los    linear energies for all default l-channels of a given species
                               This information would typically come from LINENGY.OUT
                               This ASSUMES one linear energy per l-channel.
    :param: optimised_lo_cutoff dict containing LO energy cut-off per l-channel:
                                of the form  {0: 100., 1: 100., 2: 200., 3: 80.}


    :return: List of LOEnergies, with length = n_default_l_channels.
             Each element contains an object with the l-channel, number of nodes for
             the first function in the basis, and the recommended energy parameters for
             functions: indexed from first_n_nodes (n in lorecommendations file) to n=20
             (the hard-coded limit in exciting for the max number of nodes associated
             with a function per l-channel).
    """

    assert len(lo_recommendations) == 8, "Expect 8 l-channels for lorecommendations, " \
                                         "per species, as it is hard-coded in exciting"

    assert len(lo_recommendations[0]) == 21, "Expect 21 entries per l-channel lorecommendations, " \
                                             "as it is hard-coded in exciting"

    # Don't evaluate lo recommendations for l-channels that are not in the default basis
    n_default_l_channels = len(default_los.linear_energies)

    optimised_los = []
    for l_value, lo_energies in enumerate(lo_recommendations[:n_default_l_channels]):

        optimised_los_lchannel = LOEnergies(l_value, lo_energies)

        if default_los is not None:
            optimised_los_lchannel = filter_default_functions(optimised_los_lchannel,
                                                             max(default_los.linear_energies[l_value]),
                                                             default_los.nodes[l_value],
                                                             default_los.energy_tol
                                                              )
        if optimised_lo_cutoff is not None:
            optimised_los_lchannel = filter_high_energy_optimised_functions(optimised_los_lchannel,
                                                                           optimised_lo_cutoff[l_value])

        optimised_los.append(deepcopy(optimised_los_lchannel))
        del optimised_los_lchannel

    return optimised_los


def filter_default_functions(optimised_los:LOEnergies,
                             default_max_linear_energy: float,
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
    :param default_max_linear_energy: Highest linear energy for l-channel 'l_value', in the default basis.
    :param default_max_node: lo function with the most nodes for an l-channel of the default basis.
    :return: LOEnergies with first_node set
    """
    assert optimised_los.energies.size == 21, "n(ode) index should vary from 0 to 20"
    assert optimised_los.l_value >=0, "Angular momentum cannot be less than 0"
    assert energy_tolerance >= 0., "energy tolerance must >= 0 Ha"

    lo_matches = np.amin(np.abs(optimised_los.energies - default_max_linear_energy)) <= energy_tolerance
    if lo_matches:
        i = np.argmin(np.abs(optimised_los.energies - default_max_linear_energy))
        optimised_los.first_node = i + 1

    else:
        assert default_max_node >= 0, "Function cannot have a negative number of nodes"
        i = default_max_node + 1
        if optimised_los.energies[i] <= default_max_linear_energy + energy_tolerance:
            print("First recommended energy parameter should exceed the energies of all default los")
            print("First optimised lo energy:", optimised_los.energies[i])
            print("Highest default lo energy, tolerance (Ha):", default_max_linear_energy, energy_tolerance)
            quit()

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
        <wf matchingOrder="{mo1}" trialEnergy="{te}" searchE="false"/> 
        <wf matchingOrder="{mo2}" trialEnergy="{te}" searchE="false"/>
       </lo>
        """
    string = ''
    matching_orders = generate_matching_orders(max_matching_order)
    for trial_energy in energies:
        for mo1, mo2 in matching_orders:
            string += template.format(l=l, te=round(trial_energy, 2), mo1=mo1, mo2=mo2)

    return string



def generate_optimised_basis_string(ground_state_xml: str,
                                    optmised_basis_recommendations: List[LOEnergies],
                                    max_matching_order: int) -> str:
    """
    Given some groundstate basis string (xml format), add optimised lo functions to it.

    :param ground_state_xml:
    :param optmised_basis_recommendations:
    :param max_matching_order:
    :return optimised basis string
    """
    lo_strings = {}
    for los in optmised_basis_recommendations:
        key = 'custom_l' + str(los.l_value)
        lo_strings[key] = generate_lo_l_string(los.l_value, los.get_optimised_energies(), max_matching_order)

    # Add extra lo's to the ground state basis
    return ground_state_xml.format(**lo_strings)


def parse_species_string(l_channels: list, basis_string: str):
    """
    Parse LO functions from a species.xml string.

    Every  <lo l="l-channel">  ... </lo> is counted as a unique LO in with
    orbital momentum "l-channel".

    Returns a dictionary of LO functions, of the form:
                     {0: [lo_0, lo_1, ..., lo_i],
                      1: [lo_0, lo_1, ..., lo_i],
                      .
                      .
                      .
                      l_max: [lo_0, lo_1, ..., lo_i]
                      }
    where
          lo_i = [radial_0, radial_1] and
          radial_i = {'matching_order', 'trial_energy', 'searchE'}

    :params l_channels: List of l-channels in the basis i.e. [0, 1, 3]
                        i.e. it allows l-channels to be ignored by omission.
    :param basis_string:  Expects a string from species.xml

    :return basis_los: Dictionary described above.
    """
    basis_string = basis_string.splitlines()
    basis_los = {l: [] for l in l_channels}
    lo_block = False

    for line in basis_string:

        # Start of lo block description
        if "<lo" in line:
            lo = []
            l = int(re.findall('"([^"]*)"', line)[0])

            # LO description on same line as LO tag
            if "<wf" in line:
                radial_params = re.findall('"([^"]*)"', line)[1:]
                n_radial_functions =  len(radial_params) / 3.
                assert n_radial_functions.is_integer(), 'Each radial function should be defined by 3 parameters'

                i = 0
                for ir in range(0, int(n_radial_functions)):
                    radial_function = \
                        {'matching_order': int(radial_params[i]),
                         'trial_energy': float(radial_params[i + 1]),
                         'searchE': radial_params[i + 2] == "true"}
                    i += 3
                    lo.append(radial_function)

        # LO description on separate line to the LO tag
        if ("<wf" in line) and ("<lo" not in line):
            radial_params = re.findall('"([^"]*)"', line)
            n_radial_functions = len(radial_params) / 3.
            assert n_radial_functions.is_integer(), 'Each radial function should be defined by 3 parameters'

            i = 0
            for ir in range(0, int(n_radial_functions)):
                radial_function = \
                    {'matching_order': int(radial_params[i]),
                     'trial_energy': float(radial_params[i + 1]),
                     'searchE': radial_params[i + 2] == "true"}
                lo.append(radial_function)
                i += 3

        if "</lo" in line:
            basis_los[l].append(lo)
            lo = []

    return basis_los


def restructure_energy_cutoffs(energy_cutoffs: dict) -> list:
    """
    Get energy cut-offs in a more useful structure to iterate over

    :param energy_cutoffs: dict of the form

        {'zr': {0: np.linspace(60, 120, num=4),
                1: np.linspace(60, 120, num=4),
                2: np.linspace(90, 300, num=4),
                3: np.linspace(60, 120, num=4)},
         'o':  {0: np.linspace(60, 120, num=4),
                1: np.linspace(60, 120, num=4),
                2: np.linspace(60, 120, num=4)}
         }

    which is convenient to manually define. Note, num = N must be the same N for all energy cut-offs
    else this routine will throw an error.

    :return restructured_energies: list of the form

    [ {species: {0: energy_0,  {species: {0: energy_1,    ...   {species: {0: energy_n,  ]
                 1: energy_0,             1: energy_1,                     1: energy_n,
                 2: energy_0,             2: energy_1,                     2: energy_n,
                 l: energy_0},            l: energy_1},                    l: energy_n}

    where len(restructured_energies) = n_energies_per_channel
    """

    # Assumes energy_cutoffs is small enough that one can manually inspect, if the assertion fails
    n_energies = []
    for species, l_channels in energy_cutoffs.items():
        for energies in l_channels.values():
            n_energies.append(len(energies))

    unique_n_energies = set(n_energies)
    assert len(unique_n_energies) == 1, "Number of energy cutoffs not the same for all (species, l-channels) " \
                                        "Routine 'restructure_energy_cutoffs' will not work"
    n_energies_per_channel = list(unique_n_energies)[0]

    restructured_energies = []
    for inum in range(0, n_energies_per_channel):
        data = {}
        for species, l_channels in  energy_cutoffs.items():
            data[species] = {l:energies[inum] for l, energies in l_channels.items()}
        restructured_energies.append(data)

    return restructured_energies



def create_lo_label(basis_los: dict) -> list:
    """
    Create labels for the LO basis defined in a species xml file.

    :param basis_los: Dictionary of LOs, constructed by parsing a species.xml string
    :return: basis_labels of the form ['6s ', '3d '], for example
    """
    # TODO Check if 7 = j
    channel_label = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i', 7: 'j'}
    basis_labels = []

    for l, los in basis_los.items():
        basis_labels.append(str(len(los)) + channel_label[l] + ' ')

    return basis_labels


def latex_lo_basis_labels(basis_los: dict, lorecommendations: List[np.ndarray]):
    """
    Having parsed a basis dictionary for a given species, generate latex expressions
    for each LO basis function.

    Label the principal quantum number (pqn) according to where the energy parameter
    is in the lorecommendations. Note that n in lorecommendations corresponds to
    the number of nodes, therefore pqn = n_nodes + 1

    :param basis_los: Dictionary of LO functions, of the form:
                     {l: los}, los = [ lo_0, lo_1, ..., lo_i],
                     lo_i = [radial_0, radial_1] and
                     radial_i = {'matching_order', 'trial_energy', 'searchE'}

    :param lorecommendations: List[np.ndarray]  = [species_l0, species_l1, species_l2, ...]
    where  species_li = energy parameters for a l-channel i, from lorecommendations.

    :return: Prints to stdout
    """
    # Orbital channel labels
    channel_labels = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h', 6: 'i'}
    # Radial functions
    matching_order_symbols = {0: "u", 1:"\\dot{u}", 2:"\\ddot{u}"}
    # Coefficients for linear sum of radial functions
    coefficient = {0: 'a', 1: 'b'}
    # High because LO recommendations may not match that closely to trial parameters
    tol = 1.

    for l, los in basis_los.items():
        lo_recommendations_l_channel = lorecommendations[l]
        l_label = channel_labels[l]
        print(' ------------------')
        print(" l = " +str(l))
        print(' ------------------')

        for lo in los:
            string = ''
            for ir, radial_fun in enumerate(lo):
                mo, energy, = radial_fun['matching_order'], radial_fun['trial_energy']
                n_nodes = np.where(np.abs(lo_recommendations_l_channel - energy) < tol)[0]
                assert n_nodes.size, "Matched more than one energy in LO recommendation of this l-channel"
                pqn_str = str(n_nodes[0] + 1)
                string += coefficient[ir] + matching_order_symbols[mo] + '(r;\\epsilon_{' + pqn_str + l_label + '}) + '
            string = string[:-3]
            print(string)

    return


def table_of_lo_energies(basis_los: dict, lorecommendations: List[np.ndarray]):
    """
    For a given species, output tabulated energy parameter data in the form:
        eps_nl   trial_energy (Ha)
         ε3s        −4.5
         ε4s         0.0
         ε5s        29.1
         ε6s        49.6
         ε7s        74.5

    :param basis_los: Dictionary of LO functions, of the form:
                     {l: los}, los = [ lo_0, lo_1, ..., lo_i],
                     lo_i = [radial_0, radial_1] and
                     radial_i = {'matching_order', 'trial_energy', 'searchE'}

    :param lorecommendations: List[np.ndarray]  = [species_l0, species_l1, species_l2, ...]
    where  species_li = energy parameters for a l-channel i, from lorecommendations.

    :return: Prints to stdout
    """

    # Orbital channel labels
    channel_labels = {0: 's', 1: 'p', 2: 'd', 3: 'f'}

    # High because LO recommendations may not match that closely to trial parameters
    tol = 1.

    print(' # eps_nl     # Trial energy (Ha)')
    for l, los in basis_los.items():
        lo_recommendations_l_channel = lorecommendations[l]
        l_label = channel_labels[l]

        for lo in los:
            trial_energies = [radial_function['trial_energy'] for radial_function in lo]
            trial_energy = list(set(trial_energies))
            assert len(trial_energy) == 1, "Should be one trial energy for a set of radial functions"

            n_nodes = np.where(np.abs(lo_recommendations_l_channel - trial_energy[0]) < tol)[0]
            assert n_nodes.size, "Matched more than one energy in LO recommendation for this l-channel"
            pqn_str = str(n_nodes[0] + 1)
            print('\\epsilon_{' + pqn_str + l_label + '}', trial_energy[0])

    return
