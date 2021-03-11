"""
Given reference ground state basis files and trial energy parameters from lorecommendations,
create basis files for GW calculations

Zr l-max = 6 and O l-max = 5

Highest state in groundstate/KPOINTS.OUT is ????

"""

import numpy as np
from typing import List
from basis_references import oxygen_groundstate_basis, zr_groundstate_basis
from lorecommendations_parser import parse_lorecommendations


def generate_lo_l_string(l, energies, max_matching_order):
    template = """
       <lo l="{l}">
        <wf matchingOrder="{mo1}" trialEnergy="{te}" searchE="false"/>
        <wf matchingOrder="{mo2}" trialEnergy="{te}" searchE="false"/>
       </lo>
        """
    string = ''
    matching_orders = [(i, i + 1) for i in range(0, max_matching_order)]
    for trial_energy in energies:
        for mo1, mo2 in matching_orders:
            string += template.format(l=l, te=round(trial_energy, 2), mo1=mo1, mo2=mo2)

    return string


def filter_trial_energies(trial_energies: list, energy_cutoff: float):
    """
    Expects a list of np arrays - each containing a set of trial energies
    w.r.t. a given l-value

    :param trial_energies:
    :return:
    """
    filtered_trial_energies = []
    for trial_energys_array in trial_energies:
        filtered_trial_energies.append(trial_energys_array[trial_energys_array < energy_cutoff])

    return filtered_trial_energies


def add_to_groundstate_basis(basis_xml: str, trial_energies: List[np.ndarray],
                             energy_cutoff: float, max_matching_order: int) -> str:
    """
    Given some groundstate basis xml file, add extra lo functions to it.

    TODO Think of a approach than manually adding tags to the files
    basis_xml: Groundstate xml string, with tags for where to insert basis functions
    trial_energies: List of nummpy arrays containing trial energies w.r.t. each l-value of a given species
                    Essentially parsed from lorecommendations

    :return: basis string with additional los added
    """
    filtered_trial_energies = filter_trial_energies(trial_energies,  energy_cutoff)

    # Get lo strings from l=0 up to l_max, where l_max is defined as len(trial_energies) - 1
    lo_strings = {}
    for l_value, trial_energies_array in enumerate(filtered_trial_energies):
        key = 'custom_l' + str(l_value)
        lo_strings[key] = generate_lo_l_string(l_value, trial_energies_array, max_matching_order)

    # Add extra lo's to the ground state basis
    return basis_xml.format(**lo_strings)


def remove_low_n_functions(trial_energies: List[np.ndarray], low_energy_cutoff:list, l_max: int, tol=0.1) -> List[np.ndarray]:

    filtered_trial_energies = []

    for i in range(0, l_max + 1):
        trial_energys_array = trial_energies[i]
        energy_cutoff = low_energy_cutoff[i] + tol
        filtered_trial_energies.append(trial_energys_array[trial_energys_array > energy_cutoff])

    return filtered_trial_energies


def create_basis_file():
    """
    Create basis xml files for Zr and O, using the ground state basis files as references, then adding
    lo's up to some l-max specified in basis_references (== len(zr_trial_energies) - 1) and up to some
    trial energy cut-off defined by energy_cutoffs below.

    """

    zr_l_max = 6
    o_l_max = 5

    zr_low_energy_cutoff = [-1.39,
                            -0.51,
                             0.33,
                             1.00,
                             1.00,
                             1.00,
                             1.00]

    o_low_energy_cutoff = [-0.04,
                            0.1,
                            1.0,
                            1.0,
                            1.0,
                            1.0]


    energy_cutoffs = [100]  # [21., 40., 60., 80., 100]
    max_matching_order = 1

    # Had been doing this manually up to this point
    trial_energies_per_atom = parse_lorecommendations('lorecommendations.txt', n_species=2)
    zr_trial_energies = trial_energies_per_atom[0]
    o_trial_energies = trial_energies_per_atom[1]

    # *Idea* is to remove linear dependencies in the basis (don't add functions in that are already there)
    # although this may not be the correct approach
    zr_trial_energies = remove_low_n_functions(zr_trial_energies, zr_low_energy_cutoff, zr_l_max)
    o_trial_energies = remove_low_n_functions(o_trial_energies, o_low_energy_cutoff, o_l_max)


    for energy_cutoff in energy_cutoffs:
        print("ZR. Trial energy cutoff: " + str(energy_cutoff) + " max_matching_order:" + str(max_matching_order))
        print(add_to_groundstate_basis(zr_groundstate_basis, zr_trial_energies, energy_cutoff, max_matching_order))

        print("OXYGEN. Trial energy cutoff: " + str(energy_cutoff) + " max_matching_order:" + str(max_matching_order))
        print(add_to_groundstate_basis(oxygen_groundstate_basis, o_trial_energies, energy_cutoff, max_matching_order))


create_basis_file()


