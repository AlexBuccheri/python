"""
Tools for parsing lorecommendations

"""
import numpy as np

# Path assumes script being executed from root directory
from exciting_utils.py_grep import grep


def parse_lorecommendations(file_name:str, species:list)-> dict:
    """
    Parse lorecommendations

    Note: species symbols not given in the file. Just an index.
          'n' is number of nodes, NOT the principal QN
          The energy parameters are slightly inconsistent with those returned by LINENGY.OUT
          because in one file, the basis is defined as |R|^2, and in the other as |rR|^2.

    Return list has the form:
        n_basis_species = len(basis_per_species)

        n_l_channels = len(basis_per_species[0])

        l_channel = basis_per_species[0][l_value]

        where an l_channel contains radial_solutions of an isolated atom, characterised by
        l, n(odes), energy_parameter
        and has n_nodes = len(l_channel) entries.

    basis_per_species = {'species_label1': lo_l,
                         'species_label2': lo_l
                        }

    where lo_l is a list of the form [lo_0_energies, lo_1_energies, lo_2_energies,
                                      lo_3_energies, lo_4_energies, lo_5_energies, lo_6_energies]
    and

    lo_l_energies = [-956.668857441502, -379.210225370091, -14.2958579803164, -1.39085408705473,
                     3.41674973498303,  11.8576335922056,  23.5553488286426,  38.2225139492090,  55.7183934928437,
                     75.9631434089484,  98.9062300331434,  124.506902544624,  152.736321370444,  183.572269369556,
                     216.996197372528,  252.993631242244,  291.552350804867,  332.661952981567,  376.313557393327,
                     422.499336330017,  471.212424048639]

    :param file_name: file containing lorecommendations
    :param species:  list of species characters
    :return: basis_per_species dictionary
    """

    fid = open(file=file_name, mode='r')
    lines = fid.readlines()
    fid.close()

    # Hard-coded in exciting
    l_min = 0
    l_max = 6
    n_min = 0
    n_max = 20

    lo_trial_energy = np.empty(shape=(n_max+1))

    # Skip header and first species index
    lines = lines[3:]

    i = 0
    basis_per_species = {}
    for i_atom in range(0, len(species)):
        lo_l = []
        for i_l in range(l_min, l_max + 1):
            # l_index line
            i+=1
            for i_n in range(n_min, n_max + 1):
                lo_trial_energy[i_n] = float(lines[i].split()[2])
                i+=1
            # Single line break
            lo_l.append(np.copy(lo_trial_energy))
            i+=1
        # skip species index line
        i+=1
        key = species[i_atom]
        basis_per_species[key] = lo_l

    return basis_per_species
