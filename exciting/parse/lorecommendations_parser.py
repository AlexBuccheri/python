"""
Tools for parsing lorecommendations

"""
import numpy as np

# Path assumes script being executed from root directory
from exciting_utils.py_grep import grep


def parse_lorecommendations(file_name:str, species:list)-> dict:
    """
    Parse lorecommendations

    If the number of species is not passed, the routine will
    extract it from the file and return all.

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

    :param file_name: file containing lorecommendations
    :param species:  list of species characters
    :return: basis_per_species dictionary
    """

    #if n_species == None:
    #    n_species = int(grep('species', fname=file_name).split()[-1])

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
