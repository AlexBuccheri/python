import numpy as np

def parse_lorecommendations(file_name:str, n_species:int):
    """

    Note, species symbols not given in the file. Just an index.

    :param file_name: file containing lorecommendations
    :param n_atoms:  Number of species in system
    :return:
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

    # Skip header and first species
    lines = lines[3:]

    i = 0
    basis_per_atom = []
    for i_atom in range(0, n_species):
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
        # species
        i+=1
        basis_per_atom.append(lo_l)

    return basis_per_atom
