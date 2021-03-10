"""
Basis operations
"""
from exciting_utils.py_grep import grep

def get_default_basis():
    """
    Extract the default basis from a ground state file
    :return:
    """
    return


def get_atom_labels(file_name:str) -> list:
    """
    Get a list of atom labels (in the correct order)
    TODO Alex
    Expects LINENGY.OUT

    :return: list of atom labels
    """
    atom_labels = []
    species_lines = grep("Species", file_name).splitlines()
    for line in species_lines:
        symbol = line.split(',')[0][-4:].replace("(", "").replace(")", "")
        atom_labels.append(symbol.strip())

    return atom_labels


def get_unique_atom_labels(atom_labels:list) -> list:

    #unique_atom_labels = [(0, atom_labels[0])]
    unique_atom_labels = [atom_labels[0]]

    for i, element in enumerate(atom_labels[1:]):
        if element not in unique_atom_labels:
            #unique_atom_labels.append((i, element))
            unique_atom_labels.append(element)

    return unique_atom_labels



# CHANGE TO MAKE
# Parsing LINENGY should definitely return this format
#     linear_energies = {0: [-5, -1.390000000],
#                        1: [-0.510000000],
#                        2: [0.330000000,0.330000000],
#                        3: [1.000000000],
#                        4: [1.000000000]}
# Then max_nodes_per_default_l_channel is len(set(linear_energies[l_channel ]))
#
# Would also be useful to distibguish between valence and conduction states.

# MAYBE ONLY MAKES SENSE COUNTING NODES FOR ENERGIES ABOVE ZERO

# Convert parameters per atom and l-channel into max_nodes

# atom_basis has length = n_atoms in system
# len(atom_basis['species']) = n_l_channels
# atom_basis['species'][l_value] = number of unique energy parameters in that l-channel
# max_node = atom_basis['species'][l_value] - 1


def parse_lo_linear_energies(file_name:str, filter_duplicate_species=True) -> dict:
    """

    Return a dictionary of the form:

    linear_energies = {0: [-5.12000000, -1.390000000],
                       1: [-0.51000000, -0.510000000],
                       2: [0.330000000,  0.330000000],
                       3: [1.000000000,  1.000000000],
                       4: [1.000000000,  1.000000000]}

    If filter_duplicate_species is true, only one of each species is included in linear_energies

    Valid for default and optimised basis sets
    :return:
    """

    fid = open(file_name, "r")
    file = fid.readlines()
    fid.close()

    atom_labels = get_atom_labels(file_name)
    n_atoms = len(atom_labels)
    unique_atom_labels = get_unique_atom_labels(atom_labels)

    # Get species and local-orbital line numbers
    output = grep("local-orbital functions", file_name, line_number='').splitlines()
    start_indices = [int(line.split(':')[0]) for line in output]

    # First species index not required
    output = grep("Species", file_name, line_number='').splitlines()
    # -2 moves the end index to the last lo of the prior species
    end_indices = [int(line.split(':')[0]) - 2 for line in output[1:]]
    # Add up to end of file
    end_indices.append(len(file))

    # Parse file
    linear_energies_atoms = {}
    energy_parameter = []
    prior_l_value = 0

    for iatom in range(0, n_atoms):
        atom_label = atom_labels[iatom]
        start = start_indices[iatom]
        stop = end_indices[iatom]
        linear_energies = {}

        for ilo in range(start, stop):
            line = file[ilo].split()
            l_value = int(line[5].replace(",", ""))

            if l_value != prior_l_value:
                linear_energies[prior_l_value] = energy_parameter
                prior_l_value = l_value
                energy_parameter = []

            energy_parameter.append(float(line[-1]))
        linear_energies_atoms[atom_label] = linear_energies

    if filter_duplicate_species:
        linear_energies_species = {}
        for atom_label in unique_atom_labels:
            linear_energies_species[atom_label] = linear_energies_atoms[atom_label]

        return linear_energies_species

    else:
        return linear_energies_atoms




def label_default_basis():
    """
    Label the default basis
    :return:
    """
    return


def label_optimised_basis():
    """
    Label the optimised basis
    :return:
    """
    return