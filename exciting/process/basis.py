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


def get_atom_labels(file_name:str):
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



def nodes_per_default_l_channel(file_name:str, filter_duplicate_species=True) -> dict:
    """
    Given some lo's per l-channel, quantify the number of nodes in that l-channel.
    Directly inferred from the energy parameters

    Valid for default and optimised basis sets
    :return:
    """

    #TODO This is parsing the file and should be separated
    fid = open(file_name, "r")
    file = fid.readlines()
    fid.close()

    # Get species and local-orbital line numbers
    output = grep("local-orbital functions", file_name, line_number='').splitlines()
    start_indices = [int(line.split(':')[0]) for line in output]
    n_atoms = len(output)

    atom_labels = get_atom_labels(file_name)

    # ---------------------
    unique_elements = (0, atom_labels[0])
    for i, element in enumerate(atom_labels[1:]):
        if element not in unique_elements:
            unique_elements.append((i, element))

    # --------------------


    # First species index not required
    output = grep("Species", file_name, line_number='').splitlines()
    # -2 returns end index to last lo of the prior species
    end_indices = [int(line.split(':')[0]) - 2 for line in output[1:]]
    # Add up to end of file
    end_indices.append(len(file))

    # Want nodes per l-channel

    # Want for each atom, want max_node for each l-channel
    # Doesn't help
    # Really want this per species, not atom

    prior_l_value = 0
    atom_basis = {}

    for iatom in range(0, n_atoms):
        key = atom_labels[iatom]
        start = start_indices[i]
        stop = end_indices[i]

        for ilo in range(start, stop):
            line = file[ilo].split()
            l_value = int(line[5])
            energy_parameter = float(line[-1])

            if l_value != prior_l_value:
                parameters_per_atom[prior_l_value] = len(set(parameters_per_l_channel))
                prior_l_value = l_value
                parameters_per_l_channel = []

            parameters_per_l_channel.append(energy_parameter)
        atom_basis[key] = parameters_per_atom

    if filter_duplicate_species:
        unique_atoms_basis = {}
        for i, key in unique_elements:
            unique_atoms_basis[key] = atom_basis[i]

    # MAYBE ONLY MAKES SENSE COUNTING NODES FOR ENERGIES ABOVE ZERO

    # Convert parameters per atom and l-channel into max_nodes

    # atom_basis has length = n_atoms in system
    # len(atom_basis['species']) = n_l_channels
    # atom_basis['species'][l_value] = number of unique energy parameters in that l-channel
    # max_node = atom_basis['species'][l_value] - 1

    return atom_basis


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