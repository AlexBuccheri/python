
# Write structure in xyz format
# Input = list of atoms.Atom i.e. [atoms.Atom(species,pos), atoms.Atom(species,pos), ...]
def xyz(molecule, header = None):
    hdr = ''
    if header != None:
        hdr += header

    string = str(len(molecule)) + '\n'
    string += hdr + '\n'

    for atom in molecule:
        string += atom.species + " " + "".join(str(atom.position))[1:-1] + '\n'

    return string


# Also accepted by periodic xTB
# TODO(Alex) Be able to pass lattice vectors instead of lattice_opts
# TODO(Alex) Pass cell positions in fractional coordinates
def turbomole_riper_periodic(supercell, lattice_opts):

    # $periodic
    n_peroidic_dimensions = 3
    riper_string = "$periodic " + str(n_peroidic_dimensions) + "\n"

    # $cell
    label = 'angs'
    if lattice_opts.units != 'angstrom':
        label = lattice_opts.units
    riper_string += "$cell " + label + "\n"

    lattice_opts = lattice_opts.__dict__
    keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    for key in keys:
        if lattice_opts[key]:
            riper_string += str(lattice_opts[key]) + "  "
        else:
            riper_string += " 0  "
    riper_string += "\n"

    # $coord
    riper_string += "$coord angs \n"
    for atom in supercell:
        riper_string += "".join(str(atom.position))[1:-1] + "  " + atom.species +"\n"
    riper_string += "$end"

    return riper_string