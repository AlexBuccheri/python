

class InputEntry:
    def __init__(self, command, value, unit=None):
        self.command = command
        self.value = value
        if unit !=None:
            self.unit = unit
        if unit == None:
            self.unit = ""


class SingleOption:
    def __init__(self, option):
        self.option = option


class EwaldOptions:
    def __init__(self, real_cutoff, recip_cutoff, alpha=None, relative_error=None):
        self.real_cutoff = InputEntry('ewald_real_cutoff', real_cutoff, "bohr")
        self.recip_cutoff = InputEntry('ewald_reciprocal_cutoff', recip_cutoff)
        if alpha != None:
            self.alpha = InputEntry('ewald_alpha', alpha)
        if relative_error != None:
            self.relative_error = InputEntry('ewald_relative_error', relative_error)


class TranslationCutoffOptions:
    def __init__(self, h0, overlap, repulsive, unit='bohr'):
        self.h0 = InputEntry('h0_cutoff', h0, unit)
        self.overlap = InputEntry('overlap_cutoff', overlap, unit)
        self.repulsive = InputEntry('repulsive_cutoff', repulsive, unit)


class KGridOptions:
    def __init__(self, grid_integers, symmetry_reduction=False):
        self.grid_integers = InputEntry("monkhorst_pack", grid_integers)
        self.symmetry_reduction = InputEntry("symmetry_reduction", symmetry_reduction)


class LatticeOpt:
    def __init__(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None,
                 length_unit = 'angstrom', angle_unit = 'degree'):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.length_unit = length_unit
        self.angle_unit = angle_unit


def generate_structure_string(unit: str, lattice: LatticeOpt, atoms: list, bravais: str):

    #For formatting
    string = "structure( " + unit + "=["
    atom = atoms[0]
    string += "[" + atom[0] + "," + list_to_string(atom[1:])[1:] + ',\n'
    white_space = " " * 11

    for atom in atoms[1:]:
        assert len(atom) == 4
        string += white_space + "[" + atom[0] + "," + list_to_string(atom[1:])[1:] +',\n'
    string = string[:-2] + '] \n'

    lattice_opts = lattice.__dict__
    angle_unit = lattice_opts.pop('angle_unit')
    length_unit = lattice_opts.pop('length_unit')

    string +=  "   lattice( \n"
    for key, entry in lattice_opts.items():
        # Not the cleanest but doesn't matter
        if len(key) == 1:
            unit = length_unit
        else:
            unit = angle_unit
        if entry != None:
            string += white_space + key + " = " + str(entry) + " " + unit +'\n'

    string += white_space + "bravais = " + bravais +'\n'
    string += white_space + ")\n"
    string += "        )\n"

    return string


def list_to_string(grid_integers):
    assert len(grid_integers) == 3
    return "[" + str(grid_integers[0]) + "," \
           + str(grid_integers[1]) + "," \
           + str(grid_integers[2]) + "]"


def generic_str(value) -> str:
    if isinstance(value, str):
        return value
    # bool can also be evaluated as int, hence bool comes first
    elif isinstance(value, bool):
        return str(value).lower()
    elif isinstance(value, int) or isinstance(value, float):
        return str(value)
    elif isinstance(value, list):
        return list_to_string(value)
    else:
        quit("Have not implemented type ", type(value), "in generic_str")


def options_to_string(*args):
    string = ""
    for options in args:
        for key, entry in options.__dict__.items():
            string += entry.command + " = " + generic_str(entry.value) + " " + str(entry.unit) + '\n'
    return string


def generate_xtb_string(named_result, *args):
    xtb_string = named_result + " := xtb(\n"
    for string in args:
        xtb_string += string
    xtb_string += ")"
    return xtb_string