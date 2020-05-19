import json
import subprocess
import sys

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
    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

class EwaldOptions:
    def __init__(self, real_cutoff, recip_cutoff, alpha=None, relative_error=None):
        self.real_cutoff = InputEntry('ewald_real_cutoff', real_cutoff, "bohr")
        self.recip_cutoff = InputEntry('ewald_reciprocal_cutoff', recip_cutoff)
        if alpha != None:
            self.alpha = InputEntry('ewald_alpha', alpha)
        if relative_error != None:
            self.relative_error = InputEntry('ewald_relative_error', relative_error)

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


class TranslationCutoffOptions:
    def __init__(self, h0, overlap, repulsive, unit='bohr'):
        self.h0 = InputEntry('h0', h0, unit)
        self.overlap = InputEntry('overlap', overlap, unit)
        self.repulsive = InputEntry('repulsive', repulsive, unit)

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


class KGridOptions:
    def __init__(self, grid_integers, symmetry_reduction=False):
        self.grid_integers = InputEntry("monkhorst_pack", grid_integers)
        self.symmetry_reduction = InputEntry("symmetry_reduction", symmetry_reduction)

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


def list_to_string(grid_integers):
    assert len(grid_integers) == 3
    return "[" + str(grid_integers[0]) + "," \
               + str(grid_integers[1]) + "," \
               + str(grid_integers[2]) + "]"


def generic_str(value) -> str:
    if isinstance(value, str):
        return value
    elif isinstance(value, int):
        return str(value)
    elif isinstance(value, float):
        return str(value)
    elif isinstance(value, bool):
        return str(value).lower()
    elif isinstance(value, list):
        return list_to_string(value)
    else:
        quit("Have not implemented type ", type(value), "in generic_str")


def options_to_string(*args):
    string = ""
    for options in args:
        for key, entry in dict(options).items():
            string += entry.command + " = " + generic_str(entry.value) + " " + str(entry.unit) + '\n'
    return string



# def generate_entos_command_string(named_result, cutoff_options, ewald_options, temperature, mk_grid, symmetry_reduction):
#     structure_string = \
#         """   structure( fractional=[[Si, 0,    0,    0   ],
#                                   [Si, 0.25, 0.25, 0.25]]
#                          lattice( a = 5.431 angstrom
#                                   bravais = fcc
#                                 )
#                      )
#         """
#


named_result = "si_kgrid"
cutoff_options = TranslationCutoffOptions(h0=40, overlap=40, repulsive=40)
ewald_options = EwaldOptions(real_cutoff=10, recip_cutoff=2, alpha=0.5)
kgrid_options = KGridOptions(grid_integers=[2,2,2], symmetry_reduction=False)
temperature_option = SingleOption(InputEntry(command="temperature", value=0, unit="kelvin"))
string = options_to_string(cutoff_options, ewald_options, kgrid_options, temperature_option)
print(string)




