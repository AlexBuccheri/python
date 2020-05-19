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
        self.h0 = InputEntry('h0', h0, unit)
        self.overlap = InputEntry('overlap', overlap, unit)
        self.repulsive = InputEntry('repulsive', repulsive, unit)
    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

class KGridOptions:
    def __init__(self, grid_integers, symmetry_reduction=False):
        self.grid_integers = grid_integers
        self.symmetry_reduction = symmetry_reduction

def grid_integers_string(grid_integers):
    assert len(grid_integers) == 3
    return "[" + str(grid_integers[0]) + "," \
               + str(grid_integers[1]) + "," \
               + str(grid_integers[2]) + "]"

def command_string(options):
    string = ""
    for key, entry in dict(options).items():
        string += entry.command + " = " + str(entry.value) + " " + str(entry.unit) + '\n'
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
#     cutoff_options_string = \
#         """    h0_cutoff = """ + str(cutoff_options.h0) + """ """ + cutoff_options.unit + """
#         overlap_cutoff = """ + str(cutoff_options.overlap) + """ """ + cutoff_options.unit + """
#         repulsive_cutoff = """ + str(cutoff_options.repulsive) + """ """ + cutoff_options.unit + """
#     """
#
#     ewald_options_string = \
#         """    ewald_real_cutoff = """ + str(ewald_options.real_cutoff) + """ bohr
#         ewald_reciprocal_cutoff = """ + str(ewald_options.recip_cutoff) + """
#         alpha = """ + str(ewald_options.alpha) + """
#     """
#
#     k_grid_options_string = \
#         """    monkhorst_pack = """ + grid_integers_string(mk_grid) + """
#         symmetry_reduction  = """ + str(symmetry_reduction).lower() + """
#     """
#
#     other_options_string = \
#         """    temperature = """ + str(temperature) + """ kelvin
#     """
#
#     entos_input_string = named_result + " := xtb(\n" + \
#                          structure_string + \
#                          cutoff_options_string + \
#                          ewald_options_string + \
#                          k_grid_options_string + \
#                          other_options_string + "\n)"
#     return entos_input_string
#


# Make these functions more modular

named_result = "si_kgrid"
cutoff_options = TranslationCutoffOptions(h0=40, overlap=40, repulsive=40)

string = command_string(cutoff_options)
print(string)

ewald_options = EwaldOptions(real_cutoff=10, recip_cutoff=2, alpha=0.5)
temperature = 0
mk_grid = [2,2,2]
kgrid_options = KGridOptions(grid_integers=mk_grid, symmetry_reduction=False)

# Really want to make a list from the class, then iteate over to build "command = " + "value " + "unit"

# entos string can just take **awgs
# for every other string, define a string function to return it
#Want to abstract this stuff so it's usable for each input

# entos_input_string = generate_entos_command_string(named_result,
#    cutoff_options, ewald_options, temperature, mk_grid, symmetry_reduction)
#



