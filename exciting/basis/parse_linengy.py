"""
Create basis function strings using LINENGY
"""
from typing import List
import numpy as np

from exciting_utils.py_grep import line_number


#TODO 100% create structured output for LINENGY.OUT (json) for easy parsing

def start_end_indices_per_atom(lines):
    atom_indices = [i for i, string in enumerate(lines) if 'Species' in string]
    start_end = []
    for iatom in range(0, len(atom_indices) - 1):
        start_end.append((atom_indices[iatom], atom_indices[iatom + 1]))
    # Do last (start, end) pair manually
    start_end.append((atom_indices[-1], len(lines)))
    return start_end


def basis_per_atom(lines):
    start_end = start_end_indices_per_atom(lines)
    atoms = []
    for i0, i1 in start_end:
        atoms.append(lines[i0:i1])
    return atoms


def lo_basis_per_atom(atoms: list):
    """

    :param atoms:
    :return:
    """
    lo_block_indices = []
    for atom in atoms:

        for i, basis_func in enumerate(atoms[0]):
            if 'l.o.' in basis_func:
                # Find 1st instance of 'l.o.'
                # Assumes no other basis functions after l.o.'s (per atom)
                lo_block_indices.append((i, len(atom)))
                break

    lo_blocks = []
    for i, atom in enumerate(atoms):
        start, end = lo_block_indices[i]
        lo_blocks.append(atom[start:end])

    return lo_blocks


class LoFunction:
    # Store all derivatives
    # index of lo starting containing radial functions
    # from matchingOrder = 0 to matchingOrder=max_matching_order
    xml_template = '<wf matchingOrder="{mo}" trialEnergy="{te}" searchE="false"/>'

    def __init__(self, index, l, max_file_order, trial_energy):
        self.index = index
        self.l = l
        self.trial_energy = trial_energy
        self.max_file_order = max_file_order
        self.max_matching_order = max_file_order - 1
        # self.file_order = []
        # self.matching_order = []

    def xml_string(self):
        string = '<lo l="{l}">'.format(l=self.l)
        for matching_order in range(0, self.max_matching_order + 1):
            string += self.xml_template.format(**{'mo':matching_order, 'te':self.trial_energy})
        return string + '</lo>'

    # def add_file_order(self, order):
    #     self.file_order.append(order)
    #
    # def matching_order_from_file_order(self, option='linear'):
    #      if option.lower() == 'linear':
    #          # Assume matching order linearly increases, starting at 0
    #          # This should just be printed in LINENGY.OUT
    #          self.matching_order = [i-1 for i in self.file_order]


def parse_lo_function_line(lo_string: str):
    index, l, order_and_energy = lo_string.split(',')
    order, energy = order_and_energy.split(':')
    index = index.split('=')[-1]
    l = l.split('=')[-1]
    order = order.split('=')[-1]
    return {'index': int(index), 'l': int(l), 'order': int(order), 'energy': float(energy)}


def string_LoObject(lo_basis_block_string):
    """
    lo basis for a given atom
    """

    # Remove trailing whitespace entry as it will kill the split
    lo_basis_block_string = lo_basis_block_string[:-1]

    lo_basis_functions = []
    prior_line_data = parse_lo_function_line(lo_basis_block_string[0])

    for lo_string in lo_basis_block_string[1:]:
        line_data = parse_lo_function_line(lo_string)

        # For new lo, store the data from the prior lo
        if line_data['order'] == 1:
            lo_basis_functions.append(
                LoFunction(prior_line_data['index'],
                           prior_line_data['l'],
                           prior_line_data['order'],
                           prior_line_data['energy']
                           ))
        prior_line_data = line_data

    return lo_basis_functions





def local_orbital_basis_string(lo_basis_functions: List[LoFunction], trial_energy_cutoff: float):
    """
    Generate lo basis strings from the data - see Andris's example
    :param lo_function:
    :return:
    """
    xml_string = ''
    for lo in lo_basis_functions:
        if lo.trial_energy <= trial_energy_cutoff:
            lo.xml_string()
            xml_string += lo.xml_string() + '\n'
    return xml_string



# ---------------------
# Main routine
# ---------------------
# Parse linengy.out
fid = open(file='LINENGY.OUT', mode='r')
lines = fid.readlines()
fid.close()

atoms = basis_per_atom(lines)
lo_block_strings = lo_basis_per_atom(atoms)
lo_basis_functions = string_LoObject(lo_block_strings[0])
string = local_orbital_basis_string(lo_basis_functions, 100.)
print(string)





# 6 calculations
# Include Step 20 (ha?) i.e. -ve-0, >0-20, >20-40, >40-60, >60-80, >80-100

# Plot change in self-energy, quasi-particle energy and a couple other things

# Choose one converged to 0.02 eV and plot band structure

# Increase other quantities (and plot band structure)

