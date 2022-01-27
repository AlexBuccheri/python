"""

"""
import re

import numpy as np


def parse_qcore_structure(file_name: str) -> dict:
    """ Parse structure.

    Whole routine is some HACKY GARBAGE SHIT
    -  Would be easier to get the positions from cif, using pymatgen

    Each file has the same structure:
    mgo_volume := xtb(
      structure(
          fractional= [['H', 0.1988, 0.26105, 0.14368],
                       ['H', 0.30118, 0.76105, 0.35632],
                        .
                        .
                        .
                       ['O', 0.96018, 0.57252, 0.74478]]

          lattice(
                   lattice_vectors = [[3.5263673508626314, 0.0, -0.6223031990366413], [1.4474643463714529e-15, 9.0009574, 5.511496834585842e-16], [0.0, 0.0, 6.952313179919385]]
                   lattice_vectors_unit = angstrom
                  )

    :param file_string:
    :return:
    """
    with open(file_name) as fid:
        qcore_input_lines = fid.readlines()

    # Find the structure block
    start = 0
    end = 0
    for i, line in enumerate(qcore_input_lines):
        line_str = line.strip()
        match_fractional = re.match("fractional", line_str)
        match_lattice = re.match("lattice", line_str)

        if match_fractional:
            start = i
        if match_lattice:
            end = i
            break

    # Remove fractional from first line
    # fractional= [['H', 0.1988, 0.26105, 0.14368],
    symbol_index = qcore_input_lines[start].find('=')
    qcore_input_lines[start] = qcore_input_lines[start][symbol_index + 3:]

    species = []
    fractional_positions = []
    for i in range(start, end - 1):
        # ['O', 0.96018, 0.57252, 0.74478]
        species_str, x, y, z = qcore_input_lines[i].strip()[1:-2].split(',')
        species.append(species_str.strip("\'"))
        fractional_positions.append([float(r) for r in [x, y, z]])

    # The rest can be extracted with re
    symbol_index = qcore_input_lines[end+1].find('=')
    qcore_input_lines[end + 1] = qcore_input_lines[end+1][symbol_index+1:]

    def remove_brackets(stuff):
        stuff = [r.replace('[', '') for r in stuff]
        stuff = [r.replace(']', '') for r in stuff]
        return stuff

    vectors = qcore_input_lines[end + 1].split(',')
    vectors = remove_brackets(vectors)
    a = [float(v) for v in vectors[0:3]]
    b = [float(v) for v in vectors[3:6]]
    c = [float(v) for v in vectors[6:9]]

    unit = qcore_input_lines[end + 2].split()[-1].strip()

    return {'species': species, 'fractional_positions': fractional_positions,
            'lattice': [a, b, c], 'lattice_vectors_unit': unit}




def parse_qcore_settings():
    k_points = [] # monkhorst_pack
    electronic_temperature = [] #
    return None

def apply_lattice_multipliers():
    return


def generate_tb_input():
    return None


class BinaryRunner():
    # Copy from my ASE implementation
    def __init__(self, a):
        self.a = a


def parse_tb_output():
    return None