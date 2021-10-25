"""
Parse various GW output files
"""
import subprocess
import warnings
import numpy as np
from pathlib import Path
import os

from exciting_utils.py_grep import grep


def parse_kpoints(file_path:str, file_name='KPOINTS.OUT') -> dict:
    """
    TODO(Alex) Move this routine
    Parse KPOINTS.OUT

    Looks like the grid uses symmetry reduction
    :return: number of empty states per k-point
    """
    data = np.loadtxt(file_path + '/' + file_name, skiprows = 1)
    n_empty = data[:, 5].astype(int)
    return {'n_empty':n_empty}


def get_nempty_from_evalqp(file_path: str, file_name='EVALQP.DAT') -> int:
    """
    exciting GW uses same number of nempty per k-point

    If nempty is specified as max in input, it should use the lowest
    available of all k-points (noting that nempty differs per k-point
    because the PW cut-off varies per k-point)

    :param file_path: file path
    :param file_name: file name to parse
    :return: nempty, number of empty states per k-point
    """
    lines = grep('k-point', file_path + '/' + file_name, n_lines_before=2)
    return int(lines.splitlines()[2].split()[0])


def parse_gw_evalqp(file_path: str, file_name='EVALQP.DAT') -> dict:
    """
    Parse GW output file EVALQP.DAT

    Repeating structure:
       kpoint k1 k2 l2 weight
       header line
       1
       ...
       nempty (from GW input)
       whitespace

    where Gamma always appears to be the first k-point in the file.

    file_path:
    file_name:
    nkpts: Number of irreducible q-points, I assume

    :return: dictionary of form {ik: k-point, results},
    where results[istate].keys = ['E_KS', 'E_HF', 'E_GW', 'sigma_x',' Re_sigma_c', 'Im_sigma_c', 'V_xc', 'delta_HF', 'delta_GW', 'Znk']
    """

    if not os.path.isfile(os.path.join(file_path, file_name)):
        print('File not found:', os.path.join(file_path, file_name))
        print("Skipping file")
        return {}

    # Value in input can exceed the total number of empty states.
    # The value used by exciting in GW is the smallest 'n_empty' value in KPOINTS.OUT,
    # as each k-point can differ due to the plane-wave cut-off
    # TODO(Alex) This doesn't always. As in, it may not be the minimum work as in
    # n_empty = parse_kpoints(file_path)['n_empty']
    n_empty = get_nempty_from_evalqp(file_path)

    # TODO Note, if n_empty in input is not max number, this will be used for number of entries
    # in EVALQP.DAT, not the lowest value from KPOINTS file

    #  Note, not kpoints in the KPOINTS file
    #  I assume irreducuble number of k pr q? Not sure.
    nkpts_details = grep("k-point", file_path + '/' + file_name).splitlines()
    nkpts = int(nkpts_details[-1].split()[2].replace(':', ''))

    fid = open(file_path + "/" + file_name, "r")
    file_string = fid.readlines()
    fid.close()

    keys = ['E_KS', 'E_HF', 'E_GW', 'sigma_x', 'Re_sigma_c', 'Im_sigma_c', 'V_xc', 'delta_HF', 'delta_GW', 'Znk']

    data = {}
    i = 0
    for ik in range(0, nkpts):
        k_point = [float(k) for k in file_string[i].split()[3:6]]

        # iterate past k-point and skip header
        i += 2
        results = {}
        for istate in range(0, n_empty):
            line = file_string[i].split()[1:]
            # 1-Indexing consistent with fortran
            results[istate + 1] = {keys[i]: float(line[i]) for i in range(0, len(keys))}
            i += 1

        # skips extra blank line per k-point block
        i += 1
        data[ik] = {'k_point': k_point, 'results': results.copy()}

    return data


def parse_gw_timings(file_path: str, file_name='GW_INFO.OUT'):
    """
    Get timings of each part of a GW calculation, from GW_INFO.OUT

    :param file_path: file path
    :param file_name: file name
    :return: dictionary of timings
    """
    file_path += '/' + file_name

    # Get line number GW timing info
    start_line = int(grep("GW timing info", file_path, line_number='').split(':')[0])
    #end_line = int(grep("Total", file_path, line_number='').splitlines()[-1].split(':')[0])
    fid = open(file_path, "r")
    timing_lines = fid.readlines()[start_line+2:]
    fid.close()

    timings = {}
    for line in timing_lines:
        data = line.split()

        # Skip blank lines (data[0] will throw an error)
        if not data:
            continue

        # Remove '-' prefixes
        if data[0] == '-':
            key = " ".join(data[1:-2])
        else:
            key = " ".join(data[0:-2])

        # Don't store blanks
        if len(key.strip()) != 0:
            timings[key] = float(data[-1])

    return timings


def parse_gw_info(file_path: str, file_name='GW_INFO.OUT') -> dict:
    """
    Parse variables from GW_INFO.OUT:
       max_n_lapw,
       min_n_lapw,
       n_KS,
       n_occupied,
       n_unoccupied,
       i_VBM           Index of VB top (for GW)
       i_CBm           Index of CB bottom (for GW)

    Timings returned from a separate routine.

    :param file_path: file path
    :param file_name: file name
    :return: dictionary of important variables
    """
    file_path = os.path.join(file_path, file_name)

    if not os.path.isfile(file_path):
        print('File not found:', file_path)
        return {}

    data = {'max_n_lapw': int(grep("Maximum number of LAPW states", file_path).split()[-1]),
            'min_n_lapw': int(grep("Minimal number of LAPW states", file_path).split()[-1]),
            'n_KS': int(grep("total KS", file_path).split()[-1]),
            'n_occupied': int(grep("occupied", file_path).split()[2]),
            'n_unoccupied': int(grep("occupied", file_path).split()[5]),
            'i_VBM': int(grep("Band index of VBM", file_path).split()[-1]),
            'i_CBm': int(grep("Band index of CBm", file_path).split()[-1])}

    # Save the second line in each case, which corresponds to the GW band indices
    assert data['i_CBm'] == data['i_VBM'] + 1

    return data
