"""
Parse various GW output files
"""
import subprocess
import warnings
import numpy as np
from pathlib import Path

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



def parse_gw_evalqp(file_path: str, file_name='EVALQP.DAT'):
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
    if not Path(file_path + '/' + file_name).is_file():
        print("File does not exist:", file_path)
        quit("Routine 'parse_gw_evalqp' quit")

    # Value in input can exceed the total number of empty states.
    # The value used by exciting in GW is the smallest 'n_empty' value in KPOINTS.OUT,
    # as each k-point can differ due to the plane-wave cut-off
    n_empty = parse_kpoints(file_path)['n_empty']

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
        i+=2
        results = {}
        for istate in range(0, np.min(n_empty)):
            line = file_string[i].split()[1:]
            # 1-Indexing consistent with fortran
            results[istate + 1] = {keys[i]: float(line[i]) for i in range(0, len(keys))}
            i += 1

        # skips extra blank line per k-point block
        i += 1
        data[ik] = {'k_point':k_point, 'results':results.copy()}

    return data


def parse_gw_timings(file_path: str, file_name='GW_INFO.OUT'):
    """
    Get timings of each part of a GW calculation, from GW_INFO.OUT

    :param file_path: file path
    :param file_name: file name
    :return: dictionary of timings
    """
    file_path += '/' + file_name

    # Get line number  GW timing info
    line_number = int(grep("GW timing info", file_path, line_number='').split(':')[0])
    fid = open(file_path, "r")
    timing_lines = fid.readlines()[line_number+2:]
    fid.close()

    timings = {}
    for line in timing_lines:
        data = line.split()

        # Remove '-' prefixes
        if data[0] == '-':
            key = " ".join(data[1:-2])
        else:
            key = " ".join(data[0:-2])

        # Don't store blanks
        if len(key.strip()) != 0:
            timings[key] = float(data[-1])

    return timings


def parse_gw_info(file_path: str, file_name='GW_INFO.OUT'):
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
    file_path += '/' + file_name

    data = {}
    data['max_n_lapw'] = int(grep("Maximum number of LAPW states", file_path).split()[-1])
    data['min_n_lapw'] = int(grep("Minimal number of LAPW states", file_path).split()[-1])
    data['n_KS'] = int(grep("total KS", file_path).split()[-1])
    data['n_occupied'] = int(grep("occupied", file_path).split()[2])
    data['n_unoccupied'] = int(grep("occupied", file_path).split()[5])

    # Save the second line in each case, which corresponds to the GW band indices
    data['i_VBM'] = int(grep("Band index of VBM", file_path).split()[-1])
    data['i_CBm'] = int(grep("Band index of CBm", file_path).split()[-1])
    assert data['i_CBm'] == data['i_VBM'] + 1

    return data