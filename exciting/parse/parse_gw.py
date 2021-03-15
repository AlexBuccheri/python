"""
Parse various GW output files
"""
import subprocess
import warnings

from exciting_utils.py_grep import grep


def parse_gw_evalqp(file_path: str, file_name='EVALQP.DAT', nempty=None, nkpts=None):
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
    nempty: Number of empty states
    nkpts: Number of irreducible q-points, I assume

    :return: dictionary of form {ik: k-point, results},
    where results[istate].keys = ['E_KS', 'E_HF', 'E_GW', 'sigma_x',' Re_sigma_c', 'Im_sigma_c', 'V_xc', 'delta_HF', 'delta_GW', 'Znk']
    """

    if nempty is None:
        string = grep('nempty', file_path + '/input.xml')
        nempty = float(string.split("\"")[1])

    if nkpts is None:
        #TODO Add me
        print("I assume this is the number of irredicuble q-points but need to confirm, hence must pass.")
        quit()

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
        for istate in range(0, nempty):
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