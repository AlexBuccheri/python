""" Post-process silicon G0W0 results.

* Plot indirect gap w.r.t q-points
* Plot direct gap w.r.t. q-points
* Repeat this for each emtpy state setting
"""
import os
from dataclasses import dataclass
from typing import List
import re
import matplotlib.pyplot as plt
import numpy as np

from excitingtools.dataclasses.data_structs import PointIndex, BandIndices
from excitingtools.exciting_obj_parsers.gw_eigenvalues import NitrogenEvalQPColumns, gw_eigenvalue_parser
from excitingtools.dataclasses.eigenvalues import EigenValues


@dataclass(frozen=True)
class CubicBZ:
    Gamma = [0., 0., 0.]
    # TODO(Alex) Check this is actually X
    X = [0., 0.5, 0.5]


def get_gw_bandedge_indices(path: str) -> BandIndices:
    """ Get GW band indices from GW_INFO.OUT, for exciting Nitrogen.

    Note, the strings may differ in Oxygen.

    :param path: Path to file `GW_INFO.OUT`
    :return (vbm, cbm): VBM and CBm indices.
    """
    file_name = os.path.join(path, "GW_INFO.OUT")
    try:
        with open(file=file_name) as fid:
            file_string = fid.read()
    except FileNotFoundError:
        raise FileNotFoundError(f'{file_name} cannot be found')

    # In both cases, take the last match, which corresponds to the GW indice.
    vb_max_str = re.findall(r'\s*Band index of VBM: .*$', file_string, flags=re.MULTILINE)[-1]
    cb_min_str = re.findall(r'\s*Band index of CBM: .*$', file_string, flags=re.MULTILINE)[-1]

    vbm = int(vb_max_str.split()[-1])
    cbm = int(cb_min_str.split()[-1])

    return BandIndices(VBM=vbm, CBm=cbm)


def get_gw_bandedge_k_indices(path: str) -> List[PointIndex]:
    """ Parse k-points at band edges.

    Valid for Nitrogen.

    Fundamental BandGap (eV):                 1.4935
         at k(VBM) =    0.000   0.000   0.000 ik =     1
            k(CBM) =    0.000   0.500   0.500 ik =     3

    :param path:
    :return:
    """
    file_name = os.path.join(path, "GW_INFO.OUT")
    try:
        with open(file=file_name) as fid:
            file_string = fid.read()
    except FileNotFoundError:
        raise FileNotFoundError(f'{file_name} cannot be found')

    # Last line match will correspond to GW
    vbm_line = re.findall(r'^.*k\(VBM\) = .*$', file_string, flags=re.MULTILINE)[-1]
    cbm_line = re.findall(r'^.*k\(CBM\) = .*$', file_string, flags=re.MULTILINE)[-1]

    ik_vbm = int(vbm_line.split()[-1])
    k_vbm = [float(x) for x in vbm_line.split()[3:6]]

    ik_cbm = int(cbm_line.split()[-1])
    k_cbm = [float(x) for x in cbm_line.split()[2:5]]

    return [PointIndex(k_vbm, ik_vbm), PointIndex(k_cbm, ik_cbm)]


def return_zero_if_no_file(func):
    """ Decorate `band_gaps` and return zeros
    if the file cannot be found (as the parser will throw an immediate exception).
    """
    def modified_func(path):
        full_file = os.path.join(path, 'EVALQP.DAT')
        if not os.path.isfile(full_file):
            print(f'{full_file} does not exist')
            return {'fundamental': 0.0, 'gamma': 0.0}
        return func(path)
    return modified_func


@return_zero_if_no_file
def band_gaps(path) -> dict:
    """ Parse eigenvalues, return band gaps in Ha.

    :return: Dict of band gaps
    """
    eigenvalues: EigenValues = gw_eigenvalue_parser(path, NitrogenEvalQPColumns.E_GW)

    # Indirect (or smallest)
    band_indices = get_gw_bandedge_indices(path)
    k_valence, k_conduction = get_gw_bandedge_k_indices(path)
    fundamental_gap = eigenvalues.band_gap(band_indices, k_indices=[k_valence.index, k_conduction.index])

    # Direct at Gamma
    ik = eigenvalues.get_index(CubicBZ.Gamma)
    direct_gamma_gap = eigenvalues.band_gap(band_indices, k_indices=[ik, ik])

    return {'fundamental': fundamental_gap, 'gamma': direct_gamma_gap}


def band_gaps_vs_q(root: str, q_points: list) -> tuple:
    """ Get fundamental and direct Gamma gaps as a function of q.

    :param q_points:
    :return:
    """
    direct_gaps = []
    indirect_gaps = []
    for q_point in q_points:
        q_str = "".join(str(x) for x in q_point)
        path = os.path.join(root, q_str)
        gaps: dict = band_gaps(path)
        indirect_gaps.append(gaps['fundamental'])
        direct_gaps.append(gaps['gamma'])
    return indirect_gaps, direct_gaps


if __name__ == "__main__":

    # Notebook Results directory
    root = 'GW_results/results/silicon'
    subdirectory = '_percent_empty'

    q_points = [[2, 2, 2], [4, 4, 4], [6, 6, 6], [8, 8, 8]]
    total_q_points = [np.prod(q_point) for q_point in q_points]
    n_empties = [25, 50, 100]

    for n_empty in n_empties:
        indirect_gaps, direct_gaps = band_gaps_vs_q(os.path.join(root, str(n_empty) + subdirectory), q_points)
        print(indirect_gaps)
        print(direct_gaps)

        plt.plot(total_q_points, indirect_gaps)
        plt.show()





