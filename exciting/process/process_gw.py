"""


"""
from typing import List
import sys


# For a given basis i.e. (l-max, l-max) energy cutoff == basis functions
def process_gw_gamma_point(gw_data: dict, qp_data: dict) -> dict:
    """
    Process GW data dictionaries and return the direct quasi-particle gap, KS gap,
    VBT and CBB real self-energies

    :param gw_data:
    :param qp_data:
    :return: dict result: Direct gap for GW and KS. re_sigma_VBM, re_sigma_CBm
    """
    if (not gw_data) or (not qp_data):
        return {}

    i_VBM = gw_data['i_VBM']
    i_CBm = gw_data['i_CBm']

    assert qp_data[0]['k_point'] == [0.0, 0.0, 0.0], "Gamma_point always first index in EVALQP.DAT"
    qp_gamma = qp_data[0]['results']

    result = {}
    result['E_ks'] = qp_gamma[i_CBm]['E_KS'] - qp_gamma[i_VBM]['E_KS']
    result['E_qp'] = qp_gamma[i_CBm]['E_GW'] - qp_gamma[i_VBM]['E_GW']

    result['re_sigma_VBM'] = qp_gamma[i_VBM]['Re_sigma_c']
    result['re_sigma_CBm'] = qp_gamma[i_CBm]['Re_sigma_c']

    return result


def process_gw_gap(gw_data: dict, qp_data: dict, kpoint_v: List[float], kpoint_c: List[float]) -> dict:
    """
    Given a specified k-point for the valence band and conduction band, define the
    quasi-particle gap between the two, the KS gap between the two, and the real self-energy at
    in the highest valence band at point kpoint_v, and the real self-energy in the lowest
    conduction band at kpoint_c.

    :param dict gw_data:
    :param dict qp_data: Quasi-particle data of the form
                    qp_data[ik] ={'k_point': k_point, 'results': results}
    :param List[float] kpoint_v: User-defined valence band k-point
    :param List[float] kpoint_c: User-defined conduction band k-point
    :return:
    """
    assert len(kpoint_v) == 3, "Valence band k-point should have length 3"
    assert len(kpoint_c) == 3, "Conduction band k-point should have length 3"

    if (not gw_data) or (not qp_data):
        return {}

    # Index of highest valence band
    i_VBM = gw_data['i_VBM']

    # Index of lowest conduction band
    i_CBm = gw_data['i_CBm']

    # Get energy at specific k-point
    def extract_specific_kpoint(qp_data: dict, kpoint: List[float]):
        for k in qp_data.values():
            if k['k_point'] == kpoint:
                return k['results']
        msg = 'User-specified k-point, ' + ",".join(str(k) for k in kpoint) + \
              'was not found in qp_data parsed from EVALQP.DAT'
        sys.exit(msg)

    data_valance = extract_specific_kpoint(qp_data, kpoint_v)
    data_conduction = extract_specific_kpoint(qp_data, kpoint_c)

    result = {'E_ks': data_conduction[i_CBm]['E_KS'] - data_valance[i_VBM]['E_KS'],
              'E_qp': data_conduction[i_CBm]['E_GW'] - data_valance[i_VBM]['E_GW'],
              're_sigma_VBM': data_valance[i_VBM]['Re_sigma_c'],
              're_sigma_CBm': data_conduction[i_CBm]['Re_sigma_c']}

    return result
