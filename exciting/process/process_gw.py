"""


"""

#units = ['ev', 'mev', 'ha']
#assert unit in units, "input unit not valid"


# For a given basis i.e. (l-max, l-max) energy cutoff == basis functions
def process_gw_gamma_point(gw_data, qp_data):
    """
    Process GW data dictionaries and return the quasiparticle gap - KS gap,
    VBT and CBB real self-energies

    :param gw_data:
    :param qp_data:
    :return:
    """

    i_VBM = gw_data['i_VBM']
    i_CBm = gw_data['i_CBm']

    # Gamma_point always first index in EVALQP.DAT
    assert qp_data[0]['k_point'] == [0.0, 0.0, 0.0]
    qp_gamma = qp_data[0]['results']

    result = {}
    result['E_ks'] = qp_gamma[i_CBm]['E_KS'] - qp_gamma[i_VBM]['E_KS']
    result['E_qp'] = qp_gamma[i_CBm]['E_GW'] - qp_gamma[i_VBM]['E_GW']

    result['re_sigma_VBM'] = qp_gamma[i_VBM]['Re_sigma_c']
    result['re_sigma_CBm'] = qp_gamma[i_CBm]['Re_sigma_c']

    return result