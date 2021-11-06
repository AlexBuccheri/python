"""
Run what I consider to be converged from set 8 with:
Same MT setting: (Zr, O) = (1.6, 1.6), but rgkmax =9 (rgkmax improvement o set 8)
Different MT setting: (Zr, O) = (2.0, 1.6), and rgkmax =9 (rgkmax improvement on set 9)
"""
import numpy as np

from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gamma_point, process_gw_gap
from units_and_constants.unit_conversions import ha_to_mev


def process_gw_calculation(path: str) -> dict:
    """
    :param str path: Path to GW files
    :return: dict Processed GW data
    """
    print('Reading data from ', path)
    gw_data = parse_gw_info(path)
    qp_data = parse_gw_evalqp(path)
    results = process_gw_gamma_point(gw_data, qp_data)

    if not results:
        return {'delta_E_qp': [],
                're_self_energy_VBM': [],
                're_self_energy_CBm': []
                }

    X_point = [0., 0.5, 0.5]
    Gamma_point = [0., 0., 0.]

    # Specific to the A1 system X (valence) -> Gamma (conduction)
    results_x_gamma = process_gw_gap(gw_data, qp_data, X_point, Gamma_point)
    results_x_x = process_gw_gap(gw_data, qp_data, X_point, X_point)

    return {'E_qp': np.array(results['E_qp']),
            'E_ks': np.array(results['E_ks']),
            'E_qp_X_Gamma': np.array(results_x_gamma['E_qp']),
            'E_ks_X_Gamma': np.array(results_x_gamma['E_ks']),
            'E_qp_X_X': np.array(results_x_x['E_qp']),
            'E_ks_X_X': np.array(results_x_x['E_ks']),

            'delta_E_qp': np.array(results['E_qp'] - results['E_ks']),
            're_self_energy_VBM': np.array(results['re_sigma_VBM']),
            're_self_energy_CBm': np.array(results['re_sigma_CBm'])
            }


def print_results(data: dict, message: str):

    print('Printing data for ' + message)
    print('LO cut-off (Ha), QP(G-G), QP(X-G), QP(X-X), KS(G-G)')
    try:
        qp_g_g = data['E_qp']* ha_to_mev
        qp_x_g = data['E_qp_X_Gamma'] * ha_to_mev
        qp_x_x = data['E_qp_X_X'] * ha_to_mev
        ks_g_g = data['E_ks'] * ha_to_mev
        ks_x_g = data['E_ks_X_Gamma'] * ha_to_mev
        ks_x_x = data['E_ks_X_X'] * ha_to_mev
        print(qp_g_g, qp_x_g, qp_x_x, ks_g_g, ks_x_g, ks_x_x)
    except KeyError:
        print('Index, energy not computed')


def rgkmax_and_mt_effect():
    """
    set 8 uses MT (Zr, O) = (1.6, 1.6)
    set 9 uses MT (Zr, O) = (2.0, 1.6)
    """

    l_max_pairs = [{'zr': 6, 'o': 5}]
    energy_index = 'i1'
    energy_cutoff = 100
    muffin_tins = {'equal': {'zr': 1.6, 'o': 1.6},
                   'correct_ratio': {'zr': 2.0, 'o': 1.6}
                   }

    # PWs in basis = 483, LOs = 807
    rgk8_mt_equal = process_gw_calculation('/users/sol/abuccheri/gw_benchmarks/A1_set8/zr_lmax6_o_lmax5_rgkmax8/gw_q222_omeg32_nempty2000/max_energy_i1')
    # PWs in basis = 680, LOs = 807
    rgk9_mt_equal = process_gw_calculation('/users/sol/abuccheri/gw_benchmarks/A1_set8/zr_lmax6_o_lmax5_rgkmax9/gw_q222_omeg32_nempty2000/max_energy_i1')

    # PWs in basis = 483, LOs in basis = 878
    rgk8_mt_ratio = process_gw_calculation('/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8/gw_q222_omeg32_nempty3000/max_energy_i1')
    # PWs in basis = 680, LOs in basis = 875
    # TODO Note the LO I had to remove from the basis in the GW calc
    rgk9_mt_ratio = process_gw_calculation('/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax9/gw_q222_omeg32_nempty3000/max_energy_i1')

    print_results(rgk8_mt_equal, "(6,5), 100 Ha cutoff, rgkmax8, MT (1.6, 1.6)")
    print_results(rgk9_mt_equal, "(6,5), 100 Ha cutoff, rgkmax9, MT (1.6, 1.6)")
    print_results(rgk8_mt_ratio, "(6,5), 100 Ha cutoff, rgkmax8, MT (2.0, 1.6)")
    print_results(rgk9_mt_ratio, "(6,5), 100 Ha cutoff, rgkmax9, MT (2.0, 1.6)")


if __name__ == "__main__":
    rgkmax_and_mt_effect()
