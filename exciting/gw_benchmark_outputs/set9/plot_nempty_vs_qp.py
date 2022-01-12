"""
QP gap as a function of the number of empty states.
"""
import matplotlib.pyplot as plt
import numpy as np
import os
from typing import List, Optional

from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gap
from units_and_constants.unit_conversions import ha_to_mev


class Gap:
    def __init__(self, v: list, c: list, v_label: Optional[str] = '', c_label: Optional[str] = ''):
        self.v = v
        self.c = c
        self.v_label = v_label
        self.c_label = c_label


def a1_gaps():
    """
    A1 CB-VB gaps of interest
    :return:
    """
    X = [0., 0.5, 0.5]
    Gamma = [0., 0., 0.]

    gaps = [Gap(Gamma, Gamma, 'Gamma', 'Gamma'),
            Gap(X, Gamma, 'X', 'Gamma'),
            Gap(X, X, 'X', 'X')
            ]

    return gaps


def process_gw_calculation(path: str, gaps: List[Gap]) -> dict:
    """
    Process GW calculation results for VB and CB points specified in gaps

    :param str path: Path to GW files
    :param List[Gap] gaps: List of VB and CB points

    :return: dict gw_results: Processed GW data for QP and KS energies for specified gaps
    """
    print('Reading data from ', path)
    gw_data = parse_gw_info(path)
    qp_data = parse_gw_evalqp(path)

    gw_results = {}

    for gap in gaps:
        results = process_gw_gap(gw_data, qp_data, gap.v, gap.c)
        gap_label = gap.v_label + '_' + gap.c_label
        try:
            gw_results['E_qp_' + gap_label] = np.array(results['E_qp'])
            gw_results['E_ks_' + gap_label] = np.array(results['E_ks'])
            gw_results['delta_E_qp_' + gap_label] = np.array(results['E_qp'] - results['E_ks'])
        except KeyError:
            return gw_results

    return gw_results


def process_gw_calculations(root_path: str, nempty_range: List[int]) -> dict:
    gaps = a1_gaps()
    results = {'E_qp_Gamma_Gamma': [],
               'E_ks_Gamma_Gamma': [],
               'delta_E_qp_Gamma_Gamma': [],

               'E_qp_X_Gamma': [],
               'E_ks_X_Gamma': [],
               'delta_E_qp_X_Gamma': [],

               'E_qp_X_X': [],
               'E_ks_X_X': [],
               'delta_E_qp_X_X': [],
               }

    for nempty in nempty_range:
        path = os.path.join(root_path, str(nempty))
        gw_result = process_gw_calculation(path, gaps)

        # Want it such that I store results['E_qp_Gamma_Gamma'] = Gap w.r.t nempty, etc
        for key, gap in gw_result.items():
            results[key].append(gap)

    return results


def plot_qp_vs_nempty(nempty_range, results: list, label: str, file_name='output.jpeg'):
    converted_results = [x * ha_to_mev for x in results]

    fig, ax = plt.subplots()
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})

    ax.set_xlabel('Number of Empty States', fontsize=16)
    ax.set_ylabel(label, fontsize=16)
    ax.plot(nempty_range, converted_results, 'ro-', markersize=10)
    if save_plots:
        plt.savefig(file_name, dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
    plt.show()


def print_results(nempty_range: List[int], results: list, label: str):
    print(label)
    converted_results = [x * ha_to_mev for x in results]
    for i in range(0, len(nempty_range)):
        print(nempty_range[i], converted_results[i])


def main(root_path: str, nempty_range: List[int]):
    results = process_gw_calculations(root_path, nempty_range)

    # Plots
    plot_qp_vs_nempty(nempty_range,
                      results['delta_E_qp_Gamma_Gamma'],
                      'Quasiparticle Gap - KS Gap at Gamma (meV)',
                      file_name="qp_GG_vs_nempty.jpeg")

    plot_qp_vs_nempty(nempty_range,
                      results['delta_E_qp_X_Gamma'],
                      'Quasiparticle Gap - KS Gap (X-Gamma) (meV)',
                      file_name="qp_XG_vs_nempty.jpeg")

    plot_qp_vs_nempty(nempty_range,
                      results['delta_E_qp_X_X'],
                      'Quasiparticle Gap - KS Gap (X-X) (meV)',
                      file_name="qp_XX_vs_nempty.jpeg")

    # Printing
    print_results(nempty_range,
                  results['delta_E_qp_Gamma_Gamma'],
                  'Quasiparticle Gap - KS Gap at Gamma (meV)')

    print_results(nempty_range,
                  results['delta_E_qp_X_Gamma'],
                  'Quasiparticle Gap - KS Gap (X-Gamma( (meV)')

    print_results(nempty_range,
                  results['delta_E_qp_X_X'],
                  'Quasiparticle Gap - KS Gap (X-X) (meV)')


save_plots = True
if __name__ == "__main__":
    # 1187 = max available. I.e. 1200 all - 21 occupied
    root_path = "/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8/nempty_vs_qp"
    nempty_range = [200, 400, 600, 800, 1000, 1100, 1187]
    main(root_path, nempty_range)
