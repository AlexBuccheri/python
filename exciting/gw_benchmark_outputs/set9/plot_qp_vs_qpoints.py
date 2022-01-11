"""
Extract/Plot each quasi-particle gap as a function of q-points
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
    :return: list gaps
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


def process_gw_calculations(root_path: str):
    directories = {'2,2,2': 'refined',
                   '4,4,4': '444_refined',
                   '6,6,6': '666_refined'
                   }

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

    for directory in directories.values():
        path = os.path.join(root_path, directory)
        gw_result = process_gw_calculation(path, a1_gaps())

        # Want it such that I store results['E_qp_Gamma_Gamma'] = [gap @ 2,2,2, gap @ 4,4,4, gap @ 6,6,6]
        for gap_label, gap in gw_result.items():
            results[gap_label].append(gap)

    return results


def plot_qp_vs_qpoints(q_points, results: list, label: str, file_name):
    regular_x_values = [i for i in range(1, len(q_points) + 1)]
    converted_results = [x * ha_to_mev for x in results]

    fig, ax = plt.subplots()
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})

    ax.set_xticks(regular_x_values)
    ax.set_xticklabels(q_points)

    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    ax.set_xlabel('q-point Grid', fontsize=16)
    ax.set_ylabel(label, fontsize=16)
    ax.plot(regular_x_values, converted_results, 'ro-', markersize=10)
    if save_plots:
        plt.savefig(file_name + '.jpeg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
    plt.show()


def print_results(q_points: List[str], results: list):
    converted_results = [x * ha_to_mev for x in results]

    print("Energies in meV")

    header = "# Q-point Grid, QP (Gamma), QP (X-Gamma), QP (X-X), KS (Gamma), KS (X-Gamma), KS (X-X)"
    print(header)
    for i in range(0, len(q_points)):
        print(q_points[i], converted_results['E_qp_Gamma_Gamma'][i],
              converted_results['E_qp_X_Gamma'][i],
              converted_results['E_qp_X_X'][i],
              converted_results['E_ks_Gamma_Gamma'][i],
              converted_results['E_ks_X_Gamma'][i],
              converted_results['E_ks_X_X'][i]
              )

    header = "# Q-point Grid, QP-KS (Gamma), QP-KS (X-Gamma), QP-KS (X-X)"
    print(header)
    for i in range(0, len(q_points)):
        print(q_points[i], converted_results['delta_E_qp_Gamma_Gamma'][i],
                           converted_results['delta_E_qp_X_Gamma'][i],
                           converted_results['delta_E_qp_X_X'][i],
              )


def main(root_path: str):
    """
    TODO Alex
    Set plot (and save) for each gap
    Plot the max number of irr q-points
    Push to github

    :param root_path:
    :return:
    """
    q_point_labels = ['2x2x2', '4x4x4', '6x6x6']
    # From inspection of EVALQP.DAT
    irreducible_q_points = [3, 8, 16]

    results = process_gw_calculations(root_path)

    plot_qp_vs_qpoints(q_point_labels,
                       results['delta_E_qp_Gamma_Gamma'],
                       'Quasiparticle Gap - KS Gap at Gamma (meV)')

    plot_qp_vs_qpoints(q_point_labels,
                       results['delta_E_qp_X_Gamma'],
                       'Quasiparticle Gap - KS Gap (X-Gamma) (meV)')

    plot_qp_vs_qpoints(q_point_labels,
                       results['delta_E_qp_X_X'],
                       'Quasiparticle Gap - KS Gap (X-X) (meV)')

    print_results(q_point_labels, results)


save_plots = True
if __name__ == "__main__":
    root_path = "/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8"
    main(root_path)
