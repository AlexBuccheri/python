"""
Given a variation in LOs, for each channel of the (Zr, O) = (6, 5) basis (i.e. 13 l-channels)
read and plot the effect on the QP energy w.r.t. highest-energy LOs

Tolerate ~ 1 meV change in QP direct gap got a change of basis in any given
l-channel.

Converged GW settings
 MT radius (Zr, O) = (2, 1.6)
 rgkmax = 8
 L_max=(6,5)
 LO cut-off = 100 Ha
"""
import numpy as np
import os
from typing import List, Optional

from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gap
from units_and_constants.unit_conversions import ha_to_mev

# Converged basis directory:
# /users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8/gw_q222_omeg32_nempty3000/max_energy_i1
converged_energies = {
    'E_qp_Gamma_Gamma': 5941,
    'E_ks_Gamma_Gamma': 3863,
    'E_qp_X_Gamma': 5373,
    'E_ks_X_Gamma': 3318,
    'E_qp_X_X': 5542,
    'E_ks_X_X': 3748
}


class Gap:
    def __init__(self, v: list, c: list, v_label: Optional[str]='', c_label: Optional[str]=''):
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

    # Haven't tested these, but assume same irrespective of k-point?:
    #gw_results['re_self_energy_VBM_' + gap.v_label] = np.array(results['re_sigma_VBM'])
    #gw_results['re_self_energy_CBm' + gap.c_label] = np.array(results['re_sigma_CBm'])

    return gw_results


def print_results(species: str, l_channel: int, cutoffs: list, gw_results: dict):
    """
    Print results for varying LOs in one l-channel/

    :param species:
    :param l_channel:
    :param cutoffs:
    :param gw_results:
    :return:
    """
    print('Printing data for ' + species + ' l-channel ' + str(l_channel))
    print('LO cut-off (Ha), QP(G-G), QP(X-G), QP(X-X), KS(G-G)')

    # All data is relative to this converged calculation
    print('100', converged_energies['E_qp_Gamma_Gamma'],
          converged_energies['E_qp_X_Gamma'],
          converged_energies['E_qp_X_X'],
          converged_energies['E_ks_Gamma_Gamma'],
          converged_energies['E_ks_X_Gamma'],
          converged_energies['E_ks_X_X']
          )

    for ie, energy in enumerate(cutoffs):
        try:
            # Convert to meV and around to nearest integer
            rounded_results = {key: int(round(value * ha_to_mev, 0)) for key, value in gw_results[energy].items()}

            qp_g_g = rounded_results['E_qp_Gamma_Gamma']
            qp_x_g = rounded_results['E_qp_X_Gamma']
            qp_x_x = rounded_results['E_qp_X_X']
            ks_g_g = rounded_results['E_ks_Gamma_Gamma']
            ks_x_g = rounded_results['E_ks_X_Gamma']
            ks_x_x = rounded_results['E_ks_X_X']

            print(energy, qp_g_g, qp_x_g, qp_x_x, ks_g_g, ks_x_g, ks_x_x)

        except KeyError:
            print('(Index, energy):', ie, energy, 'not computed')


def main(root_path: str):
    """
    Main

    :param str root_path: Path to calculations
    """
    channel_cutoffs = {'zr': {0: [90, 70, 50], # zr_channel0: max_energy_0, max_energy_1, max_energy_2
                              1: [90, 70, 50], # zr_channel1: max_energy_0, max_energy_1, max_energy_2
                              2: [90, 70, 50],
                              3: [80, 60, 45],
                              4: [75, 55, 40],
                              5: [85, 65, 45],
                              6: [75, 55, 40]},

                       'o': {0: [70, 45, 25],
                             1: [85, 60, 40],
                             2: [75, 50, 30],
                             3: [90, 60, 40],
                             4: [75, 50, 30],
                             5: [85, 60, 40]}
                       }

    gaps = a1_gaps()

    for l_channel, cutoffs in channel_cutoffs['zr'].items():
        channel_dir = 'zr_channel' + str(l_channel)
        results = {}
        for i in range(0, len(cutoffs)):
            full_path = os.path.join(root_path, channel_dir, "max_energy_" + str(i))
            results[cutoffs[i]] = process_gw_calculation(full_path, gaps)
        print_results('zr', l_channel, cutoffs, results)

    for l_channel, cutoffs in channel_cutoffs['o'].items():
        channel_dir = 'o_channel' + str(l_channel)
        results = {}
        for i in range(0, len(cutoffs)):
            full_path = os.path.join(root_path, channel_dir, "max_energy_" + str(i))
            results[cutoffs[i]] = process_gw_calculation(full_path, gaps)
        print_results('o', l_channel, cutoffs, results)


if __name__ == "__main__":
    root_path = "/users/sol/abuccheri/gw_benchmarks/A1_set9/zr_lmax6_o_lmax5_rgkmax8/channel_refinement/"
    main(root_path)
