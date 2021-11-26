"""
Plot set 9 data
# For each l-max pair, plot QP vs
"""
import numpy as np
import os
from typing import List
import matplotlib.pyplot as plt
from collections import OrderedDict

from parse.set_gw_input import GWInput
from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gamma_point, process_gw_gap
from units_and_constants.unit_conversions import ha_to_mev
from gw_benchmark_outputs.post_process_utils import get_basis_labels

# Energy cutoffs
from gw_benchmark_inputs.set9.basis import set_lo_channel_cutoffs, n_energies_per_channel

# GLOBAL
save_plots = True


# These path functions should be shared by input and outputs of a given set
def zro2_path(l_max: dict) -> str:
    """

    :param l_max:
    :return:
    """
    return 'zr_lmax' + str(l_max['zr']) + '_o_lmax' + str(l_max['o']) + '_rgkmax8'


def get_gw_path(gw_settings: GWInput) -> str:
    """

    :param gw_settings:
    :return:
    """
    gw_s = gw_settings
    q_str = "".join(str(q) for q in gw_s.ngridq)
    return "gw_q" + q_str + "_omeg" + str(gw_s.nomeg) + "_nempty" + str(gw_s.nempty)


def get_l_max_key(species: List[str], l_max: dict):
    """

    :param species:
    :param l_max:
    :return:
    """
    key = '('
    for x in species:
        key += str(l_max[x]) + ','
    return key[:-1] + ')'


def process_gw_calculation(path: str) -> dict:
    """
    This could be a more generic routine
    Add other gaps to it as well
    # TODO(Alex) Check parse_gw_evalqp parser is valid (and replace with one from exciting tools)

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


def process_gw_calculations(root: str,
                            species: List[str],
                            l_max_pairs: dict,
                            gw_settings: GWInput) -> dict:
    """
    Return GW calculation data, for the directory structure used by set8
    Loop over a) l_max pairs and b) energy_cutoffs

    :return: dict data: GW data
    """

    gw_path = get_gw_path(gw_settings)
    data = {}

    for l_max in l_max_pairs:
        l_max_key = get_l_max_key(species, l_max)
        energy_cutoffs = set_lo_channel_cutoffs(l_max)
        data[l_max_key] = {}

        for i in range(0, len(energy_cutoffs)):
            full_path = os.path.join(root, zro2_path(l_max), gw_path, 'max_energy_i' + str(i))
            data[l_max_key][i] = process_gw_calculation(full_path)

    return data


def print_results(data: dict, energy_cutoff, lmax_str: str):

    print('Printing data for ' + lmax_str)
    print('LO cut-off (Ha), QP(G-G), QP(X-G), QP(X-X), KS(G-G)')
    for ie, energy in enumerate(energy_cutoff):
        try:
            qp_g_g = data[lmax_str][ie]['E_qp']
            qp_g_g = qp_g_g * ha_to_mev
            qp_x_g = data[lmax_str][ie]['E_qp_X_Gamma'] * ha_to_mev
            qp_x_x = data[lmax_str][ie]['E_qp_X_X'] * ha_to_mev
            ks_g_g = data[lmax_str][ie]['E_ks'] * ha_to_mev
            ks_x_g = data[lmax_str][ie]['E_ks_X_Gamma'] * ha_to_mev
            ks_x_x = data[lmax_str][ie]['E_ks_X_X'] * ha_to_mev
            print(energy, qp_g_g, qp_x_g, qp_x_x, ks_g_g, ks_x_g, ks_x_x)
        except KeyError:
            print('Index, energy:', ie, energy, 'not computed')

    return


def plot_data(l_max_pairs, data, basis_labels):
    fig, ax = plt.subplots()
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})

    ax.set_xlabel('l_max (Zr, O)', fontsize=16)
    ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)', fontsize=16)

    # Get x data in useful form
    x = np.arange(0, len(l_max_pairs))
    x_label_to_i = {'(6,5)': 0, '(7,6)': 1}

    # Assume the l_max_pairs dict is ordered
    x_keys = ["(" + ",".join(str(l) for l in l_pair.values()) + ")" for l_pair in l_max_pairs]

    ax.set_xticks(x)
    ax.set_xticklabels(x_keys)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Plot each energy cutoff separately so missing data can easily be skipped
    energy_cutoff = [80, 100, 120, 150, 180, 200]
    for ie in range(0, n_energies_per_channel):
        x_ie = []
        y_ie = []
        for l_key in x_keys:
            qp_ks = data[l_key][ie]['delta_E_qp']
            if qp_ks:
                # print(ie, l_key, qp_ks)
                x_ie.append(x_label_to_i[l_key])
                y_ie.append(qp_ks * ha_to_mev)
        ax.plot(x_ie, y_ie, marker='o', markersize=8, label=str(energy_cutoff[ie]))

    ax.legend(loc='lower right', title='LO cutoff (Ha)')

    def make_label(basis_labels) -> str:
        label = ''
        for species in ['zr', 'o']:
            label += species.capitalize() + ': ' + basis_labels[l_key][species][ie].rstrip() + '.\n'
        return label

    # Add basis label to each point
    # Commented because the plot gets WAY too busy with the number of labels
    # for ie in range(0, n_energies_per_channel):
    #     i = 0
    #     for l_key in x_keys:
    #         qp_ks = data[l_key][ie]['delta_E_qp']
    #         if qp_ks:
    #             # Set basis label for each marker
    #             label = make_label(basis_labels)
    #             ax.annotate(label.rstrip(), (x[i], qp_ks * ha_to_mev))
    #         i += 1

    if save_plots:
        plt.savefig('qp_lmax_LO_sweep.jpeg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
    plt.show()

    return


def plot_data_for_set8_and_set9(set8, set9):
    fig, ax = plt.subplots()
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})

    ax.set_xlabel('l_max (Zr, O)', fontsize=16)
    ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)', fontsize=16)
    ax.legend(loc='lower right', title='LO cutoff (Ha)')

    # Get x data in useful form
    n_l_pairs = 4
    x = np.arange(0, n_l_pairs)
    x_label_to_i = {'(4,3)': 0, '(5,4)': 1, '(6,5)': 2, '(7,6)': 3}

    ax.set_xticks(x)
    ax.set_xticklabels(list(x_label_to_i))
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Plot each energy cutoff separately so missing data can easily be skipped
    def plot_data(ax, data, x_keys, energy_cutoff):
        for ie in range(0, n_energies_per_channel):
            x_ie = []
            y_ie = []
            for l_key in x_keys:
                qp_ks = data[l_key][ie]['delta_E_qp']
                if qp_ks:
                    # print(ie, l_key, qp_ks)
                    x_ie.append(x_label_to_i[l_key])
                    y_ie.append(qp_ks * ha_to_mev)
            ax.plot(x_ie, y_ie, marker='o', markersize=8, label=str(energy_cutoff[ie]))
        return ax

    ax = plot_data(ax, set8, ['(4,3)', '(5,4)', '(6,5)', '(7,6)'], [80, 100, 120, 150, 180, 200, 250])
    ax = plot_data(ax, set9, ['(6,5)', '(7,6)'], [80, 100, 120, 150, 180, 200])

    if save_plots:
        plt.savefig('qp_lmax_LO_sweep.jpeg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
    plt.show()



def plot_65_data(data, basis_labels):

    fig, ax = plt.subplots()
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})

    ax.set_xlabel('LO Energy Cutoff (Ha)', fontsize=16)
    ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)', fontsize=16)

    plt.xlim(75, 275)
    plt.ylim(2100, 2108)

    def make_label(basis_labels: OrderedDict) -> str:
        label = ''
        for species in ['zr', 'o']:
            label += species.capitalize() + ': ' + basis_labels[l_key][species][ie].rstrip() + '.\n'
        return label

    # Plot each energy cutoff separately so missing data can easily be skipped
    l_key = '(6,5)'
    energy_cutoffs = [80, 100, 120, 150, 180, 200]

    x_ie = []
    y_ie = []
    for ie, energy_cutoff in enumerate(energy_cutoffs):
        qp_ks = data[l_key][ie]['delta_E_qp']
        if qp_ks:
            # Set basis label for each marker
            label = make_label(basis_labels)
            ax.annotate(label.rstrip(), (energy_cutoff, qp_ks * ha_to_mev), fontsize=12)
            # Store points
            x_ie.append(energy_cutoff)
            y_ie.append(qp_ks * ha_to_mev)

    ax.plot(x_ie, y_ie, marker='o', markersize=8)
    if save_plots:
        plt.savefig('qp_lmax65.jpeg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
    plt.show()

    return


def get_set8_data(root_set8):
    """
    rgkmax = 8
    MT (Zr, O) = (1.6, 1.6)
    Energies in meV

    :param root_set8:
    :return:
    """

    data = {}

    energies_43 = {80:  [],
                   100: [],
                   120: [],
                   150: [],
                   180: [],
                   200: [],
                   250: []}

    energies_54 = {80:  [],
                   100: [],
                   120: [],
                   150: [],
                   180: [],
                   200: [],
                   250: []}

    #                  QP(G-G) | QP(X-G) | QP(X-X) | KS(G-G) | KS(X-G) | KS(X-X) |
    energies_65 = {80:  [5972,  5403, 5560, 3865, 3321, 3748],
                   100: [5968,  5399, 5557, 3865, 3321, 3748],
                   120: [5967,  5398, 5555, 3865, 3321, 3748],
                   150: [5968,  5399, 5557, 3865, 3321, 3748],
                   180: [5967,  5398, 5556, 3865, 3321, 3748],
                   200: [5972,  5403, 5561, 3865, 3321, 3748],
                   250: [5967,  5398, 5557, 3865, 3321, 3748]}

    for ie in [80, 100, 120, 150, 180, 200, 250]:
        data['(6,5)'][ie]['delta_E_qp'] = energies_65[ie][0] - energies_65[ie][2]

    energies_76 = {80:  [5975, 5407, 5564, 3865, 3321, 3748],
                   100: [5972, 5403, 5561, 3865, 3321, 3748],
                   120: [5974, 5404, 5562, 3865, 3321, 3748]}

    for ie in [80, 100, 120]:
        data['(7,6)'][ie]['delta_E_qp'] = energies_76[ie][0] - energies_76[ie][2]

    return data


def basis_convergence(root):
    """
    Main routine for plotting basis convergence
    :param root:
    :return:
    """
    # Settings
    species = ['zr', 'o']
    l_max_pairs = [{'zr': 6, 'o': 5}, {'zr': 7, 'o': 6}]
    gw_settings = GWInput(taskname="g0w0",
                          nempty=3000,
                          ngridq=[2, 2, 2],
                          skipgnd=False,
                          n_omega=32,
                          freqmax=1.0)

    data = process_gw_calculations(root, species, l_max_pairs, gw_settings)

    # TODO API for generating basis labels needs to change
    l_max_values = [OrderedDict([('zr', 6), ('o', 5)]),
                    OrderedDict([('zr', 7), ('o', 6)])]

    energy_indices = ['i0', 'i1', 'i2', 'i3', 'i4', 'i5']
    energy_cutoffs = [80, 100, 120, 150, 180, 200]

    settings = {'rgkmax': 8,
                'l_max_values': l_max_values,
                'n_img_freq': 32,
                'q_grid': [2, 2, 2],
                'n_empty_ext': [gw_settings.nempty] * len(l_max_values),
                # Extensions, not energies
                'max_energy_cutoffs': energy_indices}

    # Post-process
    basis_labels = get_basis_labels(root, settings)

    print_results(data, energy_cutoffs, '(6,5)')
    print_results(data, energy_cutoffs, '(7,6)')

    # Plot just for set 9
    # plot_data(l_max_pairs, data, basis_labels)
    # TODO(Alex) Fix plotting for this
    # plot_65_data(data, basis_labels)

    # Plot for set 8 and set 9
    # TODO(Alex) Test they both load on the same plot
    # If so, fix the colour scheme and plot
    set8_data = get_set8_data("/users/sol/abuccheri/gw_benchmarks/A1_set8/")
    plot_data_for_set8_and_set9(set8_data, data)

    return


def main():
    basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_set9/")


if __name__ == "__main__":
    main()
