"""
Plot set 8 data
# For each l-max pair, plot QP vs
"""
import numpy as np
import os
from typing import List
import matplotlib.pyplot as plt
from collections import OrderedDict

from parse.set_gw_input import GWInput
from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gamma_point
from units_and_constants.unit_conversions import ha_to_mev
from gw_benchmark_outputs.post_process_utils import get_basis_labels

# Energy cutoffs
from gw_benchmark_inputs.set8.basis import set_lo_channel_cutoffs, n_energies_per_channel


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

    if results:
        return {'delta_E_qp': np.array(results['E_qp'] - results['E_ks']),
                're_self_energy_VBM': np.array(results['re_sigma_VBM']),
                're_self_energy_CBm': np.array(results['re_sigma_CBm'])
                }
    else:
        return {'delta_E_qp': [],
                're_self_energy_VBM': [],
                're_self_energy_CBm': []
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


def print_results(species: List[str], l_max_pairs: dict, data: dict):
    i = 3  # 0,1,2,3
    l_channel = 0

    for l_max in l_max_pairs:
        l_max_key = get_l_max_key(species, l_max)
        energy_cutoffs = set_lo_channel_cutoffs(l_max)
        # Zr and O cut-offs are the same, and consistent per l_channel
        # energy_cutoffs['zr'][l_channel][i]

        qp = data[l_max_key][i]['delta_E_qp'] * ha_to_mev if data[l_max_key][i]['delta_E_qp'] else None
        print(l_max_key, i, qp)

    return


def plot_data(l_max_pairs, data, basis_labels):
    fig, ax = plt.subplots()
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})

    ax.set_xlabel('l_max (Zr, O)', fontsize=16)
    ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)', fontsize=16)

    # Get x data in useful form
    x = np.arange(0, len(l_max_pairs))
    x_label_to_i = {'(4,3)': 0, '(5,4)': 1, '(6,5)': 2, '(7,6)': 3}

    # Assume the l_max_pairs dict is ordered
    x_keys = ["(" + ",".join(str(l) for l in l_pair.values()) + ")" for l_pair in l_max_pairs]

    ax.set_xticks(x)
    ax.set_xticklabels(x_keys)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Plot each energy cutoff separately so missing data can easily be skipped
    energy_cutoff = [80, 100, 120, 150, 180, 200, 250]
    for ie in range(0, n_energies_per_channel):
        x_ie = []
        y_ie = []
        for l_key in x_keys:
            qp_ks = data[l_key][ie]['delta_E_qp']
            if qp_ks:
                print(ie, l_key, qp_ks)
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

    plt.show()

    return


def main():
    # Settings
    root = "/users/sol/abuccheri/gw_benchmarks/A1_set8/"
    species = ['zr', 'o']
    l_max_pairs = [{'zr': 4, 'o': 3}, {'zr': 5, 'o': 4}, {'zr': 6, 'o': 5}, {'zr': 7, 'o': 6}]
    gw_settings = GWInput(taskname="g0w0",
                          nempty=2000,
                          ngridq=[2, 2, 2],
                          skipgnd=False,
                          n_omega=32,
                          freqmax=1.0)

    data = process_gw_calculations(root, species, l_max_pairs, gw_settings)

    print_results(species, l_max_pairs, data)

    # TODO API for generating basis labels needs to change
    l_max_values = [OrderedDict([('zr', 4), ('o', 3)]),
                    OrderedDict([('zr', 5), ('o', 4)]),
                    OrderedDict([('zr', 6), ('o', 5)]),
                    OrderedDict([('zr', 7), ('o', 6)])]

    settings = {'rgkmax': 8,
                'l_max_values': l_max_values,
                'n_img_freq': 32,
                'q_grid': [2, 2, 2],
                'n_empty_ext': [2000] * len(l_max_values),
                # Extensions, not energies
                'max_energy_cutoffs': ['i0', 'i1', 'i2', 'i3', 'i4', 'i5', 'i6']}

    basis_labels = get_basis_labels(root, settings)

    plot_data(l_max_pairs, data, basis_labels)

    return None


if __name__ == "__main__":
    main()
