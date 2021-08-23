"""
Plots
"""
# External libraries
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

# My modules
from units_and_constants.unit_conversions import ha_to_mev
from gw_benchmark_outputs.post_process_utils import parse_gw_results, get_basis_labels, combine_species_basis_labels, \
    n_local_orbitals


def process_basis_numbers(delta_E_qp, max_energy_exts:list):
    """
    Print change in QP gap for each max energy parameter (i.e. max LO)
    """
    assert delta_E_qp.size == len(max_energy_exts), 'Should be same number of QP energies as there are basis energy cutoffs'

    print("For basis (Zr,O) = (3,2)")
    for ie, energy in enumerate(max_energy_exts[1:], start=1):
        change_in_delta_E_qp = (delta_E_qp[ie] - delta_E_qp[ie-1])
        print("Change in Delta E_QP from max energy param: ")
        print(max_energy_exts[ie-1].replace("\n", " "), "to ")
        print(energy.replace("\n", " "), ":")
        print(change_in_delta_E_qp, "(meV)")
    return


def sum_los_per_species(n_los_species_resolved: dict) ->list:
    """
    Sum N LOs for each species, per calculation.

    :return: List, n_los, containing the number of LOs per calculation.
    """
    species = [key for key in n_los_species_resolved.keys()]
    n_energy_cutoffs = len(n_los_species_resolved[species[0]])

    n_los = []
    for i in range(0, n_energy_cutoffs):
        n = 0
        for element in species:
            n += n_los_species_resolved[element][i]
        n_los.append(n)

    return n_los


def gw_basis_convergence(root: str):
    """
    Plot the quasi particle - KS gap w.r.t to the LO basis

    Have one data point for converged QP gap with l-max = (3, 2)

    Then plot for several LOs in the channels
        l-max = (4, 2)
        l-max = (3, 3)
        l-max = (5, 3)
        l-max = (4, 4)
    """
    D = OrderedDict

    # --------------------------------------------
    # Extra Zr l-Channel. l_max(Zr, O) = (4, 2)
    # zr_lmax4_o_lmax2_rgkmax8
    # --------------------------------------------
    settings_zr = {'rgkmax': 8,
                   'l_max_values': [D([('zr', 4), ('o', 2)])],
                   'n_img_freq': 32,
                   'q_grid': [2, 2, 2],
                   'n_empty_ext': [2000],
                   'max_energy_cutoffs': ['i0', 'i1', 'i2', 'i3', 'i4']
                   }

    data_set_zr = parse_gw_results(root, settings_zr)
    E_qp_zr = data_set_zr['E_qp'][:, 0] * ha_to_mev
    E_ks_zr = data_set_zr['E_ks'][:, 0] * ha_to_mev
    delta_E_qp_zr = data_set_zr['delta_E_qp'][:, 0] * ha_to_mev
    assert np.allclose(E_ks_zr, 3865.10677893), "KS gap shouldn't change"

    # Get LO labels
    basis_labels = get_basis_labels(root, settings_zr, verbose=True)
    basis_labels_zr = combine_species_basis_labels(basis_labels, species_per_line=True)['(4,2)']
    n_los_zr = n_local_orbitals(basis_labels)['(4,2)']
    print('Zr Calculation:', basis_labels_zr)


    # --------------------------------------------
    # Extra O l-channel. l_max(Zr, O) = (3, 3)
    # zr_lmax3_o_lmax3_rgkmax8
    # --------------------------------------------
    max_energy_exts_o = ['i0', 'i1', 'i2']
    settings_o = {'rgkmax': 8,
                  'l_max_values': [D([('zr', 3), ('o', 3)])],
                  'n_img_freq': 32,
                  'q_grid': [2, 2, 2],
                  'n_empty_ext': [2000],
                  'max_energy_cutoffs': max_energy_exts_o
                  }

    data_set_o = parse_gw_results(root, settings_o)
    E_qp_o = data_set_o['E_qp'][:, 0] * ha_to_mev
    E_ks_o = data_set_o['E_ks'][:, 0] * ha_to_mev
    delta_E_qp_o = data_set_o['delta_E_qp'][:, 0] * ha_to_mev
    assert np.allclose(E_ks_o, 3865.10677893), "KS gap shouldn't change"

    # Get LO labels
    basis_labels = get_basis_labels(root, settings_o, verbose=True)
    basis_labels_o = combine_species_basis_labels(basis_labels, species_per_line=True)['(3,3)']
    n_los_o = n_local_orbitals(basis_labels)['(3,3)']
    print('O Calculation:', basis_labels_o)

    # Report Details
    print('Quasi-particle gap w.r.t LOs in Zr l=4 channel (O l=2):', E_qp_zr)
    print('Quasi-particle minus gap w.r.t LOs in Zr l=4 channel (O l=2):', delta_E_qp_zr)

    print('Quasi-particle gap w.r.t LOs in O l=3 channel (Zr l=3):', E_qp_o)
    print('Quasi-particle minus gap w.r.t LOs in O l=3 channel (Zr l=3):', delta_E_qp_o)


    # ------------------------------------------------
    # Using prior Zr and O bases, increase Zr l-channel
    # l_max(Zr, O) = (5, 3)
    # zr_lmax5_o_lmax3_rgkmax8
    # ------------------------------------------------
    settings_zr = {'rgkmax': 8,
                   'l_max_values': [D([('zr', 5), ('o', 3)])],
                   'n_img_freq': 32,
                   'q_grid': [2, 2, 2],
                   'n_empty_ext': [2000],
                   'max_energy_cutoffs': ['i1', 'i2', 'i3', 'i4']  # i0 failed => 'i0',ignore it for now
                   }

    data_set_zr = parse_gw_results(root, settings_zr)
    E_qp_zr = data_set_zr['E_qp'][:, 0] * ha_to_mev
    E_ks_zr = data_set_zr['E_ks'][:, 0] * ha_to_mev
    print("KS gap has changed by 0.3 meV when moving to the (5, 3) basis\n"
          "I assume this is ok, BUT should not")
    delta_E_qp_zr = data_set_zr['delta_E_qp'][:, 0] * ha_to_mev
    assert np.allclose(E_ks_zr, 3865.3788929), "KS gap shouldn't change"

    # Get LO labels
    basis_labels = get_basis_labels(root, settings_zr, verbose=True)
    basis_labels_zr = combine_species_basis_labels(basis_labels, species_per_line=True)['(5,3)']
    n_los_zr = n_local_orbitals(basis_labels)['(5,3)']
    print('Zr Calculation (5,3):', basis_labels_zr)


    # ------------------------------------------------
    # Using prior Zr and O bases, increase O l-channel
    # l_max(Zr, O) = (4, 4)
    # zr_lmax4_o_lmax4_rgkmax8
    # ------------------------------------------------
    settings_o = {'rgkmax': 8,
                   'l_max_values': [D([('zr', 4), ('o', 4)])],
                   'n_img_freq': 32,
                   'q_grid': [2, 2, 2],
                   'n_empty_ext': [2000],
                   'max_energy_cutoffs': ['i0', 'i1', 'i2']
                   }

    data_set_o = parse_gw_results(root, settings_o)
    E_qp_o = data_set_o['E_qp'][:, 0] * ha_to_mev
    E_ks_o = data_set_o['E_ks'][:, 0] * ha_to_mev
    delta_E_qp_o = data_set_o['delta_E_qp'][:, 0] * ha_to_mev
    print(E_ks_o)
    assert np.allclose(E_ks_o, 3865.3788929), "KS gap shouldn't change"

    # Get LO labels
    basis_labels = get_basis_labels(root, settings_o, verbose=True)
    basis_labels_o = combine_species_basis_labels(basis_labels, species_per_line=True)['(4,4)']
    n_los_o = n_local_orbitals(basis_labels)['(4,4)']
    print('O Calculation:', basis_labels_o)

    print('Quasi-particle gap w.r.t LOs in Zr l=5 channel (O l=3):', E_qp_zr)
    print('Quasi-particle minus gap w.r.t LOs in Zr l=5 channel (O l=3):', delta_E_qp_zr)

    print('Quasi-particle gap w.r.t LOs in O l=4 channel (Zr l=4):', E_qp_o)
    print('Quasi-particle minus gap w.r.t LOs in O l=4 channel (Zr l=4):', delta_E_qp_o)


    # Do plot of Oxygen only, as that effect seems to be large
    # TODO Add set 6 converged data point
    # This has opened the gap a fair bit - should check tbe basis and see if another channel also affects results

    # if plot_delta_E_qp:
    #     fig, ax = plt.subplots()
    #     ax.set_xlabel('N LOs')
    #     ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)')
    #
    #     print(n_los_o)
    #     print(delta_E_qp_o)
    #
    #     ax.plot(n_los_o['o'], delta_E_qp_o, color='blue', marker='o', markersize=8)
    #
    #     for i, txt in enumerate(basis_labels_o):
    #        ax.annotate(txt, (n_los_o['o'][i], delta_E_qp_o[i]))
    #
    #     if save_plots:
    #         plt.savefig('qp_convergence_lmax43.ps', dpi=300, facecolor='w', edgecolor='w',
    #                     orientation='portrait', papertype=None, transparent=True, bbox_inches=None, pad_inches=0.1)
    #     plt.show()



    return


gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set7")
