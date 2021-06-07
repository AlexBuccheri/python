"""
Plots
"""
# External libraries
from collections import OrderedDict
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
    Plot the quasi particle - KS gap w.r.t bs LO basis, for l-max = (4, 3)
    """
    D = OrderedDict

    info = """  Converge QP gap w.r.t the s-channel LOs of Zr
    """

    print(info)

    lmax_32 = 0
    plot_delta_E_qp = True
    plot_sigma = False
    save_plots = False

    max_energy_exts_set4_spd = ['i0', 'i1', 'i2', 'i3']  # 'i4'

    settings_set4_spd = {'rgkmax': 8,
                         'l_max_values': [D([('zr', 3), ('o', 2)])],
                         'n_img_freq': 32,
                         'q_grid': [2, 2, 2],
                         'n_empty_ext': [2000],
                         'max_energy_cutoffs': max_energy_exts_set4_spd
                        }

    # S out of core
    data_set4 = parse_gw_results(root, settings_set4_spd)
    delta_E_qp_set4 = data_set4['delta_E_qp'][:, lmax_32] * ha_to_mev
    E_qp_set4 = data_set4['E_qp'][:, lmax_32] * ha_to_mev
    E_ks_set4 = data_set4['E_ks'][:, lmax_32] * ha_to_mev


    basis_labels = get_basis_labels(root, settings_set4_spd, verbose=True)
    basis_labels_set4 = combine_species_basis_labels(basis_labels, species_per_line=True)['(3,2)']
    n_los_set4 = n_local_orbitals(basis_labels)['(3,2)']

    # Give some changes in QP gap w.r.t. calculations
    print('s out of core')
    print('KS:', E_ks_set4)
    print('QP:', E_qp_set4)
    print('QP-KS:', delta_E_qp_set4)
    process_basis_numbers(delta_E_qp_set4, basis_labels_set4)

    # Total number of LOs per calculation (i.e. sum Zr and O basis sizes together)
    n_los_set4 = sum_los_per_species(n_los_set4)

    # S in core
    data_set4_in_core = parse_gw_results(root, settings_set4_spd, dir_prefix='s_in_core_')
    delta_E_qp_s_core = data_set4_in_core['delta_E_qp'][:, lmax_32] * ha_to_mev
    E_qp_set4_s_core = data_set4_in_core['E_qp'][:, lmax_32] * ha_to_mev
    E_ks_set4_s_core = data_set4_in_core['E_ks'][:, lmax_32] * ha_to_mev

    print('s in core')
    print('KS:', E_ks_set4_s_core)
    print('QP:', E_qp_set4_s_core)
    print('QP-KS:', delta_E_qp_s_core)
    quit()


    # --------------------------------------
    # QP vs LO basis size, for l-max = (4,3)
    # --------------------------------------
    if plot_delta_E_qp:

        fig, ax = plt.subplots()
        ax.set_xlabel('N LOs')
        ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)')

        ax.plot(n_los_set4, delta_E_qp_set4, color='blue', marker='o', markersize=8)
        for i, txt in enumerate(basis_labels_set4):
           ax.annotate(txt, (n_los_set4[i], delta_E_qp_set4[i]))

        if save_plots:
            plt.savefig('qp_convergence_lmax43.ps', dpi=300, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, transparent=True, bbox_inches=None, pad_inches=0.1)
        plt.show()


    # ----------------------------------------------------------------------
    # Plot VBT and CBB real self-energies w.r.t. LO energy cut-off
    # for set 3, only
    # ----------------------------------------------------------------------
    if plot_sigma:
        print("Requires rewriting")
        pass
        # # What we expect: VBT converges with the oxygen LOs, and is well-converged
        # re_self_energy_VBM = data_set4['re_self_energy_VBM'][:, lmax_32] * ha_to_mev
        # fig, ax = plt.subplots()
        #
        # fig.set_size_inches(12, 10)
        # plt.rcParams.update({'font.size': 11})
        # ax.tick_params(axis='both', which='major', labelsize=14)
        #
        # ax.set_xlim(82, 120)
        # ax.set_xlabel('N LOs', fontsize=16)
        # ax.set_ylabel('Re{Self Energy} at VBM (meV)', fontsize=16)
        # ax.plot(n_los_set4, re_self_energy_VBM, 'ro', markersize=10)
        # for i, txt in enumerate(basis_labels_set4):
        #     ax.annotate(txt, (n_los_set4[i], re_self_energy_VBM[i]))
        #
        # if save_plots:
        #     plt.savefig('VBM_convergence_set3_lmax43.jpeg', dpi=300, facecolor='w', edgecolor='w',
        #                 orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
        #
        # plt.show()
        #
        # # CBm is strongly dependent on the number of Zr LOs, and is still not converged.
        # re_self_energy_CBm = data_set4['re_self_energy_CBm'][:, lmax_32] * ha_to_mev
        # fig, ax = plt.subplots()
        #
        # fig.set_size_inches(12, 10)
        # plt.rcParams.update({'font.size': 11})
        # ax.tick_params(axis='both', which='major', labelsize=14)
        #
        # ax.set_xlim(82, 120)
        # ax.set_xlabel('N LOs', fontsize=16)
        # ax.set_ylabel('Re{Self Energy} at CBm (meV)', fontsize=16)
        # ax.plot(n_los_set4, re_self_energy_CBm, 'ro', markersize=10)
        # for i, txt in enumerate(basis_labels_set4):
        #     ax.annotate(txt, (n_los_set4[i], re_self_energy_CBm[i]))
        #
        # if save_plots:
        #     plt.savefig('CBm_convergence_set3_lmax43.jpeg', dpi=300, facecolor='w', edgecolor='w',
        #                 orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
        #
        # plt.show()

    return


gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set4_spd")
