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

    print("For basis (Zr,O) = (4,3)")
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



def gw_basis_convergence(root:str):
    """
    Plot the quasi particle - KS gap w.r.t bs LO basis, for l-max = (4, 3)
    """
    D = OrderedDict

    lmax_43 = 0
    plot_delta_E_qp = True
    plot_sigma = False
    save_plots = False

    max_energy_exts_set4_5 = ['i0', 'i1', 'i2', 'i3', 'i4', 'i5', 'i6', 'i7', 'i8', 'i9']

    settings_set4_5 = {'rgkmax': 8,
                     'l_max_values': [D([('zr', 4), ('o', 3)])],
                     'n_img_freq': 32,
                     'q_grid': [2, 2, 2],
                     'n_empty_ext': [2000],
                     'max_energy_cutoffs': max_energy_exts_set4_5
                     }

    data_set4 = parse_gw_results(root, settings_set4_5)
    delta_E_qp_set4 = data_set4['delta_E_qp'][:, lmax_43] * ha_to_mev
    basis_labels = get_basis_labels(root, settings_set4_5, verbose=True)
    basis_labels_set4 = combine_species_basis_labels(basis_labels,  species_per_line=True)['(4,3)']
    n_los_set4 = n_local_orbitals(basis_labels)['(4,3)']

    # Give some changes in QP gap w.r.t. calculations
    print(delta_E_qp_set4)
    process_basis_numbers(delta_E_qp_set4, basis_labels_set4)

    # Total number of LOs per calculation (i.e. sum Zr and O basis sizes together)
    n_los_set4= sum_los_per_species(n_los_set4)

    # --------------------------------------
    # QP vs LO basis size, for l-max = (4,3)
    # --------------------------------------
    if plot_delta_E_qp:

        fig, ax = plt.subplots()
        ax.set_xlabel('N LOs')
        ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)')

        ax.plot(n_los_set4, delta_E_qp_set4, color='blue', marker='o', linestyle='solid', linewidth=3, markersize=8)
        for i, txt in enumerate(basis_labels_set4):
           ax.annotate(txt, (n_los_set4[i], delta_E_qp_set4[i]))

        if save_plots:
            plt.savefig('qp_convergence_lmax43.ps', dpi=300, facecolor='w', edgecolor='w',
                        orientation='portrait', papertype=None, transparent=True, bbox_inches=None, pad_inches=0.1)
        plt.show()

        # Set 3 by itself, where one is primarily looking at the l=2 channel of Zr.
        # for i in range(0, delta_E_qp_set3.size):
        #     print(max_energy_exts_set3[i], delta_E_qp_set3[i])
        #
        # fig, ax = plt.subplots()
        # fig.set_size_inches(12, 10)
        # plt.rcParams.update({'font.size': 11})   # Affects fonts of the labels only
        # ax.tick_params(axis='both', which='major', labelsize=14)
        #
        # ax.set_xlim(80, 125)
        # ax.set_ylim(1580, 1660)
        # ax.set_xlabel('N LOs', fontsize=16)
        # ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)', fontsize=16)
        # ax.plot(n_los_set3, delta_E_qp_set3, 'ro', markersize=10)
        # for i, txt in enumerate(basis_labels_set3):
        #     ax.annotate(txt, (n_los_set3[i], delta_E_qp_set3[i]))

        # Note, can use plt. or fig. to save and show
        if save_plots:
            # ALso note that ps cuts off the edges, hence using jpeg
            plt.savefig('qp_convergence_set3_lmax43.jpeg', dpi=300, facecolor='w', edgecolor='w',
                        orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)

        plt.show()

    # ----------------------------------------------------------------------
    # Plot VBT and CBB real self-eneries w.r.t. LO energy cut-off
    # for set 3, only
    # ----------------------------------------------------------------------
    if plot_sigma:

        # What we expect: VBT converges with the oxygen LOs, and is well-converged
        re_self_energy_VBM = data_set4['re_self_energy_VBM'][:, lmax_43] * ha_to_mev
        fig, ax = plt.subplots()

        fig.set_size_inches(12, 10)
        plt.rcParams.update({'font.size': 11})
        ax.tick_params(axis='both', which='major', labelsize=14)

        ax.set_xlim(82, 120)
        ax.set_xlabel('N LOs', fontsize=16)
        ax.set_ylabel('Re{Self Energy} at VBM (meV)', fontsize=16)
        ax.plot(n_los_set4, re_self_energy_VBM, 'ro', markersize=10)
        for i, txt in enumerate(basis_labels_set4):
            ax.annotate(txt, (n_los_set4[i], re_self_energy_VBM[i]))

        if save_plots:
            plt.savefig('VBM_convergence_set3_lmax43.jpeg', dpi=300, facecolor='w', edgecolor='w',
                        orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)

        plt.show()

        # CBm is strongly dependent on the number of Zr LOs, and is still not converged.
        re_self_energy_CBm = data_set4['re_self_energy_CBm'][:, lmax_43] * ha_to_mev
        fig, ax = plt.subplots()

        fig.set_size_inches(12, 10)
        plt.rcParams.update({'font.size': 11})
        ax.tick_params(axis='both', which='major', labelsize=14)

        ax.set_xlim(82, 120)
        ax.set_xlabel('N LOs', fontsize=16)
        ax.set_ylabel('Re{Self Energy} at CBm (meV)', fontsize=16)
        ax.plot(n_los_set4, re_self_energy_CBm, 'ro', markersize=10)
        for i, txt in enumerate(basis_labels_set4):
            ax.annotate(txt, (n_los_set4[i], re_self_energy_CBm[i]))

        if save_plots:
            plt.savefig('CBm_convergence_set3_lmax43.jpeg', dpi=300, facecolor='w', edgecolor='w',
                        orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)

        plt.show()

    return

gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set4_5")



# def gw_q_convergence(root:str):
#     """
#     TODO Clean up
#
#     :param root:
#     """
#
#     # max LO energy
#     E_LO = 140
#     # Uniform q-grids
#     # q_points = [[2, 2, 2], [3, 3, 3], [4, 4, 4]]
#     q_points = [[2, 2, 2], [4, 4, 4]]
#
#     # ---------------
#     # Parse
#     # ---------------
#     delta_E_qp = []
#     for q in q_points:
#         q_str = "".join(str(entry) for entry in q)
#         directory = root + '/q' + q_str + '_max_energy_' + str(E_LO)
#         print(directory)
#         gw_data = parse_gw_info(directory)
#         qp_data = parse_gw_evalqp(directory, nkpts=3)
#         results = process_gw_gamma_point(gw_data, qp_data)
#         delta_E_qp.append(results['E_qp'] - results['E_ks'])
#
#     # ---------------
#     # Plot
#     # ---------------
#     q = np.arange(0, len(q_points))
#     q_labels = []
#     for q_point in q_points:
#         q_labels.append("(" + ",".join(str(entry) for entry in q_point) + ")")
#
#     fig, ax = plt.subplots()
#     plt.xlabel("q-grid")
#     plt.ylabel("Quasiparticle Gap - KS Gap at Gamma (meV)")
#
#     ax.set_xticks(q)
#     ax.set_xticklabels(q_labels)
#
#     plt.plot(q, np.asarray(delta_E_qp) * ha_to_mev, 'bo--', linewidth=2, markersize=12)
#     plt.show()
#
#     # ---------------
#     # Process
#     # ---------------
#     for i in range(1, len(q_points)):
#         change_in_delta_E_qp = (delta_E_qp[i] - delta_E_qp[i-1]) * ha_to_mev
#         print("Change in delta_E_qp from q-grid", q_points[i-1], "to", q_points[i],
#               ":", change_in_delta_E_qp, 'meV')
#
#     return
#
#gw_q_convergence("/users/sol/abuccheri/gw_benchmarks/A1/zr_lmax4_o_lmax3_rgkmax7/q_convergence")

