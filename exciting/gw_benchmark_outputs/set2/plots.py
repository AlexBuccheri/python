"""
Generate Plots and Post-process Data
"""

# External libraries
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

# My modules
from units_and_constants.unit_conversions import ha_to_mev
from ex_plot.plot import Plot
from gw_benchmark_outputs.post_process_utils import parse_gw_results, get_basis_labels, combine_species_basis_labels


def process_basis_numbers(l_max_values, delta_E_qp, max_energy_exts:list, max_energy):
    """
    Print information of interest
    """
    # For max energy, print the change in E_QP w.r.t. the l_max of the basis
    print('For max energy parameter cutoff = ' +str(max_energy) + ' Ha')
    ienergy = max_energy_exts.index(max_energy)
    prior_lmax = [3, 2]

    for i in range(1, len(l_max_values)):
        current_lmax = [l + 1 for l in prior_lmax]
        change_in_delta_E_qp = (delta_E_qp[ienergy, i] - delta_E_qp[ienergy, i-1]) * ha_to_mev
        print("Change in Delta E_QP from basis ", prior_lmax, "to ", current_lmax, ":", change_in_delta_E_qp, "(meV)")
        prior_lmax = current_lmax

    # Print change in QP gap for each max energy parameter (i.e. max LO)
    print("For basis (Zr,O) = (4,3)")
    lmaxs_to_index = {'(3,2)':0, '(4,3)':1, '(5,4)':2, '(6,5)':3}
    lmax_index = lmaxs_to_index['(4,3)']

    for ie, energy in enumerate(max_energy_exts[1:], start=1):
        change_in_delta_E_qp = (delta_E_qp[ie, lmax_index] - delta_E_qp[ie-1, lmax_index]) * ha_to_mev
        print("Change in Delta E_QP from max energy param ", max_energy_exts[ie-1], "to ", energy, ":", change_in_delta_E_qp, "(meV)")

    return


def gw_basis_convergence(root:str):
    """
    Plot the quasipartilce-KS gap w.r.t basis l-max, for
    several max lo energies.
    """
    D = OrderedDict

    # Details required for generating directory extension
    l_max_values = [D([('zr', 3), ('o', 2)]),
                    D([('zr', 4), ('o', 3)]),
                    D([('zr', 5), ('o', 4)]),
                    D([('zr', 6), ('o', 5)])]

    # Makes more sense to specifiy the extensions rather than use:
    # max_energy_exts = max_energy_ext_per_directory(energy_cutoffs)
    # Same with n_empty_ext
    max_energy_exts = [90, 160, 230, 300]

    settings = {'rgkmax': 7,
                'l_max_values': l_max_values,
                'n_img_freq': 32,
                'q_grid': [2, 2, 2],
                'n_empty_ext': [800, 1000, 1300, 2000],
                'max_energy_cutoffs': max_energy_exts}

    # Parse data
    data = parse_gw_results(root, settings)
    delta_E_qp = data['delta_E_qp']

    # Print some data
    process_basis_numbers(l_max_values, delta_E_qp, max_energy_exts, max_energy = 300)

    # Generate plots
    save_plots = False
    plot_delta_E_qp = True
    plot_sigma = False


    # -----------------
    # Plot delta_E_qp
    # -----------------
    if plot_delta_E_qp:

        # --------------------------------------------------------------
        # Plot 1
        # --------------------------------------------------------------
        x = np.arange(0, len(l_max_values))
        x_labels = ["(" + ",".join(str(l) for l in l_pair.values()) + ")" for l_pair in l_max_values]

        gw_plot = Plot(x, [], x_label='l_max (Zr, O)', y_label='Quasiparticle Gap - KS Gap at Gamma (meV)',
                       xticklabels=x_labels, legend_title='Max LO Energy (Ha)', linewidth=3)

        for i, energy in enumerate(max_energy_exts):
            gw_plot.plot_data(x, delta_E_qp[i, :] * ha_to_mev, label=str(energy))

        if save_plots:
            save_options = {'file_name':'qp_convergence_lmax_and_los.jpeg',
                            'dpi':300,
                            'facecolor':'w',
                            'edgecolor':'w',
                            'orientation':'landscape',
                            'bbox_inches':'tight'}
            gw_plot.save(save_options)

        gw_plot.show()

        # --------------------------------------------------------------
        # Second plot. Same as above but label every point with the basis
        # --------------------------------------------------------------
        fig, ax = plt.subplots()
        fig.set_size_inches(14, 10)
        plt.rcParams.update({'font.size': 16})

        ax.set_xlabel('l_max (Zr, O)', fontsize=16)
        ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)', fontsize=16)

        ax.set_xticks(x)
        ax.set_xticklabels(x_labels)
        ax.tick_params(axis='both', which='major', labelsize=16)

        # Overlay lines for each LO energy cutoff
        for i, energy in enumerate(max_energy_exts):
            ax.plot(x, delta_E_qp[i, :] * ha_to_mev, marker='o', markersize=8)

        basis_labels = get_basis_labels(root, settings)

        # Set basis label for each marker
        for j, l_max in enumerate(['(3,2)', '(4,3)', '(5,4)', '(6,5)']):
            for i, energy in enumerate(max_energy_exts):
                label = ''
                for species in ['zr', 'o']:
                    label += species.capitalize() + ': ' + basis_labels[l_max][species][i].rstrip()+'.\n'
                ax.annotate(label.rstrip(), (x[j],  delta_E_qp[i, j] * ha_to_mev))

        if save_plots:
            plt.savefig('qp_convergence_set2_basis_labels.jpeg', dpi=300, facecolor='w', edgecolor='w',
                        orientation='landscape', transparent=True, bbox_inches='tight', pad_inches=1.)

        plt.show()



    # ----------------------------------------------------------------------
    # Plot VBT and CBB real self-eneries w.r.t. basis and energy parameter
    # ----------------------------------------------------------------------
    if plot_sigma:
        re_self_energy_VBM = data['re_self_energy_VBM']
        r_sigma_VBM_plot = Plot(x, [], x_label='l_max (Zr, O)', y_label='Re{Self Energy} at VBM (meV) ',
                                xticklabels=x_labels, legend_title='Max LO Energy Param (Ha)', linewidth=3)

        for i, energy in enumerate(max_energy_exts):
            r_sigma_VBM_plot.plot_data(x, re_self_energy_VBM[i, :] * ha_to_mev, label=str(energy))
        r_sigma_VBM_plot.show()

        re_self_energy_CBm = data['re_self_energy_CBm']
        r_sigma_CBm_plot = Plot(x, [], x_label='l_max (Zr, O)', y_label='Re{Self Energy} at CBm (meV) ',
                                xticklabels=x_labels, legend_title='Max LO Energy Param (Ha)', linewidth=3)

        for i, energy in enumerate(max_energy_exts):
            r_sigma_CBm_plot.plot_data(x, re_self_energy_CBm[i, :] * ha_to_mev, label=str(energy))
        r_sigma_CBm_plot.show()

    return

gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_var_cutoff")

