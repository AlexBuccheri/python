"""
Plots
"""
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gamma_point
from units_and_constants.unit_conversions import ha_to_mev


def directory_string(l_max:OrderedDict) -> str:
    """
    Top level basis convergence directory string
    """
    path = ''
    for species, l in l_max.items():
        path += species + '_lmax' + str(l) + '_'

    return path

# TODO(Alex) Duplicate of routine in all A1_zr_lmax_o_lmax scripts
def restructure_energy_cutoffs(energy_cutoffs: dict) -> list:
    """
    Get in a more useful structure to iterate over

    """
    restructured_energies = []

    for inum in range(0, 4):
        data = {}
        for species, l_channels in  energy_cutoffs.items():
            data[species] = {l:energies[inum] for l, energies in l_channels.items()}
        restructured_energies.append(data)

    return restructured_energies


def max_energy_ext_per_directory(energy_cutoffs):
    """
    Get the energy extension string for each run directory,
    for a given set of calculations

    :param energy_cutoffs: of the form

    {'zr': {0: np.linspace(60, 120, num=4),
            1: np.linspace(60, 120, num=4),
            2: np.linspace(90, 300, num=4),
            3: np.linspace(60, 120, num=4)},
     'o':  {0: np.linspace(60, 120, num=4),
            1: np.linspace(60, 120, num=4),
            2: np.linspace(60, 120, num=4)}
                       }
    :return: max_energy_exts
    """
    max_energy_exts = []
    for energy in restructure_energy_cutoffs(energy_cutoffs):
        max_energy_per_species = [max(energy_per_l_channel.values()) for energy_per_l_channel in
                                  energy.values()]
        max_energy_exts.append( str(int(max(max_energy_per_species))) )
    return max_energy_exts


def parse_gw_results(root:str, settings:dict) -> dict:
    """

    QP direct-gap (relative to the KS gap)
    and self-energies of the band edges at Gamma.

    Almost easier to just path full directive strings.

    :return: dictionary containing the above.
    """

    # Basis settings
    rgkmax = settings['rgkmax']
    l_max_values = settings['l_max_values']
    max_energy_exts = settings['max_energy_cutoffs']

    # GW settings
    n_empty_ext = settings['n_empty_ext']
    q_grid = settings['q_grid']
    n_img_freq = settings['n_img_freq']

    # Data to parse and return
    delta_E_qp = np.empty(shape=(len(max_energy_exts), len(l_max_values)))
    re_self_energy_VBM = np.empty(shape=delta_E_qp.shape)
    re_self_energy_CBm = np.empty(shape=delta_E_qp.shape)
    q_str = "".join(str(q) for q in q_grid)

    # Lmax in LO basis
    for i, l_max in enumerate(l_max_values):
        basis_root = root + '/' + directory_string(l_max) + 'rgkmax' + str(rgkmax)
        gw_root = basis_root + "/gw_q" + q_str + "_omeg" + str(n_img_freq) + "_nempty" + str(n_empty_ext[i])

        # Max energy cut-off of LOs in each l-channel
        for ienergy, energy in enumerate(max_energy_exts):
            file_path = gw_root + '/max_energy_' + str(energy)

            gw_data = parse_gw_info(file_path)
            qp_data = parse_gw_evalqp(file_path)

            print('Reading data from ', file_path)
            results = process_gw_gamma_point(gw_data, qp_data)
            delta_E_qp[ienergy, i] = results['E_qp'] - results['E_ks']
            re_self_energy_VBM[ienergy, i] = results['re_sigma_VBM']
            re_self_energy_CBm[ienergy, i] = results['re_sigma_CBm']

    return {'delta_E_qp': delta_E_qp,
            're_self_energy_VBM': re_self_energy_VBM,
            're_self_energy_CBm': re_self_energy_CBm
            }



class PlotGW:

    def __init__(self, x, y, x_label='', y_label='', xticklabels=None, yticklabels=None, legend_loc='upper right', legend_title=None):
        self.x = x
        self.y = y
        self.x_label = x_label
        self.y_label = y_label
        self.xticklabels = xticklabels
        self.yticklabels = yticklabels
        self.legend_loc = legend_loc
        self.legend_title = legend_title
        self.initialise_plot()

    def initialise_plot(self):
        """
        Initialise plot axes and labels
        :return:
        """
        self.fig , self.ax = plt.subplots()
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel(self.y_label)

        if self.xticklabels is not None:
            self.ax.set_xticks(self.x)
            self.ax.set_xticklabels(self.xticklabels)

        if self.yticklabels is not None:
            self.ax.set_yticks(self.y)
            self.ax.set_yticklabels(self.yticklabels)

        return

    def plot_data(self, x, y, label=None):
        """
        Basic plot wrapper
        :param x: x data
        :param y: y data
        :param label: optional label for legend
        """
        if label is None:
            self.ax.plot(x, y)
        else:
            self.ax.plot(x, y, label=str(label))
        return

    def show(self, file_name=None):
        if file_name is not None:
            self.ax.legend(loc=self.legend_loc, title=self.legend_title)
            plt.savefig(file_name, dpi=300)
        else:
            self.ax.legend(loc=self.legend_loc, title=self.legend_title)
            plt.show()
        return



def process_basis_numbers(l_max_values, delta_E_qp, max_energy_exts:list, max_energy):
    """
    Print information of interest
    """

    # For max energy, print the change in E_QP w.r.t. the l_max of the basis
    ienergy = max_energy_exts.index(str(max_energy))
    prior_lmax = [3, 2]

    print('For max energy parameter cutoff = ' +str(max_energy) + ' Ha')
    for i in range(1, len(l_max_values)):
        current_lmax = [l + 1 for l in prior_lmax]
        change_in_delta_E_qp = (delta_E_qp[ienergy, i] - delta_E_qp[ienergy, i-1]) * ha_to_mev

        print("Change in Delta E_QP from basis ", prior_lmax, "to ", current_lmax, ":", change_in_delta_E_qp, "(meV)")
        prior_lmax = current_lmax

    # For (4,3), change in energy for each energy parameter
    lo_43 = 1
    print("For basis (Zr,O) = (4,3)")
    for ie, energy in enumerate(max_energy_exts[1:], start=1):
        change_in_delta_E_qp = (delta_E_qp[ie, lo_43] - delta_E_qp[ie-1, lo_43]) * ha_to_mev
        print("Change in Delta E_QP from max energy param ", max_energy_exts[ie-1], "to ", energy, ":", change_in_delta_E_qp, "(meV)")

    return


def gw_basis_convergence(root:str):
    """
    Plot the quasipartilce-KS gap w.r.t basis l-max, for
    several max lo energies.
    NOTE: n_empty_ext corresponds to the directory extensions
    """
    D = OrderedDict

    # Details required for generating directory extension
    l_max_values = [D([('zr', 3), ('o', 2)]),
                    D([('zr', 4), ('o', 3)]),
                    D([('zr', 5), ('o', 4)]),
                    D([('zr', 6), ('o', 5)])]

    # Makes more sense to just specifiy the extensions, rather than generate them
    # in the same way I've done with n_empty
    # max_energy_exts = max_energy_ext_per_directory(energy_cutoffs)
    max_energy_exts = [60, 80, 100, 120, 140, 160, 180, 200]

    data = parse_gw_results(root, {'rgkmax': 7,
                                   'l_max_values': l_max_values,
                                   'n_img_freq': 32,
                                   'q_grid': [2, 2, 2],
                                   'n_empty_ext': [800, 1000, 1300, 2000],
                                   'max_energy_cutoffs': max_energy_exts})

    delta_E_qp = data['delta_E_qp']
    process_basis_numbers(l_max_values, delta_E_qp, max_energy_exts, max_energy = 300)

    plot_delta_E_qp = True
    plot_sigma = True

    # -----------------
    # Plot delta_E_qp
    # -----------------
    if plot_delta_E_qp:
        x = np.arange(0, len(l_max_values))
        x_labels = ["(" + ",".join(str(l) for l in l_pair.values()) + ")" for l_pair in l_max_values]

        gw_plot = PlotGW(x, [], x_label='l_max (Zr, O)', y_label='Quasiparticle Gap - KS Gap at Gamma (meV)',
                         xticklabels=x_labels, legend_title='Max LO Energy Param (Ha)')

        for i, energy in enumerate(max_energy_exts):
            gw_plot.plot_data(x, delta_E_qp[i, :] * ha_to_mev, label=str(energy))

        gw_plot.show()


    # ----------------------------------------------------------------------
    # Plot VBT and CBB real self-eneries w.r.t. basis and energy parameter
    # ----------------------------------------------------------------------
    if plot_sigma:
        re_self_energy_VBM = data['re_self_energy_VBM']
        r_sigma_VBM_plot = PlotGW(x, [], x_label='l_max (Zr, O)', y_label='Re{Self Energy} at VBM (meV) ',
                                  xticklabels=x_labels, legend_title='Max LO Energy Param (Ha)')

        for i, energy in enumerate(max_energy_exts):
            r_sigma_VBM_plot.plot_data(x, re_self_energy_VBM[i, :] * ha_to_mev, label=str(energy))
        r_sigma_VBM_plot.show()

        re_self_energy_CBm = data['re_self_energy_CBm']
        r_sigma_CBm_plot = PlotGW(x, [], x_label='l_max (Zr, O)', y_label='Re{Self Energy} at CBm (meV) ',
                                  xticklabels=x_labels, legend_title='Max LO Energy Param (Ha)')

        for i, energy in enumerate(max_energy_exts):
            r_sigma_CBm_plot.plot_data(x, re_self_energy_CBm[i, :] * ha_to_mev, label=str(energy))
        r_sigma_CBm_plot.show()


    return

gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_var_cutoff")



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

