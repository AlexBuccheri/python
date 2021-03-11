import matplotlib.pyplot as plt
import numpy as np

ha_to_ev = 27.2114

# q = 2x2x2. All unoccupied states. Note other settings
# keys = l_max(Zr, O)
# Waiting on (6,5) results
# For tabulating KS values, I looked at 20 Ha and 100 Ha basis cut-off results. If they were consistent, I didn't check the intermediate calcs
# Ultimately, need to move to parsing these values

gw_results = {(2,1): {'trial_energy': [20, 40, 60, 80, 100],
                      'E_GW_VBT': np.array([-0.10974, -0.10987, -0.11013, -0.11000, -0.11007]),
                      'real_sigma_VBT': np.array([0.16228, 0.16211, 0.16177, 0.16194, 0.16186]),
                      'E_GW_CBB': np.array([0.10079, 0.09905, 0.09829, 0.09660, 0.09508]),
                      'real_sigma_CBB': np.array([-0.21359, -0.21588, -0.21690, -0.21907, -0.22106]),
                      'E_ks_VBT': np.array([-0.08109, -0.08109, -0.08109, -0.08109, -0.08109]),
                      'E_ks_CBB': np.array([0.06091, 0.06091, 0.06091, 0.06091, 0.06091])
                      },

              (3,2):  {'trial_energy'  : [20, 40, 60, 80, 100],
                       'E_GW_VBT'      : np.array([-0.11878, -0.11984, -0.11985, -0.12007, -0.11999]),
                       'real_sigma_VBT': np.array([0.15031, 0.14894, 0.14890, 0.14862, 0.14873]),
                       'E_GW_CBB'      : np.array([0.06091, 0.09519, 0.09539,  0.09368, 0.09120]),
                       'real_sigma_CBB': np.array([-0.21767, -0.22067, -0.22041, -0.22265, -0.22586]),
                       'E_ks_VBT': np.array([-0.08109, -0.08109, -0.08109, -0.08109, -0.08109]),
                       'E_ks_CBB': np.array([0.06091, 0.06091, 0.06091, 0.06091, 0.06091])
                       },

              (4,3):  {'trial_energy'  : [20, 40, 60, 80, 100],
                       'E_GW_VBT'      : np.array([-0.12295, -0.12388, -0.12397, -0.12394, -0.12407]),
                       'real_sigma_VBT': np.array([0.14478, 0.14356, 0.14344, 0.14347, 0.14333]),
                       'E_GW_CBB'      : np.array([0.09482, 0.09397, 0.09421, 0.09369, 0.08954]),
                       'real_sigma_CBB': np.array([-0.22102, -0.22206, -0.22175, -0.22243, -0.22783]),
                       'E_ks_VBT': np.array([-0.08108, -0.08108, -0.08108, -0.08108, -0.08108]),
                       'E_ks_CBB': np.array([0.06091, 0.06091, 0.06091, 0.06091, 0.06091])
                       },

              (5,4):  {'trial_energy'  : [20, 40, 60, 80, 100],
                       'E_GW_VBT'      : np.array([-0.12564, -0.12641, -0.12656, -0.12647, -0.12663]),
                       'real_sigma_VBT': np.array([0.14123, 0.14019, 0.14001, 0.14011, 0.13995]),
                       'E_GW_CBB'      : np.array([0.09236, 0.09152, 0.09184, 0.09105, 0.08663]),
                       'real_sigma_CBB': np.array([-0.22408, -0.22512, -0.22472, -0.22575, -0.23147]),
                       'E_ks_VBT': np.array([-0.08108, -0.08108, -0.08108, -0.08108, -0.08108]),
                       'E_ks_CBB': np.array([0.06091, 0.06091, 0.06091, 0.06091, 0.06091])
                       },

             # Did I change the basis? Apparently state 12 is still the high-occupied but the KS degeneracy
              # differs in the GW output
              # (6, 5): {'trial_energy': [20, 40, 60, 80, 100],
              #          'E_GW_VBT': np.array([0.08145]),
              #          'real_sigma_VBT': np.array([]),
              #          'E_GW_CBB': np.array([-0.05516]),
              #          'real_sigma_CBB': np.array([]),
              #          'E_ks_VBT': np.array([-0.02697]),
              #          'E_ks_CBB': np.array([])
              #          }

              }


trial_energy_to_index = {20:  0,
                         40:  1,
                         60:  2,
                         80:  3,
                         100: 4}


def quasiparticle_gap_wrt_lmax(gw_results, trial_energy, fig=None, ax=None):

    if fig is None and ax is None:
        # Initialise figure and axes
        fig, ax = plt.subplots(figsize=(6, 9))


    i_energy = trial_energy_to_index[trial_energy]
    x_values = np.arange(0, len(gw_results.values()))

    x_labels = []
    y_values = []

    for l_max, value in gw_results.items():
        Eg_GW = (value['E_GW_CBB'][i_energy] - value['E_GW_VBT'][i_energy]) * ha_to_ev
        x_labels.append(str(l_max))
        y_values.append(Eg_GW)

    plt.xlabel("(Zr l_max, O l_max)")
    plt.ylabel("Quasiparticle Gap at Gamma (eV)")

    plt.plot(x_values, y_values, label=str(trial_energy))
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_labels)

    return fig, ax


def delta_quasiparticle_gap_wrt_lmax(gw_results, trial_energy, fig=None, ax=None):

    if fig is None and ax is None:
        # Initialise figure and axes
        fig, ax = plt.subplots(figsize=(6, 9))

    i_energy = trial_energy_to_index[trial_energy]
    x_values = np.arange(0, len(gw_results.values()))

    Eg_GW_least_converged = \
        (gw_results[(2, 1)]['E_GW_CBB'][i_energy] - gw_results[(2, 1)]['E_GW_VBT'][i_energy]) * ha_to_ev * 1000.

    x_labels = []
    y_values = []

    for l_max, value in gw_results.items():
        Eg_GW = (value['E_GW_CBB'][i_energy] - value['E_GW_VBT'][i_energy]) * ha_to_ev * 1000.
        x_labels.append(str(l_max))
        y_values.append(Eg_GW - Eg_GW_least_converged)

    plt.xlabel("(Zr l_max, O l_max)")
    plt.ylabel("Delta Quasiparticle Gap at Gamma (meV)")

    plt.plot(x_values, y_values, 'ro-', label=str(trial_energy))
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_labels)

    return fig, ax


def quasiparticle_gap_wrt_lmax_and_trial_energies(gw_results, fig, ax, trial_energies):
    for trial_energy in trial_energies:
        fig, ax = quasiparticle_gap_wrt_lmax(gw_results, trial_energy, fig=fig, ax=ax)
    return fig, ax


def quasiparticle_gap_minus_ks_gap(gw_results, trial_energy, fig=None, ax=None):

    if fig is None and ax is None:
        fig, ax = plt.subplots()

    i_energy = trial_energy_to_index[trial_energy]

    x_values = np.arange(0, len(gw_results.values()))
    x_labels = []
    y_values = []

    for l_max, value in gw_results.items():
        Eg_GW = (value['E_GW_CBB'][i_energy] - value['E_GW_VBT'][i_energy]) * ha_to_ev
        KS_gap = (value['E_ks_CBB'][i_energy] - value['E_ks_VBT'][i_energy]) * ha_to_ev
        x_labels.append(str(l_max))
        y_values.append(Eg_GW - KS_gap)

    plt.xlabel("(Zr l_max, O l_max)")
    plt.ylabel("Quasiparticle Gap - KS Gap at Gamma (eV)")

    plt.plot(x_values, y_values, label=str(trial_energy))
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_labels)

    return fig, ax




def real_sigma_cbb_wrt_lmax(gw_results, trial_energy, fig=None, ax=None):

    if fig is None and ax is None:
        # Initialise figure and axes
        fig, ax = plt.subplots(figsize=(6, 9))

    i_energy = trial_energy_to_index[trial_energy]
    x_values = np.arange(0, len(gw_results.values()))

    x_labels = []
    y_values = []

    for l_max, value in gw_results.items():
        re_sigma_cbb = value['real_sigma_CBB'][i_energy] * ha_to_ev
        x_labels.append(str(l_max))
        y_values.append(re_sigma_cbb)

    plt.xlabel("(Zr l_max, O l_max)")
    plt.ylabel("Re(Self-energy) CBm at Gamma (eV)")

    plt.plot(x_values, y_values, label=str(trial_energy))
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_labels)

    return fig, ax


def real_sigma_cbb_lmax_and_trial_energies(gw_results, fig, ax, trial_energies):
    for trial_energy in trial_energies:
        fig, ax = real_sigma_cbb_wrt_lmax(gw_results, trial_energy, fig=fig, ax=ax)
    return fig, ax


def real_sigma_vbt_wrt_lmax(gw_results, trial_energy, fig=None, ax=None):

    if fig is None and ax is None:
        # Initialise figure and axes
        fig, ax = plt.subplots(figsize=(6, 9))

    i_energy = trial_energy_to_index[trial_energy]
    x_values = np.arange(0, len(gw_results.values()))

    x_labels = []
    y_values = []

    for l_max, value in gw_results.items():
        re_sigma_vbt = value['real_sigma_VBT'][i_energy] * ha_to_ev
        x_labels.append(str(l_max))
        y_values.append(re_sigma_vbt)

    plt.xlabel("(Zr l_max, O l_max)")
    plt.ylabel("Re(Self-energy) VBM at Gamma (eV)")

    plt.plot(x_values, y_values, label=str(trial_energy))
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_labels)

    return fig, ax


def real_sigma_vbt_lmax_and_trial_energies(gw_results, fig, ax, trial_energies):
    for trial_energy in trial_energies:
        fig, ax = real_sigma_vbt_wrt_lmax(gw_results, trial_energy, fig=fig, ax=ax)
    return fig, ax


