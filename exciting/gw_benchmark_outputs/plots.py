"""
Plots
"""
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gamma_point


def directory_string(l_max:OrderedDict):
    path = ''
    for species, l in l_max.items():
        path += species + '_lmax' + str(l) + '_'

    return path


def gw_basis_convergence(root:str):
    """
    Plot the quasipartilce-KS gap w.r.t basis l-max, for
    several max lo energies
    """
    D = OrderedDict

    # Calculation settings
    l_max_values = [D([('zr', 3), ('o', 2)]),
                    D([('zr', 4), ('o', 3)]),
                    D([('zr', 5), ('o', 4)]),
                    D([('zr', 6), ('o', 5)])]

    n_img_freq = 32
    q_grid = [2, 2, 2]
    n_empty = [450, 500, 550, 600]
    #max_energy_cutoff = [60, 80, 100, 120, 140, 160]
    max_energy_cutoff = [60]

    delta_E_qp = np.empty(shape=(len(max_energy_cutoff), len(l_max_values)))

    for i, l_max in enumerate(l_max_values):
        basis_root = root + '/' + directory_string(l_max) + 'rgkmax7'
        q_str = "".join(str(q) for q in q_grid)
        gw_root = basis_root + "/gw_q" + q_str + "_omeg" + str(n_img_freq) + "_nempty" + str(n_empty[i])

        for energy in max_energy_cutoff:
            file_path = gw_root + '/max_energy_' + str(energy)
            gw_data = parse_gw_info(file_path)
            qp_data = parse_gw_evalqp(file_path, nempty=n_empty[i], nkpts=3)
            results = process_gw_gamma_point(gw_data, qp_data)

            delta_E_qp[ienergy, i] = results['E_qp'] - results['E_KS']




    # Produce the plots


    x = np.arange(0, len(l_max_values))
    fig, ax = plt.subplots()
    plt.xlabel("l_max (Zr, O)")
    plt.ylabel("Quasiparticle Gap - KS Gap at Gamma (eV)")

    for i, energy in enumerate(max_energy_cutoff):
        plt.plot(x, delta_E_qp[i, :], label=str(energy))

    ax.set_xticks(x)
    #ax.set_xticklabels(x_labels)
    legend = ax.legend(loc='upper right', title='Energy Param Cutoff (Ha)')
    # plt.savefig('Re_sigma_VBM.pdf', dpi=300)
    plt.show()

    return

gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1")



