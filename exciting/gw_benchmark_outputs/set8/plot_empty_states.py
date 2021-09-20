import os
from typing import List
import matplotlib.pyplot as plt

from units_and_constants.unit_conversions import ha_to_mev

from plots import process_gw_calculation


def get_empty_states_vs_qp(root: str, nempty_states: List[int]):
    """
    Change in QP energy as a function of empty states.
    Converged basis and q=2x2x2

    Expect data[i] to be a dict with keys:
    'E_qp' 'E_ks' 'E_qp_X_Gamma' 'E_ks_X_Gamma' 'E_qp_X_X' 'E_ks_X_X'
    'delta_E_qp' 're_self_energy_VBM' 're_self_energy_CBm'
    """
    data = []
    for nempty in nempty_states:
        full_path = os.path.join(root, str(nempty))
        data.append(process_gw_calculation(full_path))
    return data


def plot_empty_states_vs_qp(nempty_states: List[int], data: List[dict]):
    """

    :param List[int] nempty_states:
    :param List[dict] data:
    """
    # Direct gap as a function of nempty states
    E_GG = [entry['E_qp'] * ha_to_mev for entry in data]

    fig, ax = plt.subplots()
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})

    ax.set_xlabel('Nempty', fontsize=16)
    ax.set_ylabel('Direct Quasiparticle Gap (meV)', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.plot(nempty_states, E_GG, marker='o', markersize=8)
    plt.show()


if __name__ == "__main__":
    l_max_pairs = {'zr': 6, 'o': 5}
    # Where 1245 = maximum number of states for this basis
    nempty_states = [250, 500, 750, 1000, 1100, 1245]
    root = "/users/sol/abuccheri/gw_benchmarks/A1_set8/zr_lmax6_o_lmax5_rgkmax8/nempty_222"
    data = get_empty_states_vs_qp(root, nempty_states)
    plot_empty_states_vs_qp(nempty_states, data)
