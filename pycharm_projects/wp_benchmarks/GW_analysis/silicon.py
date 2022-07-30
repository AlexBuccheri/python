""" Post-process silicon G0W0 results.
"""
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
from typing import Tuple

from src.conversions import ha_to_ev
from src.parsers import band_gaps


def band_gaps_vs_q(root: str, q_points: list) -> dict:
    """ Get fundamental and direct Gamma gaps as a function of q.

    :param path: Path to file.
    :param q_points: List of q-points.
    :return: Fundamental and Gamma-Gamma gaps.
    """
    fundamental = np.empty(shape=len(q_points))
    gamma_gamma = np.empty(shape=len(q_points))

    for iq, q_point in enumerate(q_points):
        q_str = "".join(str(x) for x in q_point)
        path = os.path.join(root, q_str)
        gaps: dict = band_gaps(path)
        fundamental[iq] = gaps['fundamental']
        gamma_gamma[iq] = gaps['gamma']

    return {'fundamental': fundamental, 'gamma_gamma': gamma_gamma}


def initialise_indirect_gap_plot() -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    ax.set_title('Convergence in Indirect Gap of Silicon')
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})
    ax.set_xlabel('q-points', fontsize=16)
    ax.set_ylabel('Quasiparticle Band Gap (eV)', fontsize=16)
    ax.set_ylim(1.0, 1.6)
    return fig, ax


def initialise_gamma_gap_plot() -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    ax.set_title('Convergence in Gamma-Gamma Gap of Silicon')
    fig.set_size_inches(14, 10)
    plt.rcParams.update({'font.size': 16})
    ax.set_xlabel('q-points', fontsize=16)
    ax.set_ylabel('Quasiparticle Band Gap (eV)', fontsize=16)
    ax.set_ylim(3.0, 3.6)
    return fig, ax


if __name__ == "__main__":
    # Notebook Results directory
    root = 'GW_results/results/silicon'
    subdirectory = '_percent_empty'

    q_points = [[2, 2, 2], [4, 4, 4], [6, 6, 6], [8, 8, 8]]
    total_q_points = [np.prod(q_point) for q_point in q_points]
    n_empties = [25, 50, 100]

    fig, ax = initialise_indirect_gap_plot()
    for n_empty in n_empties:
        gaps = band_gaps_vs_q(os.path.join(root, str(n_empty) + subdirectory), q_points)
        ax.plot(total_q_points, gaps['fundamental'] * ha_to_ev, 'o', markersize=10, label=str(n_empty))

    ax.legend(title='Number of Emtpy States (%)')
    plt.show()

    # Must be a smarter way of doing this
    fig, ax = initialise_gamma_gap_plot()
    for n_empty in n_empties:
        gaps = band_gaps_vs_q(os.path.join(root, str(n_empty) + subdirectory), q_points)
        ax.plot(total_q_points, gaps['gamma_gamma'] * ha_to_ev, 'o', markersize=10, label=str(n_empty))

    ax.legend(title='Number of Emtpy States (%)')
    plt.show()

    #     if save_plots:
    #         plt.savefig(file_name, dpi=300, facecolor='w', edgecolor='w',
    #                     orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)

