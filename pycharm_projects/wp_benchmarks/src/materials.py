"""
WP2 Benchmarks, written in Python data structures
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class MoS2WS2Bilayer:
    """ Tabulated E1 System.

    Ref: https://github.com/nomad-coe/greenX-wp2/blob/main/Benchmarks/Heterobilayers/MoS2-WS2.xyz
    :return:
    """
    # lattice = np.array([[3.168394160510246, -3.3524146226426544e-10,  0.0],
    #                     [-1.5841970805453853, 2.7439098312114987, 0.0],
    #                     [6.365167633892244e-18 -3.748609667104484e-30, 39.58711265]])

    # Lattice with ~ zero terms rounded
    lattice = np.array([[3.168394160510246, 0.0, 0.0],
                        [-1.5841970805453853, 2.7439098312114987, 0.0],
                        [0.0, 0.0, 39.58711265]])

    # Positions in angstrom
    positions = np.array([[0.00000000, 0.00000000, 16.68421565],
                          [1.58419708, 0.91463661, 18.25982194],
                          [1.58419708, 0.91463661, 15.10652203],
                          [1.58419708, 0.91463661, 22.90251866],
                          [0.00000000, 0.00000000, 24.46831689],
                          [0.00000000, 0.00000000, 21.33906353]])
    elements = ['W', 'S', 'S', 'Mo', 'S', 'S']
    atomic_numbers = [74, 16, 16, 42, 16, 16]
