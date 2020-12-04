
import numpy as np

from modules.electronic_structure.structure import bravais
from modules.electronic_structure.structure import atoms


class System():
    def __init__(self, system_label, unit_cell, lattice, al=None):
        self.system_label = system_label
        self.unit_cell = unit_cell
        self.lattice = lattice
        self.al = al


def simple_cubic():
    system_label = 'cubic_silicon'
    al = 10.33
    species = ['Si'] * 8
    positions = al * np.array([
        [0.25000000, 0.75000000, 0.25000000],
        [0.00000000, 0.00000000, 0.50000000],
        [0.25000000, 0.25000000, 0.75000000],
        [0.00000000, 0.50000000, 0.00000000],
        [0.75000000, 0.75000000, 0.75000000],
        [0.50000000, 0.00000000, 0.00000000],
        [0.75000000, 0.25000000, 0.25000000],
        [0.50000000, 0.50000000, 0.50000000]
    ])
    unit_cell = atoms.Atoms(species, positions)
    lattice = bravais.simple_cubic(al)
    return System(system_label, unit_cell, lattice)


def bcc():
    system_label = 'bcc_potassium'
    al = 8.612
    unit_cell = atoms.Atoms(['K'], [[0, 0, 0]])
    lattice = bravais.body_centred_cubic(al)
    return System(system_label, unit_cell, lattice)


def fcc():
    system_label = 'fcc_silicon'
    al = 10.26
    unit_cell = atoms.Atoms(['Si', 'Si'], [[0, 0, 0], [0.25 * al] * 3])
    lattice = bravais.face_centred_cubic(al)
    return System(system_label, unit_cell, lattice)


def simple_tetragonal():
    system_label = 'simple_tet_tio2'
    a = 8.71
    c = 5.65
    lattice = bravais.simple_tetragonal(a, c)
    fractional_positions = np.array([
       [0.000000, 0.000000, 0.000000],
       [0.500000, 0.500000, 0.500000],
       [0.695526, 0.695526, 0.000000],
       [0.304474, 0.304474, 0.000000],
       [0.195526, 0.804474, 0.500000],
       [0.804474, 0.195526, 0.500000]
    ]).transpose()
    positions = np.matmul(lattice, fractional_positions).transpose()
    unit_cell = atoms.Atoms(['Ti', 'Ti', 'O', 'O', 'O', 'O'], positions.tolist())
    return System(system_label, unit_cell, lattice)


# def body_centred_tetragonal():
#     system_label = 'bc_tet_tio2'
#     a =
#     c =
#     lattice = bravais.body_centred_tetragonal(a, c)
#     fractional_positions = np.array([
#        [0.500000, 0.500000, 0.000000],
#        [0.250000, 0.750000, 0.500000],
#        [0.456413, 0.956413, 0.500000],
#        [0.706413, 0.706413, 0.000000],
#        [0.043587, 0.543587, 0.500000],
#        [0.293587, 0.293587, 0.000000],
#     ]).transpose()
#     positions = np.matmul(lattice, fractional_positions).transpose()
#     unit_cell = atoms.Atoms(['Ti', 'Ti', 'O','O','O','O'], positions.tolist())
#     return System(system_label, unit_cell, lattice)


def hexagonal():
    system_label = 'hex_bn'
    a = 4.75
    c = 14.56
    lattice = bravais.hexagonal(a, c)
    fractional_positions = np.array([
        [0.333333, 0.666667, 0.250000],
        [0.666667, 0.333333, 0.750000],
        [0.333333, 0.666667, 0.750000],
        [0.666667, 0.333333, 0.250000]
    ]).transpose()
    positions = np.matmul(lattice, fractional_positions).transpose()
    unit_cell = atoms.Atoms(['B', 'B', 'N', 'N'], positions.tolist())
    return System(system_label, unit_cell, lattice)


def rhomohedral():
    system_label = 'rhom_mos'
    a = 12.96
    alpha = 0.470
    lattice = bravais.rhombohedral(a, alpha, 'rhom')
    fractional_positions = np.array([
        [0.000007, 0.000007, 0.000007],
        [0.254176, 0.254176, 0.254176],
        [0.412458, 0.412458, 0.412458]
    ]).transpose()
    positions = np.matmul(lattice, fractional_positions).transpose()
    unit_cell = atoms.Atoms(['Mo', 'S', 'S'], positions.tolist())
    return System(system_label, unit_cell, lattice)


def simple_orthorhombic():
    system_label = 'simple_orth_cdau'
    a = 5.83
    b = 9.30
    c = 9.74
    lattice = bravais.simple_orthorhombic(a, b, c)
    fractional_positions = np.array([
       [0.00000000, 0.75000000, 0.69415500],
       [0.00000000, 0.25000000, 0.30584500],
       [0.50000000, 0.75000000, 0.18925900],
       [0.50000000, 0.25000000, 0.81074100]
    ]).transpose()
    positions = np.matmul(lattice, fractional_positions).transpose()
    unit_cell = atoms.Atoms(['Cd', 'Cd', 'Au', 'Au'], positions.tolist())
    return System(system_label, unit_cell, lattice)


def simple_monoclinic():
    system_label = 'simpe_mono_lisn'
    a = 9.79
    b = 6.09
    c = 14.75
    beta = 1.85
    lattice = bravais.simple_monoclinic(a, b, c, beta)
    fractional_positions = np.array([
        [0.268026, 0.500000, 0.331701],
        [0.500000, 0.500000, 0.000000],
        [0.731974, 0.500000, 0.668299],
        [0.241273, 0.000000, 0.661946],
        [0.758727, 0.000000, 0.338054],
        [0.000000, 0.000000, 0.000000]
    ]).transpose()
    positions = np.matmul(lattice, fractional_positions).transpose()
    unit_cell = atoms.Atoms(['Li', 'Li', 'Li', 'Sn', 'Sn', 'Sn'], positions.tolist())
    return System(system_label, unit_cell, lattice)
