""" For each of the Bravais lattices, generate a set of maximum
translation integers corresponding to a translation summation cut-off, output the
sueprcells resulting from using these integers in conjunction with a
radial lattice sum and confirm the outputs are as close to spherical as
possible"""

import numpy as np

from modules.electronic_structure.structure import bravais
from modules.electronic_structure.structure import lattice as lattice_module
from modules.electronic_structure.structure import supercell
from modules.electronic_structure.structure import atoms
from modules.fileio import write


class System():
    def __init__(self, system_label, unit_cell, lattice):
        self.system_label = system_label
        self.unit_cell = unit_cell
        self.lattice = lattice


# Use proper systems as examples: Finding nearest neighbours will make sense
# Lattice constants are not as important
# Cover most systems. If these all look fine, we're probably ok


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


def lattice_sum(lattice, n, cutoff):
    """
    Sum over translation vectors, using a radial criterion
    :param lattice:
    :param n:
    :param cutoff:
    :return:
    """
    cutoff_squared = cutoff * cutoff
    translations = []
    for k in range(-n[2], n[2]):
        for j in range(-n[1], n[1]):
            for i in range(-n[0], n[0]):
                r = np.matmul(lattice, np.array([i, j, k]))
                if np.dot(r, r) <= cutoff_squared:
                    translations.append(r)
    return translations


def output_cell_from_cubic_criterion(system, cutoff):
    """
    Works for simple cubic cells, but
    is not appropriate as cell angles deviation from 90 degrees
    This typically underestimates max_integers by 1 (see BCC or FCC), such
    that one can improve the situation by doing:
    max_integers = [i + 1 for i in max_integers]
    """

    name = system.system_label
    lattice = system.lattice
    unit_cell = system.unit_cell
    max_integers = lattice_module.simple_cubic_cell_translation_integers(lattice, cutoff)
    print("Cubic criterion for " + name, "Lattice integers:", max_integers)

    translations = lattice_sum(lattice, max_integers, cutoff)
    super_cell = supercell.build_supercell(unit_cell, translations)
    write.xyz(name + "_cubic_criterion" , super_cell)
    return


def output_cell_from_general_criterion(system, cutoff):
    """
    (Hopefully) works for any cell
    Expect it to give the same max integers as cubic cells
    """
    name = system.system_label
    lattice = system.lattice
    unit_cell = system.unit_cell
    max_integers = lattice_module.translation_integers_for_radial_cutoff(lattice, cutoff)
    print("General criterion for " + name, "Lattice integers:", max_integers)

    translations = lattice_sum(lattice, max_integers, cutoff)
    super_cell = supercell.build_supercell(unit_cell, translations)
    write.xyz(name + "_general_criterion" , super_cell)
    return

# Different approve that kind of works but one ultimately doesn't need to use
# def test(a, cutoff):
#     b = np.array([cutoff, cutoff, cutoff])
#     return 2 * np.floor(np.linalg.solve(a, b))
#
# print("test:", test(lattice, cutoff))


cutoff = 40

for system in [simple_cubic(), bcc(), fcc(), simple_tetragonal(), hexagonal(), rhomohedral(),
               simple_orthorhombic(), simple_monoclinic()]:
    output_cell_from_cubic_criterion(system, cutoff)
    output_cell_from_general_criterion(system, cutoff)
