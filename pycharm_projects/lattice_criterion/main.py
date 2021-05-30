"""
For each of the Bravais lattices, generate a set of maximum
translation integers corresponding to a translation summation cut-off, output the
sueprcells resulting from using these integers in conjunction with a
radial lattice sum and confirm the outputs are as close to spherical as
possible

Tested on:
cubic, BCC, FCC, simple_tetragonal, hexagonal, rhomohedral,
simple_orthorhombic, simple_monoclinic
"""

from modules.electronic_structure.structure import lattice as lattice_module
from modules.electronic_structure.structure import supercell
from modules.fileio import write

from reference_systems import *


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
