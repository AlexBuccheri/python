import numpy as np
from scipy import spatial
from typing import List

from modules.electronic_structure.structure import supercell
from modules.electronic_structure.structure import atoms
from modules.fileio.write import xyz

#TODO(Alex) Should remove default for upper_bound_length
def find_atoms_neighbouring_central_cell(unit_cell: atoms.Atoms, translations: list,
                                         upper_bound_length=1.8, visualise=False):
    """
    :brief: For a given unit cell, list all atoms in adjacent cells that are
    neighbours with atoms in the (central) unit cell

    :param unit_cell: A list of atoms in the unit cell
    :param translations: List of translation vectors
    :param upper_bound_length: Upper bound for NN bond length of AEI framework (in angstrom)
    :param visualise: Optional, output intermediate structures
    :return coordinating_atoms: A list of atoms in adjacent unit cells that neighbour
    atoms in the central unit_cell.
    """
    assert len(unit_cell) > 0, "len(unit_cell) = 0"
    assert isinstance(unit_cell, list), "unit_cell should be list[atoms.Atom]"
    assert isinstance(unit_cell[0], atoms.Atom), "unit_cell should be list[atoms.Atom]"
    n_atoms_prim = len(unit_cell)

    # Move central cell translation vector to start of list
    # Makes central-cell atoms entries [0: n_atoms_prim]
    zero_translation = np.array([0, 0, 0])
    for i, translation in enumerate(translations):
        if np.array_equal(translation, zero_translation):
            translations.insert(0, translations.pop(i))
            break

    # Super cell = Unit cell plus all coordinating cells
    super_cell = supercell.build_supercell(unit_cell, translations)
    assert len(super_cell) == len(translations) * n_atoms_prim

    positions = [atom.position for atom in super_cell]
    d_supercell = spatial.distance_matrix(positions, positions)
    assert d_supercell.shape[0] == len(super_cell), "d_supercell.shape[0] != len(super_cell)"
    assert d_supercell.shape[1] == d_supercell.shape[0], "distance matrix is not square"

    # For atoms (ia) in central cell, check for neighbours inside and outside of central cell
    coordinating_atom_indices = []
    for ia in range(0, n_atoms_prim):
        # This should only check for neighbours outside of the central cell
        # but doesn't give the correct answer and I can't see why
        # indices = np.where((d_supercell[ia, n_atoms_prim:] > 0.) &
        #                    (d_supercell[ia, n_atoms_prim:] <= upper_bound_length))[0]
        indices = np.where((d_supercell[ia, :] > 0.) &
                           (d_supercell[ia, :] <= upper_bound_length))[0]
        coordinating_atom_indices.extend(indices)

    # Remove dups
    coordinating_atom_indices = list(set(coordinating_atom_indices))
    # Remove indices associated with atoms in central cell
    coordinating_atom_indices = np.asarray(coordinating_atom_indices)
    coordinating_atom_indices = coordinating_atom_indices[coordinating_atom_indices >= n_atoms_prim]

    # Store coordinating atoms
    coordinating_atoms = [super_cell[ia] for ia in coordinating_atom_indices]

    if visualise:
        #TODO(Alex) Resolve this to make general
        output_dir = 'aei_outputs'
        xyz(output_dir + '/' + "aei_supercell", super_cell)
        xyz(output_dir + '/' + "aei_central_cell", unit_cell)
        xyz(output_dir + '/' + "aei_neighbours", coordinating_atoms)

    return coordinating_atoms


def index_loose_atoms(unit_cell_positions: list, upper_bound_length: float):
    """
    Index loose or uncoordinated atoms

    :param unit_cell: Positions of atoms in the unit cell
    :param upper_bound_length: Upper bound for NN bond length (in angstrom)
    :return loose_atom_indices: Indices of atoms without any neighbours in central cell
    """
    d = spatial.distance_matrix(unit_cell_positions, unit_cell_positions)

    loose_atom_indices = []
    for ia in range(0, len(unit_cell_positions)):
        indices = np.where((d[ia, :] > 0.) & (d[ia, :] <= upper_bound_length))[0]
        no_neighbours = len(indices) == 0
        if no_neighbours:
            loose_atom_indices.append(ia)

    return loose_atom_indices


class BondLengthBounds:
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper


def find_equivalent_position(loose_atom_index: int, unit_cell: list, translations: list,
                             bond_bounds: BondLengthBounds) -> list:
    """
    Find the position of the loose_atom_index in an adjacent cell, for which it's
    a neighbour of an atom in the central cell (T=0)

    Target bond length is ~ 1.529 - 1.6 angstrom for AEI

    :param loose_atom_index: Index of uncoordinated atom
    :param unit_cell:
    :param translations:
    :param bond_bounds Lower and upper limits for bond length
    :return: equivalent position for an atom indexed by loose_atom_index, which
    gives it a coordination
    """
    lower_bound = bond_bounds.lower
    upper_bound = bond_bounds.upper

    unit_cell_positions = [atom.position for atom in unit_cell]
    loose_atom_position = unit_cell[loose_atom_index].position
    loose_atom_positions = [loose_atom_position + translation for translation in translations]
    d_loose = spatial.distance_matrix(loose_atom_positions, unit_cell_positions)

    valid_equivalent = []
    for ia in range(0, len(loose_atom_positions)):
        # neighbours of loose atom in central cell
        indices = np.where((d_loose[ia, :] > lower_bound) & (d_loose[ia, :] <= upper_bound))[0]
        if len(indices) > 0:
            # print("distance:", d_loose[ia, indices])
            valid_equivalent.append(loose_atom_positions[ia].tolist())

    assert len(valid_equivalent) == 1, "Should only be one equivalent position in adjacent " \
                                       "cells that is a nearest neighbour to an atom in the central" \
                                       "cell"
    return valid_equivalent[0]


def replace_loose_atoms(unit_cell: list, loose_atom_indices: List[int], equivalent_positions: list,
                        visualise_parts=False) -> list:
    """
    Remove uncoordinated atoms and their equivalents to the unit cell

    :param unit_cell: List of atoms
    :param loose_atom_indices: Indices for uncoordinated atoms
    :param equivalent_positions: List of equivalent positions for uncoordinated atoms.
    of size len(loose_atom_indices)
    :param visualise_parts: bool, output .xyz of original unit_cell, removed atoms and
     replacement atoms
    :return new_unit_cell with no uncoordinated atoms
    """
    assert len(loose_atom_indices) == len(equivalent_positions), "Should be one new (equivalent) " \
                                                                 "position per uncoordinated atom"

    # Can't use pop as it changes the indexing each time
    new_unit_cell = []
    removed_atoms = []
    for ia in range(0, len(unit_cell)):
        if ia not in loose_atom_indices:
            new_unit_cell.append(unit_cell[ia])
        else:
            removed_atoms.append(unit_cell[ia])

    #TODO(Alex) Generalise
    # I know they're all oxygens but should treat this correctly if making general
    replacements = [atoms.Atom(position=position, species='O') for position in equivalent_positions]

    if visualise_parts:
        output_dir = 'aei_outputs'
        xyz(output_dir + '/' + "aei_central_cell", unit_cell)
        xyz(output_dir + '/' + "aei_replacements_cell", replacements)
        xyz(output_dir + '/' + "aei_removed_atoms", removed_atoms)

    return new_unit_cell + replacements


def replace_uncoordinated_atoms(unit_cell: list, translations: list) -> List[atoms.Atom]:
    """
    Identify non-bonded atoms in the unit cell and replace with their equivalent positions
    in adjacent unit cells that result in a fully connected structure

    Translate each loose atom by all translation vectors and see if any
    equivalent position puts each loose atom as a NN of an atom in the central cell
    Retain these and dump the uncoordinated atoms, producing a fully-coordinated structure.

    This is the structure one works with for boron substitutions. Once the boron framework
    has been constructed, one wraps all atoms outside the cell to inside.

    :param unit_cell: atoms.Atoms Unit cell of atoms
    :param translations: List of translation vectors
    :return unit_cell: Unit cell with no uncoordinated atoms
    """

    # AEI Si and O bond bounds (angstrom)
    bond_bounds = BondLengthBounds(1.4, 1.8)
    n_atoms = len(unit_cell)

    unit_cell_positions = [atom.position for atom in unit_cell]
    loose_atom_indices = index_loose_atoms(unit_cell_positions, bond_bounds.upper)
    assert len(loose_atom_indices) == 4, \
        "For the AEI structure, expect 4 Oxy atoms with no connections in the central cell"

    # One equivalent position per loose atom
    equivalent_positions = [find_equivalent_position(ia, unit_cell, translations, bond_bounds)
                            for ia in loose_atom_indices]

    unit_cell = replace_loose_atoms(unit_cell, loose_atom_indices, equivalent_positions)
    assert len(unit_cell) == n_atoms

    output = True
    if output:
        output_dir = 'aei_outputs'
        xyz(output_dir + '/' + "aei_primitive_cell_fullcon", unit_cell)
        print("Output a cell with full connectivity. "
              "Atoms that lie outside the primitive cell can be folded back in")

    return unit_cell
