import numpy as np
from scipy import spatial

from modules.electronic_structure.structure import supercell
from modules.electronic_structure.structure import atoms
from modules.fileio.write import xyz


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
    for i,translation in enumerate(translations):
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
        output_dir = 'aei_outputs'
        xyz(output_dir + '/' + "aei_supercell", super_cell)
        xyz(output_dir + '/' + "aei_central_cell", unit_cell)
        xyz(output_dir + '/' + "aei_neighbours", coordinating_atoms)

    return coordinating_atoms
