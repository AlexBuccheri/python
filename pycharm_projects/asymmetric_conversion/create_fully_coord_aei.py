""" In order to apply substitutions to the system, one needs a system where all atoms are connected.
    Unlike in create_aei.cell, I'm:
    a) Deleting all uncoordinated atoms
    b) Finding nearest neighbours of remaining atoms in the central cell.

    One can then perform substitutions on this system
    and fold any positions that lie outside the central cell of the final structure, back inside
    (assuming I can get the new lattice vectors correct)
"""

from modules.fileio.write import xyz

import cell_operations

# Main Routine
directories = cell_operations.Directories(structure='aei',
                                          input='inputs',
                                          output='aei_outputs')

unit_cell, lattice_vectors = cell_operations.get_primitive_unit_cell(directories, visualise=True)
xyz("fully_coordinated/unit_cell", unit_cell)
translations = cell_operations.translations_for_fully_coordinated_unit(unit_cell, lattice_vectors)

# Ah, the below won't work properly if I pass unit_cell_no_dangling to find_neighbour_cell_oxygens
# because the supercell is then constructed from a unit cell that doesn't contain ALL oxygen atoms
#
# Find and delete atoms that are not coordinated
# unit_cell_no_dangling = cell_operations.delete_uncoordinated_atoms(unit_cell)
# xyz("fully_coordinated/unit_no_dangling", unit_cell_no_dangling)

# Find oxygens from other cells that neighbour Si in the central cell
coordinating_atoms = cell_operations.find_neighbour_cell_oxygens(unit_cell, translations)
xyz("fully_coordinated/coordinating_oxy", coordinating_atoms)

