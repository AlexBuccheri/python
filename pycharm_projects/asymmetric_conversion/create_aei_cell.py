""" In order to apply substitutions to the system, one needs a system where all atoms are connected.
    As such, I'm translating all loose atoms into positions in other unit cells BUT are connected as
    nearest neighbours to atoms in the central cell. One can then perform substitutions on this system
    and fold any positions that lie outside the central cell of the final structure, back inside
"""

# External libraries
import cell_operations

# Main Routine
directories = cell_operations.Directories(structure='aei',
                                          input='inputs',
                                          output='aei_outputs')

unit_cell, lattice_vectors = cell_operations.get_primitive_unit_cell(directories, visualise=True)
translations = cell_operations.translations_for_fully_coordinated_unit(unit_cell, lattice_vectors)

# Find atoms neighbouring central cell
coordinating_atoms = cell_operations.find_atoms_neighbouring_central_cell(unit_cell, translations, visualise=True)

# Replace uncoordinated atoms in central cell (doesn't utilise coordinating_atoms)
unit_cell_replaced = cell_operations.replace_uncoordinated_atoms(unit_cell, translations)

# Check moved atoms can be restored to positions in the central cell
unit_cell_restored = cell_operations.ensure_atoms_in_central_cell(
    unit_cell_replaced, lattice_vectors)


# -------------------------------------------------------
# Other ideas/notes
# -------------------------------------------------------
# Modify rings by doing one B-O-B
# Structurally relax the lattice maintaining bond angles
# OR
# Replace around the entirety of each ring before relaxing

# Once one has the supercell, could also randomly apply the snip (whilst retaining global structural connectivity)