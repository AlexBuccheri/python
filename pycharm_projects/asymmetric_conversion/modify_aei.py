""" Convert an AEI zeolite structure to a boron oxide structure
"""

# Modules for this project
import ring_mods
import cell_operations
from modules.fileio.write import xyz


#TODO(Alex) Check that all uncoordinated oxy occupy positions in 'coordinating_atoms'
def full_coordinated_structure(directories):
    """
    :brief Construct a fully-coordinated AEI structure
    Some atoms won't exist in the central unit cell

    :return
    """

    # Unit cell from input cif
    unit_cell, lattice_vectors = cell_operations.get_primitive_unit_cell(directories)
    translations = cell_operations.translations_for_fully_coordinated_unit(unit_cell, lattice_vectors)

    # Find oxygens from other cells that neighbour those in the central cell
    coordinating_atoms = cell_operations.find_neighbour_cell_oxygens(unit_cell, translations)

    # Find and delete atoms that are not coordinated - assume these MUST occupy sites in adjacent cells
    # captured by 'coordinating_atoms', hence ok to do do
    unit_cell_no_dangling = cell_operations.delete_uncoordinated_atoms(unit_cell)

    return coordinating_atoms + unit_cell_no_dangling



def get_structural_components(unit_cell, visualise=False):
    """
    :brief For an input unit cell, split into main structural components
        upper ring, lower ring and connecting chains

     ghost atoms are duplicates present in the two connecting
     structures (i.e. ring + chain), and act to indicate the bonding direction

     NOTE: those returned are erroenous. I've tabulated Si but should
     have looked at oxy

    :return structures :  List of structures
    """

    upper_ring = ring_mods.get_ring(unit_cell, lambda atom: atom.position[0] > 11.5,
                                  erroneous_indices=[4, 11])
    lower_ring = ring_mods.get_ring(unit_cell, lambda atom: atom.position[0] <= 7,
                                    erroneous_indices=[3, 12])

    # Don't need these ghosts: the terminating oxygens are the actual ghost/duplicate atoms
    chain1, chain2 = ring_mods.get_connecting_chains(unit_cell)

    structures = {'upper_ring': upper_ring, "lower_ring": lower_ring,
                  'chain1': chain1, 'chain2': chain2}

    if visualise:
        xyz('fully_coordinated/upper_ring', structures['upper_ring'])
        xyz('fully_coordinated/lower_ring', structures['lower_ring'])
        xyz('fully_coordinated/chain1', structures['chain1']['main'])
        xyz('fully_coordinated/chain2', structures['chain2']['main'])
        xyz('fully_coordinated/ghost1', structures['chain1']['ghost'])
        xyz('fully_coordinated/ghost2', structures['chain2']['ghost'])

    return structures


def main():
    """
    :brief Convert an AEI zeolite structure to a borate structure
           Main routine

    :return borate: borate structure
    """
    directories = cell_operations.Directories(structure='aei',
                                              input='inputs',
                                              output='aei_outputs')
    silicate = full_coordinated_structure(directories)
    xyz('fully_coordinated/fully_coord', silicate)

    print("If atomic ordering changes => Different starting structure,"
          "I break get_structural_components")
    structures = get_structural_components(silicate, visualise=False)

    # Boron framework components
    upper_ring = ring_mods.ring_substitutions(structures['upper_ring'])
    xyz('fully_coordinated/upper_ring_swap', upper_ring)

    # lower_ring = ring_substitutions(structures['lower_ring'])
    # chains = chain_substitutions(structures['chains'])


    # Piece the borate components together

    # Some guess at the new lattice vectors

    # Fold atoms back into central cell and remove ~ duplicates

    # Return the final structure for relaxing
    return

main()