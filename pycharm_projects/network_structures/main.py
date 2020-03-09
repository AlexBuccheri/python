
import numpy as np

import ase.io
from ase.atoms import Atoms
import spglib

from modules.electronic_structure.structure import atoms, bravais
from modules.fileio import write

import utils
import operations



# Beta Cristobolite structure from CIF to ASE
ase_data= ase.io.read("beta_crist_EntryWithCollCode77458.cif", store_tags=False)
# ase.io.write("cell.xyz", ase_data)

# Initial silicate system
beta_crist = utils.ase_atom_my_atom(ase_data)


# NN bond length in fractional
nn_radius = np.sqrt(3.) / 4.
# For some reason with this system, NN distance is actually 0.np.sqrt(3.) / 8. in fractional
nn_radius = np.sqrt(3.) / 8.

# Alt approach to this would be to construct a distance matrix, then apply the cutoff
neighbours = utils.neighbour_list(beta_crist, nn_radius)

print("N atoms: ",len(beta_crist))
print(" Should visualise this system - it's not what I thought from vesta view. Smaller")

corner_sharing_oxy = utils.find_corner_sharing_oxy(beta_crist, neighbours)

print("corner_sharing_oxy:", corner_sharing_oxy)


class SnipClass:
    def __init__(self, old_atom_index, new_atom_index = None, old_atom = None, new_atom = None, modifier = 'replace'):
        self.old_atom_index = old_atom_index
        if new_atom_index is not None: self.new_atom_index = new_atom_index
        if old_atom is not None: self.old_atom = old_atom
        if new_atom is not None: self.new_atom = new_atom
        self.modifier = modifier
        if modifier not in ['replace','add','remove']:
            exit("Invalid choice of modifier in SnipClass")


# Pass index of one corner-sharing oxygen
def local_snip(molecule, neighbours, shared_oxy):

    silicons = neighbours[shared_oxy]
    assert(len(silicons) == 2)

    # List corresponding tetrahedrons wth common corner-sharing
    # oxygen removed
    tetrahedra = []
    for silicon in silicons:
        tetrahedron = [molecule[silicon]]
        for oxy in neighbours[silicon]:
            if oxy != shared_oxy: tetrahedron.append(molecule[oxy])
        tetrahedra.append(tetrahedron)

    molecule = []
    for tetrahedron in tetrahedra:
        for atom in tetrahedron:
            molecule.append(atom)
    write.xyz("two_tet.xyz", molecule, "Two tetra")

    # tetrahedron -> triangle
    triangles = []
    for tetrahedron in tetrahedra:
        # If the tetrahedron has one si and 3 oxy
        if len(tetrahedron) == 4:
            triangles.append(operations.triangle_from(tetrahedron))
        else:
            triangles.append([])


    # Remove shared oxygen last, due to indexing
    changes = [SnipClass(old_atom_index=shared_oxy, modifier='remove')]
    # Extract B positions and labels to B
    for i,silicon in enumerate(silicons):
        new_atom = [atom for atom in molecule if atom.species.lower()=='si'][0]
        change = SnipClass(old_atom_index=silicon, new_atom=new_atom, modifier='replace')
        changes.append(change)

    return operations

changes = local_snip(beta_crist, neighbours, corner_sharing_oxy[6])


# Apply replacements to structure
for change in changes:
    if change.modifier == 'replace':
        beta_crist[change.old_atom_index] = change.new_atom

#Apply changes that affect number of atoms, hence atomic indexing
for change in changes:
    if change.modifier == 'remove':
        beta_crist.pop(change.old_atom_index)



# Apply to all tetrahedra in the structure

# Write another test to apply at random







# for ia,atom in enumerate(beta_crist):
#
#     if atom.species.lower() == 'o':
#         silicons = neighbours[ia]
#
#         # Oxygen corner-shares with two tetrahedra
#         if len(silicons) == 2:
#             tetrahedra = []
#             for silicon in silicons:
#                 tetrahedron = []
#                 for oxy in neighbours[silicon]:
#                     tetrahedron.append(beta_crist[oxy])
#                 tetrahedra.append(tetrahedron)
#             assert(len(tetrahedra) == 2)




# Go through list. If a given atom has two Si NN => bridging oxygen
# Save the two corresponding triangles and disregard the connecting oxy
# Apply triangle routine above (needs a few mods)
# Add atoms back into structure (most likely only the moved silicons)

# Once shown to work, do at random































quit("Before reading in ASE data")

# # Convert from fractional to ang, with guess for lattice constant
# inv_lattice = bravais.body_centred_cubic(5.)
# triangle = []
# for atom in triangle_fractional:
#     fractional = atom.position
#     xyz = np.matmul(inv_lattice, fractional)
#     triangle.append(atoms.Atom(atom.species, xyz))



# Extend into super cell:  super_cell = ase_data.repeat(2)


# Find nearest neighbours of each atom
def neighbour_list(ase_data, neighbour_radius):
    n_atoms = len(ase_data)
    neighbours = []

    for ia in range(0, n_atoms):
        #atom_i = ase_data[ia]
        pos_i = ase_data[ia].position
        neighbours_of_atom_i = []
        for ja in range(0, n_atoms):
            #atom_j = np.asarray(ase_data[ja])
            pos_j = ase_data[ja].position
            separation = np.linalg.norm(pos_i - pos_j)
            if (separation <= neighbour_radius) and (ia != ja):
                neighbours_of_atom_i.append(ja)
        neighbours.append(neighbours_of_atom_i)

    return neighbours


NN_radius = 2.55  #Assume ang
NN_neighbours = neighbour_list(ase_data, NN_radius)

print("Does not look correct. Should not have more than 4 NN")
for ia, neighbours in enumerate(NN_neighbours):
    print(ia, neighbours)

# try modifying the cell to demonstrate that my functions for insertion and removing work

# Try minimising that cell with periodic B/C

# Try applying operations on a framework









# ASE to SPG Don't require here
#beta_crist = utils.ase_atom_to_spg_atom(ase_data)

# Apply operation throughout and at each, run the recommended MM simulation. If this is ok, add the 2nd operation