
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


for ia,atom in enumerate(beta_crist):

    if atom.species.lower() == 'o':
        silicons = neighbours[ia]

        # Oxygen corner-shares with two tetrahedra
        if len(silicons) == 2:
            tetrahedra = []
            for silicon in silicons:
                tetrahedron = []
                for oxy in neighbours[silicon]:
                    tetrahedron.append(beta_crist[oxy])
                tetrahedra.append(tetrahedron)
            assert(len(tetrahedra) == 2)




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