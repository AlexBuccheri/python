import numpy as np
import ase.io
from ase.atoms import Atoms
import numpy as np
import spglib

from modules.electronic_structure.structure import atoms, bravais
from modules.fileio import write

import utils

# Central silicon with 4 surrounding oxygen
al = 7.12637
molecule = [atoms.Atom('si', al * np.array([ 0,     0,     0   ])),
            atoms.Atom('o', al * np.array([ 0.25,  0.25,  0.25])),
            atoms.Atom('o', al * np.array([ 0.25, -0.25, -0.25])),
            atoms.Atom('o', al * np.array([-0.25,  0.25, -0.25])),
            atoms.Atom('o', al * np.array([-0.25, -0.25,  0.25]))]



# For a given tetrahedron, remove one atom (at random once working)
# and move the silicon into the plane of the remaining 3 atoms
#
def bob(tetrahedron):

    # Separate silicons and oxygens
    triangle = []
    silicon = []
    for atom in tetrahedron:
        if atom.species == 'si':
            silicon.append(atom)
        elif atom.species == 'o':
            triangle.append(atom)

    assert(len(triangle) == 3)
    assert(len(silicon) == 1)

    # Remove one oxygen, although already done in test case
    # Add me

    # Define 2 of 3 vertices that form the oxygen triangle
    v01 = triangle[1].position - triangle[0].position
    v12 = triangle[2].position - triangle[1].position

    # Normal to triangle surface
    normal = np.cross(v01, v12)
    unit_normal = normal / np.linalg.norm(normal)
    #print(normal, unit_normal)

    # Silicon - oxygen distances
    si_o_0 = np.linalg.norm(silicon[0].position - triangle[0].position)
    si_o_1 = np.linalg.norm(silicon[0].position - triangle[1].position)
    si_o_2 = np.linalg.norm(silicon[0].position - triangle[2].position)


    # Silicon atom is above the plane and one needs to reverse the norm
    if np.dot(silicon[0].position - triangle[0].position, unit_normal) > 1:
        unit_normal = -unit_normal

    # Progagate si incrementally along unit normal, until distances from all oxy are minimised
    bond_length = si_o_0
    si_0 = np.asarray(silicon[0].position)

    # Alt way of stepping
    # sep = 1000
    # pos = np.zeros(shape=(3))
    # for dr in np.linspace(0, bond_length, sep):

    # ang
    step = 0.001
    pos = np.zeros(shape=(3))
    for dr in np.arange(0 + step, bond_length + step, step):
        pos += si_0 + (dr * unit_normal)
        print(np.linalg.norm(pos - triangle[0].position),
              np.linalg.norm(pos - triangle[1].position),
              np.linalg.norm(pos - triangle[2].position))

        if (np.linalg.norm(pos - triangle[0].position) <= si_o_0) and \
           (np.linalg.norm(pos - triangle[1].position) <= si_o_1) and \
           (np.linalg.norm(pos - triangle[2].position) <= si_o_2):
            si_o_0 = np.linalg.norm(pos- triangle[1].position)
            si_o_1 = np.linalg.norm(pos- triangle[1].position)
            si_o_2 = np.linalg.norm(pos - triangle[2].position)
            final_pos = pos
        else:
         break


    triangle.append(atoms.Atom('B',final_pos))
    return triangle



# Tetrahedron with last oxygen removed
test_molecule = [atoms.Atom('si',al * np.array([ 0,     0,     0   ])),
                 atoms.Atom('o', al * np.array([ 0.25,  0.25,  0.25])),
                 atoms.Atom('o', al * np.array([ 0.25, -0.25, -0.25])),
                 atoms.Atom('o', al * np.array([-0.25,  0.25, -0.25]))]

triangle = bob(test_molecule)
write.xyz("triangle", triangle)



# triangle2 = []
# for atom in triangle:
#     xyz = al * atom.position
#     triangle2.append(atoms.Atom(atom.species, xyz))
# write.xyz("triangle", triangle2)

#atom0_si should move along 0.25*unit_n1
#for each step, store norm(atom0_si - atom1), norm(atom0_si - atom2) and norm(atom0_si - atom3)
#these should get smaller. Once they start getting larger, stop. atom_si is now in the plane.
# Maybe a more mathematical way of doing it but should be quite fast


#For insertion, just bang in B-O-B with O centred on Si and visualise - might then be able to see what the sensible
#operation on the remaining oxygens are - will be harder as attached to the rest of the network







quit("Before reading in ASE data")

# # Convert from fractional to ang, with guess for lattice constant
# inv_lattice = bravais.body_centred_cubic(5.)
# triangle = []
# for atom in triangle_fractional:
#     fractional = atom.position
#     xyz = np.matmul(inv_lattice, fractional)
#     triangle.append(atoms.Atom(atom.species, xyz))


# Beta Cristobolite structure from CIF to ASE
ase_data= ase.io.read("beta_crist_EntryWithCollCode77458.cif", store_tags=False)
ase.io.write("cell.xyz", ase_data)

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