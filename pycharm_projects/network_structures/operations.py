
import numpy as np
import random

from modules.electronic_structure.structure import atoms, bravais


# Central silicon with 4 surrounding oxygen
al = 7.12637
molecule = [atoms.Atom('si', al * np.array([ 0,     0,     0   ])),
            atoms.Atom('o', al * np.array([ 0.25,  0.25,  0.25])),
            atoms.Atom('o', al * np.array([ 0.25, -0.25, -0.25])),
            atoms.Atom('o', al * np.array([-0.25,  0.25, -0.25])),
            atoms.Atom('o', al * np.array([-0.25, -0.25,  0.25]))]

# Test tetrahedron to triangle operation
#triangle = operations.triangle_from(molecule)
#write.xyz("triangle", triangle)



# For a given tetrahedron, remove one oxygen at random
# and move the silicon into the plane of the remaining 3 atoms
#
#for each step, store norm(atom0_si - atom1), norm(atom0_si - atom2) and norm(atom0_si - atom3)
#these should get smaller. Once they start getting larger, stop. atom_si is now in the plane.
# Maybe a more mathematical way of doing it but should be quite fast
#
# Make more robust to deal with fractional or angstrom
#
def triangle_from(tetrahedron):

    # Separate silicons and oxygens
    triangle = []
    silicon = []
    for atom in tetrahedron:
        if atom.species.lower() == 'si':
            silicon.append(atom)
        elif atom.species.lower() == 'o':
            triangle.append(atom)

    assert len(triangle) == 3, "triangle_from function expects tetrahedron containing 3 oxygen"
    assert len(silicon) == 1, "triangle_from function expects tetrahedron containing 1 silicon"

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

    print("bond_length",bond_length)

    # Alt way of stepping
    # sep = 1000
    # pos = np.zeros(shape=(3))
    # for dr in np.linspace(0, bond_length, sep):

    # ang
    step = 0.001 * bond_length
    pos = np.zeros(shape=(3))
    for dr in np.arange(0 + step, bond_length + step, step):
        pos += si_0 + (dr * unit_normal)
        # print("pos,",pos)
        # print(np.linalg.norm(pos - triangle[0].position), \
        #       np.linalg.norm(pos - triangle[1].position), \
        #       np.linalg.norm(pos - triangle[2].position))
        # print(si_o_0, si_o_1, si_o_2)

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



# For a given tetrahedron, remove one oxygen at random
# and move the silicon into the plane of the remaining 3 atoms
#
#for each step, store norm(atom0_si - atom1), norm(atom0_si - atom2) and norm(atom0_si - atom3)
#these should get smaller. Once they start getting larger, stop. atom_si is now in the plane.
# Maybe a more mathematical way of doing it but should be quite fast
#
# Make more robust to deal with fractional or angstrom
#
# def triangle_from(tetrahedron):
#
#     # Separate silicons and oxygens
#     triangle = []
#     silicon = []
#     for atom in tetrahedron:
#         if atom.species.lower() == 'si':
#             silicon.append(atom)
#         elif atom.species.lower() == 'o':
#             triangle.append(atom)
#
#     assert len(triangle) == 4, "triangle_from function expects tetrahedron containing 4 oxygen"
#     assert len(silicon) == 1, "triangle_from function expects tetrahedron containing 1 silicon"
#
#     # Remove one oxygen
#     index = random.randrange(0, 4)
#     triangle.pop(index)
#
#     # Define 2 of 3 vertices that form the oxygen triangle
#     v01 = triangle[1].position - triangle[0].position
#     v12 = triangle[2].position - triangle[1].position
#
#     # Normal to triangle surface
#     normal = np.cross(v01, v12)
#     unit_normal = normal / np.linalg.norm(normal)
#     #print(normal, unit_normal)
#
#     # Silicon - oxygen distances
#     si_o_0 = np.linalg.norm(silicon[0].position - triangle[0].position)
#     si_o_1 = np.linalg.norm(silicon[0].position - triangle[1].position)
#     si_o_2 = np.linalg.norm(silicon[0].position - triangle[2].position)
#
#
#     # Silicon atom is above the plane and one needs to reverse the norm
#     if np.dot(silicon[0].position - triangle[0].position, unit_normal) > 1:
#         unit_normal = -unit_normal
#
#     # Progagate si incrementally along unit normal, until distances from all oxy are minimised
#     bond_length = si_o_0
#     si_0 = np.asarray(silicon[0].position)
#
#     # Alt way of stepping
#     # sep = 1000
#     # pos = np.zeros(shape=(3))
#     # for dr in np.linspace(0, bond_length, sep):
#
#     # ang
#     step = 0.001
#     pos = np.zeros(shape=(3))
#     for dr in np.arange(0 + step, bond_length + step, step):
#         pos += si_0 + (dr * unit_normal)
#
#         if (np.linalg.norm(pos - triangle[0].position) <= si_o_0) and \
#            (np.linalg.norm(pos - triangle[1].position) <= si_o_1) and \
#            (np.linalg.norm(pos - triangle[2].position) <= si_o_2):
#             si_o_0 = np.linalg.norm(pos- triangle[1].position)
#             si_o_1 = np.linalg.norm(pos- triangle[1].position)
#             si_o_2 = np.linalg.norm(pos - triangle[2].position)
#             final_pos = pos
#         else:
#          break
#
#
#     triangle.append(atoms.Atom('B',final_pos))
#     return triangle



#For insertion, just bang in B-O-B with O centred on Si and visualise - might then be able to see what the sensible
#operation on the remaining oxygens are - will be harder as attached to the rest of the network
# def insert_bob(tetrahedron):
#
#     two_triangles = []
#
#     return two_triangles
