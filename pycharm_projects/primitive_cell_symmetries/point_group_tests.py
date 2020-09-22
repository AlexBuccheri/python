import numpy as np
import spglib

from modules.electronic_structure.structure import bravais
from modules.electronic_structure.crystal_point_groups import crystal_point_groups as wykcoff
from modules.electronic_structure.structure import lattice


# Sample real-space lattice

# Sample reciprocal space lattice


# All points mus map onto another point under wyckoff operations
# If so, that point group confirmed

# Use the MP scheme for centreing the grid - Ah, the BZ is always centred on Gamma
# whereas a primitive cell won't be
# => Just centre my girds accordingly


def find_centre_coefficients(coefficients):
    centre_cofficients = []
    for i_coefficients in coefficients:
        centre_index = int(0.5 * len(i_coefficients)) #python indexing starts at 0
        centre_value = i_coefficients[centre_index]
        centre_cofficients.append(centre_value)
    return centre_cofficients


def sample_grid(lattice, grid_sampling):

    # Have to use odd numbers to ensure grid is always centred
    # (or could use the MP scheme)
    assert grid_sampling[0] % 2 != 0
    assert grid_sampling[1] % 2 != 0
    assert grid_sampling[2] % 2 != 0

    x_coefficients = np.linspace(0, 1, num=grid_sampling[0], endpoint=True)
    y_coefficients = np.linspace(0, 1, num=grid_sampling[1], endpoint=True)
    z_coefficients = np.linspace(0, 1, num=grid_sampling[2], endpoint=True)

    coefficients = [x_coefficients, y_coefficients, z_coefficients]

    centre_coefficients = find_centre_coefficients(coefficients)
    centre = np.matmul(lattice, np.array([centre_coefficients[0],
                                          centre_coefficients[1],
                                          centre_coefficients[2]]))
    grid = []
    for k in z_coefficients:
        for j in y_coefficients:
            for i in x_coefficients:
                # Subtracting centre centres grid on (0,0,0)
                point = np.matmul(lattice, np.array([i, j, k])) - centre
                grid.append(point)
    return grid


def check_equivalent_points_are_in_grid(equivalent_points, grid):
    n_points = len(equivalent_points)
    found = np.full((n_points), False)

    for i,equivalent_point in enumerate(equivalent_points):
        for grid_point in grid:
            if np.all(equivalent_point == grid_point):
                found[i] = True
                break

    return np.all(found == True)

# Get this working. Must be faster:
# def check_equivalent_points_are_in_grid(equivalent_points, grid):
#     n_points = len(equivalent_points)
#     found = np.full((n_points), False)
#     print(grid)
#     for point in equivalent_points:
#         print(point)
#         found += (grid == point).all(-1)
#     quit()
#     return np.all(found == True)


# Main Routine

theta = np.pi   #1 * (np.pi / 180)
r_z = np.array( [[np.cos(theta), -np.sin(theta), 0],
                 [np.sin(theta),  np.cos(theta), 0],
                 [0,              0,             1]])

print(np.matmul(r_z, np.array([1,1,1])))
quit()

# Test wykcoff positions

equivalent_points = wykcoff.wyckoff_positions(wykcoff.PointGroup.C_1, np.array([1,2,5.5]))
print(equivalent_points)

quit()

mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])

a = 1
lattice = bravais.face_centred_cubic(a)
grid_sampling = [3, 3, 3]
grid = sample_grid(lattice, grid_sampling)

found = []
for point in grid:

    print(equivalent_points)
    found.append(check_equivalent_points_are_in_grid(equivalent_points, grid))

print(found)
print(np.all(found))

quit()

a = 1
lattice = bravais.simple_cubic(a)
grid_sampling = [3, 3, 3]
grid = sample_grid(lattice, grid_sampling)

# Apply symmetry operation to each point in grid
# The resulting point also be a point in the grid
# If true for all points, the grid has the corresponding point group
# therefore the lattice has that point group
found = []
for point in grid:
    equivalent_points = wykcoff.wyckoff_positions(wykcoff.PointGroup.O_h, point)
    found.append(check_equivalent_points_are_in_grid(equivalent_points, grid))

if np.all(found):
    print("Test passed! Cubic is O_h")
else:
    print("Test failed! Cubic is not O_h")



found = []
for point in grid:
    equivalent_points = wykcoff.wyckoff_positions(wykcoff.PointGroup.C_2v, point)
    found.append(check_equivalent_points_are_in_grid(equivalent_points, grid))

# This passes as all wykcoff positions in C_2v are in O_h.
# That is, all point groups are a subgroup of O_h
if np.all(found):
    print("Test passed! Cubic contains C_2v")
else:
    print("Test failed! Cubic is not C_2v")











