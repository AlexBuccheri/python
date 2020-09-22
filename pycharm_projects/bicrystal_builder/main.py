
import os

# My modules
from modules.dftb import hsd_types

# List of xyz data with no header or atom type (as all Cu)
def read_plane(fname: str):
    if not os.path.exists(fname):
        exit("File: "+fname+" does not exist")

    fid = open(fname, "r")
    data = fid.readlines()
    fid.close()

    positions = []
    for line in data:
        xyz = [float(i) for i in line.split()]
        positions.append(xyz)
    return positions

# Assign layer indices, centred on central_layer
# Expects an odd number of layers
def assign_layer_indices(n_layers, central_layer = 6):
    assert (n_layers % 2 != 0)

    indices = [central_layer]
    for i in range(1, int(n_layers - 1 / 2)):
        indices += [central_layer + i, central_layer - i]

    indices.sort()
    return indices

# --------------
# Main Routine
# --------------

root = "/Users/alexanderbuccheri/Desktop/Work 2020/bicrystal/cu_planes/"
n_layers = 5
central_layer = 6


plane_8 = read_plane(root + "zplane8.dat")
plane_10 = read_plane(root + "zplane10.dat")

for i in range(len(plane_8)):
    print(plane_8[i][0] - plane_10[i][0], plane_8[i][1] - plane_10[i][1], plane_8[i][2] - plane_10[i][2])


# layer_indices = assign_layer_indices(n_layers)
# bicrystal_positions = []
# for layer in layer_indices:
#     layer_positions = read_plane(root + "zplane" + str(layer) + ".dat")
#     bicrystal_positions += layer_positions
#
# n_atoms = len(bicrystal_positions)

# Write as xyz



# Need lattice vectors for (x,y)

# Convert to DFTB+ format - see what prior routines I have for that!!! OR use ASE
#hsd_types.Geometry()
