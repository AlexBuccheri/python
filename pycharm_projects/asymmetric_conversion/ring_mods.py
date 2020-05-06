# -----------------------
# Read in ring and expand
# -----------------------

import numpy as np

# My modules
from modules.fileio import read, write
from modules.electronic_structure.structure import atoms

species, positions = read.xyz("inputs/ring.xyz")

# https://stackoverflow.com/questions/18305712/how-to-compute-the-center-of-a-polygon-in-2d-and-3d-space
# Find ring's centre
# points = [ [x1,y1,z1], [x2,y2,z2], ...]
def find_centre(points):
    s = np.zeros(shape=(3))
    sL = 0
    for i in range(len(points)):
        r0 = np.asarray(points[i - 1])
        r1 = np.asarray(points[i])
        r01 = np.linalg.norm(r1-r0)
        s += 0.5 * r01 * (r0 + r1)
        sL += r01
    return s/sL

centre = find_centre(positions)

# Write something to indexc each tetrahedron, and each corresponding oxygen position per tetrahedron
# Probably want to treat each si with the 2 non-ring oxy as a single unit and move its central position

# Move each atom out along vector with ring centre
l = 5
new_positions = []
for ia in range(0, len(species)):
    position = positions[ia]
    # Point outwards from the centre
    radial_vector = np.asarray(position) - centre
    unit_radial_vector = radial_vector / np.linalg.norm(radial_vector)
    scaled_radial_vector = l * unit_radial_vector
    assert( np.isclose(np.linalg.norm(scaled_radial_vector), l) )
    new_positions.append(scaled_radial_vector + position)

molecule = atoms.Atoms(species, new_positions)
write.xyz("test.xyz", molecule)