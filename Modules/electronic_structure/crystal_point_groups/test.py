import numpy as np
import enum 

class PointGroup(enum.Enum):
    C_1 = enum.auto
    C_i = enum.auto

def c_1(point):
    x = point[0]
    y = point[0]
    z = point[0]
    return np.array([x, y, z])


def wyckoff_positions(point_group, point):    
    switch = {
        PointGroup.C_1 : c_1
        }
    return switch.get(point_group)(point)

# Simple test
symmetry_equivalent_positions = wyckoff_positions(PointGroup.C_1, np.array([1,1,1]))
print(symmetry_equivalent_positions)
