import numpy as np 
import enum 

def wyckoff_positions_c_1(x, y, z): 
    positions = np.array([[x, y, z]]) 
    return positions 

def wyckoff_positions_c_i(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, -z]]) 
    return positions 

def wyckoff_positions_c_2(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, y, -z]]) 
    return positions 

def wyckoff_positions_c_s(x, y, z): 
    positions = np.array([[x, y, z],
                          [x, -y, z]]) 
    return positions 

def wyckoff_positions_c_2h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, y, -z],
                          [-x, -y, -z],
                          [x, -y, z]]) 
    return positions 

def wyckoff_positions_d_2(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-x, y, -z],
                          [x, -y, -z]]) 
    return positions 

def wyckoff_positions_c_2v(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [x, -y, z],
                          [-x, y, z]]) 
    return positions 

def wyckoff_positions_d_2h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [-x, -y, -z],
                          [x, y, -z],
                          [x, -y, z],
                          [-x, y, z]]) 
    return positions 

def wyckoff_positions_c_4(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-y, x, z],
                          [y, -x, z]]) 
    return positions 

def wyckoff_positions_s_4(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [y, -x, -z],
                          [-y, x, -z]]) 
    return positions 

def wyckoff_positions_c_4h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-y, x, z],
                          [y, -x, z],
                          [-x, -y, -z],
                          [x, y, -z],
                          [y, -x, -z],
                          [-y, x, -z]]) 
    return positions 

def wyckoff_positions_d_4(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-y, x, z],
                          [y, -x, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [y, x, -z],
                          [-y, -x, -z]]) 
    return positions 

def wyckoff_positions_c_4v(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-y, x, z],
                          [y, -x, z],
                          [x, -y, z],
                          [-x, y, z],
                          [-y, -x, z],
                          [y, x, z]]) 
    return positions 

def wyckoff_positions_d_2d(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [y, -x, -z],
                          [-y, x, -z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [-y, -x, z],
                          [y, x, z]]) 
    return positions 

def wyckoff_positions_d_4h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-y, x, z],
                          [y, -x, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [y, x, -z],
                          [-y, -x, -z],
                          [-x, -y, -z],
                          [x, y, -z],
                          [y, -x, -z],
                          [-y, x, -z],
                          [x, -y, z],
                          [-x, y, z],
                          [-y, -x, z],
                          [y, x, z]]) 
    return positions 

def wyckoff_positions_c_3(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z]]) 
    return positions 

def wyckoff_positions_c_3i(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-x, -y, -z],
                          [y, -x+y, -z],
                          [x-y, x, -z]]) 
    return positions 

def wyckoff_positions_d_3(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-y, -x, -z],
                          [-x+y, y, -z],
                          [x, x-y, -z]]) 
    return positions 

def wyckoff_positions_c_3v(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-y, -x, z],
                          [-x+y, y, z],
                          [x, x-y, z]]) 
    return positions 

def wyckoff_positions_d_3d(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-y, -x, -z],
                          [-x+y, y, -z],
                          [x, x-y, -z],
                          [-x, -y, -z],
                          [y, -x+y, -z],
                          [x-y, x, -z],
                          [y, x, z],
                          [x-y, -y, z],
                          [-x, -x+y, z]]) 
    return positions 

def wyckoff_positions_c_6(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-x, -y, z],
                          [y, -x+y, z],
                          [x-y, x, z]]) 
    return positions 

def wyckoff_positions_c_3h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [x, y, -z],
                          [-y, x-y, -z],
                          [-x+y, -x, -z]]) 
    return positions 

def wyckoff_positions_c_6h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-x, -y, z],
                          [y, -x+y, z],
                          [x-y, x, z],
                          [-x, -y, -z],
                          [y, -x+y, -z],
                          [x-y, x, -z],
                          [x, y, -z],
                          [-y, x-y, -z],
                          [-x+y, -x, -z]]) 
    return positions 

def wyckoff_positions_d_6(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-x, -y, z],
                          [y, -x+y, z],
                          [x-y, x, z],
                          [y, x, -z],
                          [x-y, -y, -z],
                          [-x, -x+y, -z],
                          [-y, -x, -z],
                          [-x+y, y, -z],
                          [x, x-y, -z]]) 
    return positions 

def wyckoff_positions_c_6v(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-x, -y, z],
                          [y, -x+y, z],
                          [x-y, x, z],
                          [-y, -x, z],
                          [-x+y, y, z],
                          [x, x-y, z],
                          [y, x, z],
                          [x-y, -y, z],
                          [-x, -x+y, z]]) 
    return positions 

def wyckoff_positions_d_3h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [x, y, -z],
                          [-y, x-y, -z],
                          [-x+y, -x, -z],
                          [-y, -x, z],
                          [-x+y, y, z],
                          [x, x-y, z],
                          [-y, -x, -z],
                          [-x+y, y, -z],
                          [x, x-y, -z]]) 
    return positions 

def wyckoff_positions_d_6h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-y, x-y, z],
                          [-x+y, -x, z],
                          [-x, -y, z],
                          [y, -x+y, z],
                          [x-y, x, z],
                          [y, x, -z],
                          [x-y, -y, -z],
                          [-x, -x+y, -z],
                          [-y, -x, -z],
                          [-x+y, y, -z],
                          [x, x-y, -z],
                          [-x, -y, -z],
                          [y, -x+y, -z],
                          [x-y, x, -z],
                          [x, y, -z],
                          [-y, x-y, -z],
                          [-x+y, -x, -z],
                          [-y, -x, z],
                          [-x+y, y, z],
                          [x, x-y, z],
                          [y, x, z],
                          [x-y, -y, z],
                          [-x, -x+y, z]]) 
    return positions 

def wyckoff_positions_t(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [z, x, y],
                          [z, -x, -y],
                          [-z, -x, y],
                          [-z, x, -y],
                          [y, z, x],
                          [-y, z, -x],
                          [y, -z, -x],
                          [-y, -z, x]]) 
    return positions 

def wyckoff_positions_t_h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [z, x, y],
                          [z, -x, -y],
                          [-z, -x, y],
                          [-z, x, -y],
                          [y, z, x],
                          [-y, z, -x],
                          [y, -z, -x],
                          [-y, -z, x],
                          [-x, -y, -z],
                          [x, y, -z],
                          [x, -y, z],
                          [-x, y, z],
                          [-z, -x, -y],
                          [-z, x, y],
                          [z, x, -y],
                          [z, -x, y],
                          [-y, -z, -x],
                          [y, -z, x],
                          [-y, z, x],
                          [y, z, -x]]) 
    return positions 

def wyckoff_positions_o(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [z, x, y],
                          [z, -x, -y],
                          [-z, -x, y],
                          [-z, x, -y],
                          [y, z, x],
                          [-y, z, -x],
                          [y, -z, -x],
                          [-y, -z, x],
                          [y, x, -z],
                          [-y, -x, -z],
                          [y, -x, z],
                          [-y, x, z],
                          [x, z, -y],
                          [-x, z, y],
                          [-x, -z, -y],
                          [x, -z, y],
                          [z, y, -x],
                          [z, -y, x],
                          [-z, y, x],
                          [-z, -y, -x]]) 
    return positions 

def wyckoff_positions_t_d(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [z, x, y],
                          [z, -x, -y],
                          [-z, -x, y],
                          [-z, x, -y],
                          [y, z, x],
                          [-y, z, -x],
                          [y, -z, -x],
                          [-y, -z, x],
                          [y, x, z],
                          [-y, -x, z],
                          [y, -x, -z],
                          [-y, x, -z],
                          [x, z, y],
                          [-x, z, -y],
                          [-x, -z, y],
                          [x, -z, -y],
                          [z, y, x],
                          [z, -y, -x],
                          [-z, y, -x],
                          [-z, -y, x]]) 
    return positions 

def wyckoff_positions_o_h(x, y, z): 
    positions = np.array([[x, y, z],
                          [-x, -y, z],
                          [-x, y, -z],
                          [x, -y, -z],
                          [z, x, y],
                          [z, -x, -y],
                          [-z, -x, y],
                          [-z, x, -y],
                          [y, z, x],
                          [-y, z, -x],
                          [y, -z, -x],
                          [-y, -z, x],
                          [y, x, -z],
                          [-y, -x, -z],
                          [y, -x, z],
                          [-y, x, z],
                          [x, z, -y],
                          [-x, z, y],
                          [-x, -z, -y],
                          [x, -z, y],
                          [z, y, -x],
                          [z, -y, x],
                          [-z, y, x],
                          [-z, -y, -x],
                          [-x, -y, -z],
                          [x, y, -z],
                          [x, -y, z],
                          [-x, y, z],
                          [-z, -x, -y],
                          [-z, x, y],
                          [z, x, -y],
                          [z, -x, y],
                          [-y, -z, -x],
                          [y, -z, x],
                          [-y, z, x],
                          [y, z, -x],
                          [-y, -x, z],
                          [y, x, z],
                          [-y, x, -z],
                          [y, -x, -z],
                          [-x, -z, y],
                          [x, -z, -y],
                          [x, z, y],
                          [-x, z, -y],
                          [-z, -y, x],
                          [-z, y, -x],
                          [z, -y, -x],
                          [z, y, x]]) 
    return positions 

class PointGroup(enum.Enum):
    C_1 = enum.auto
    C_i = enum.auto
    C_2 = enum.auto
    C_s = enum.auto
    C_2h = enum.auto
    D_2 = enum.auto
    C_2v = enum.auto
    D_2h = enum.auto
    C_4 = enum.auto
    S_4 = enum.auto
    C_4h = enum.auto
    D_4 = enum.auto
    C_4v = enum.auto
    D_2d = enum.auto
    D_4h = enum.auto
    C_3 = enum.auto
    C_3i = enum.auto
    D_3 = enum.auto
    C_3v = enum.auto
    D_3d = enum.auto
    C_6 = enum.auto
    C_3h = enum.auto
    C_6h = enum.auto
    D_6 = enum.auto
    C_6v = enum.auto
    D_3h = enum.auto
    D_6h = enum.auto
    T = enum.auto
    T_h = enum.auto
    O = enum.auto
    T_d = enum.auto
    O_h = enum.auto

def wyckoff_positions(point_group: PointGroup, point):
    assert len(point) == 3
    x = point[0]
    y = point[1]
    z = point[2]

    print(point_group, type(point_group))
    print()

    switch = {
        PointGroup.C_1: wyckoff_positions_c_1,
        PointGroup.C_i: wyckoff_positions_c_i,
        PointGroup.C_2: wyckoff_positions_c_2,
        PointGroup.C_s: wyckoff_positions_c_s,
        PointGroup.C_2h: wyckoff_positions_c_2h,
        PointGroup.D_2: wyckoff_positions_d_2,
        PointGroup.C_2v: wyckoff_positions_c_2v,
        PointGroup.D_2h: wyckoff_positions_d_2h,
        PointGroup.C_4: wyckoff_positions_c_4,
        PointGroup.S_4: wyckoff_positions_s_4,
        PointGroup.C_4h: wyckoff_positions_c_4h,
        PointGroup.D_4: wyckoff_positions_d_4,
        PointGroup.C_4v: wyckoff_positions_c_4v,
        PointGroup.D_2d: wyckoff_positions_d_2d,
        PointGroup.D_4h: wyckoff_positions_d_4h,
        PointGroup.C_3: wyckoff_positions_c_3,
        PointGroup.C_3i: wyckoff_positions_c_3i,
        PointGroup.D_3: wyckoff_positions_d_3,
        PointGroup.C_3v: wyckoff_positions_c_3v,
        PointGroup.D_3d: wyckoff_positions_d_3d,
        PointGroup.C_6: wyckoff_positions_c_6,
        PointGroup.C_3h: wyckoff_positions_c_3h,
        PointGroup.C_6h: wyckoff_positions_c_6h,
        PointGroup.D_6: wyckoff_positions_d_6,
        PointGroup.C_6v: wyckoff_positions_c_6v,
        PointGroup.D_3h: wyckoff_positions_d_3h,
        PointGroup.D_6h: wyckoff_positions_d_6h,
        PointGroup.T: wyckoff_positions_t,
        PointGroup.T_h: wyckoff_positions_t_h,
        PointGroup.O: wyckoff_positions_o,
        PointGroup.T_d: wyckoff_positions_t_d,
        PointGroup.O_h: wyckoff_positions_o_h,
        }

    print(switch[point_group])

    return switch.get(point_group)(x,y,z)

