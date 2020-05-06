# Primitive lattice vectors for each Bravais lattice
# 14 Bravais lattices with 7 crystal classes
# All defined column-wise

import numpy as np


#TODO(Alex) Add body-centred orthorhomic
#TODO(Alex) Add face-centred orthorhomic
#TODO(Alex) Add other missing lattices and make consistent with the aflow paper
#Add wrappers to accept object of lattice type




class LatticeOpt:
    def __init__(self, bravais_type, a=None, b=None, c=None, alpha=None, beta=None, gamma=None, units='angstrom'):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.bravais_type = bravais_type
        # Convert to ENUM class
        if units != 'angstrom':
            self.units = units
 

# ----------------------------------------------------------------------
# Triclinic
# Most general system: Primitive lattice vectors are completely general
# http://aflowlib.org/CrystalDatabase/triclinic_lattice.html
# ----------------------------------------------------------------------
def simple_triclinic(a,b,c, alpha, beta, gamma):
    #TODO(Alex) Add assert used in entos
   bx = b * np.cos(gamma)
   by = b * np.sin(gamma)
   cx = c * np.cos(beta)
   cy = c * (np.cos(alpha) - (np.cos(beta) * np.cos(gamma))) / np.sin(gamma)
   cz = np.sqrt(c*c - cx*cx - cy*cy)
   return   np.array([[a, bx, cx],
                     [0, by, cy],
                     [0,  0, cz]])

# ----------------------------------------------
# Monoclinic
# ----------------------------------------------
def simple_monoclinic(a,b,c, beta):
    cx = c * np.cos(beta)
    cz = c * np.sin(beta)
    cell = np.array( [[a, 0, cx],
                      [0, b,  0],
                      [0, 0, cz]])
    return cell

def base_centred_monoclinic(a, b, c, beta):
    cx = c * np.cos(beta)
    cz = c * np.sin(beta)
    cell = np.array([[ 0.5*a,-0.5*a, cx],
                     [-0.5*b, 0.5*b,  0],
                     [ 0,     0,     cz]])
    return cell

# ----------------------------------------------
# Orthorhombic
# http://aflowlib.org/CrystalDatabase/orthorhombic_lattice.html
# ----------------------------------------------
def simple_orthorhombic(a,b,c):
    cell = np.array( [[a, 0, 0],
                      [0, b, 0],
                      [0, 0, c]])
    return cell

def base_centred_orthorhombic_A(a,b,c):
    cell = np.array( [[a,   0.     , 0.     ],
                      [0.,  0.5 * b, 0.5 * b],
                      [0., -0.5 * c, 0.5 * c]]);
    return cell

def base_centred_orthorhombic_C(a,b,c):
    cell = np.array( [[ 0.5 * a,  0.5 * a,  0.],
                      [-0.5 * b,  0.5 * b,  0.],
                      [ 0.     ,  0.     ,  c]]);
    return cell

def base_centred_orthorhombic(spacegroup, a,b,c):
    if spacegroup[0] =='A':
        return base_centred_orthorhombic_A(a,b,c)
    elif spacegroup[0] =='C':
        return base_centred_orthorhombic_C(a,b,c)
    else:
        exit('Space group must start with A or C')


# def body_centred_orthorhombic(self, a,b,c):
#     return cell

# def face_centred_orthorhombic(self, a,b,c):
#     return cell


# -----------------------
# Tetragonal
# -----------------------
def simple_tetragonal(a, c):
    cell = np.array([[a, 0, 0],
                     [0, a, 0],
                     [0, 0, c]])
    return cell

# Body-centred tetragonal
# http://aflowlib.org/CrystalDatabase/tetragonal_lattice.html
def body_centred_tetragonal(a, c):
    cell = 0.5 * a * np.array([[-1,     1,    1],
                               [ 1,    -1,    1],
                               [ c/a, c/a, -c/a]])
    return cell

# -----------------------
# Trigonal or Hexagonal
# -----------------------
# Adapted from Electronic structure (2004). Martin pg
# 2 triangular lattices
def graphene(a, c):
    sqrt3 = np.sqrt(3.)
    cell = a * np.array([[1, 0.5,          0],
                         [0, 0.5 * sqrt3,  0],
                         [0, 0,          c/a]])
    return cell


# Simple hexagonal or trigonal
# http://aflowlib.org/CrystalDatabase/trigonal_lattice.html
# 3 triangular lattices


# Simple rhombohedral
# http://aflowlib.org/CrystalDatabase/trigonal_lattice.html
def simple_rhombohedral(a, c):
    sqrt3 = np.sqrt(3.)
    cell = np.array( [[ 0.5*a,       0,       -0.5*a       ],
                      [-0.5*a/sqrt3, a/sqrt3, -0.5*a/sqrt3 ],
                      [ c/3,         c/sqrt3,      c/sqrt3 ]])
    return cell

# -----------------------
# Cubic
# -----------------------
# Simple or primitive cubic cell
def simple_cubic(a):
    cell = a * np.array( [[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
    return cell

# Body-centred cubic
def body_centred_cubic(a):
    cell = 0.5 * a * np.array([[-1,  1,  1],
                               [ 1, -1,  1],
                               [ 1,  1, -1]])
    return cell

# Face-centred cubic
def face_centred_cubic(a):
    cell = 0.5 * a * np.array([[0, 1, 1],
                               [1, 0, 1],
                               [1, 1, 0]])
    return cell
