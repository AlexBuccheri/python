# -----------------------------------------------------
# Primitive lattice vectors for each Bravais lattice
# 14 Bravais lattices with 7 crystal classes
# All defined column-wise:
# lattice = [[ax bx cx]
#            [ay by cy]
#            [az bz cz]]
#
# TODO(Alex) Add wrappers to accept object of lattice type
# Document LatticeOpt
# -----------------------------------------------------

import numpy as np


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
 

# ------------------------------------
# Triclinic
# ------------------------------------
def simple_triclinic(a,b,c, alpha, beta, gamma):
   """
    Triclinic lattice

    Given cell constants (a,b,c) and angles (alpha, beta, gamma), define lattice vectors
    for the triclinic lattice, stored columnwise.

    Parameters
    ----------
    a, b, c : double
              lattice constants
    alpha, beta, gamma  : double
                          cell angle, in radians
    Returns
    -------
    lattice : ndarray
              Matrix containing lattice vectors .shape(3,3)

    Notes
    -----
    Most general system: Primitive lattice vectors are completely general
    http://aflowlib.org/CrystalDatabase/triclinic_lattice.html
   """
   bx = b * np.cos(gamma)
   by = b * np.sin(gamma)
   cx = c * np.cos(beta)
   cy = c * (np.cos(alpha) - (np.cos(beta) * np.cos(gamma))) / np.sin(gamma)
   cz = np.sqrt(c*c - cx*cx - cy*cy)
   if c * c <= cx * cx + cy * cy:
       raise Exception("erroneous choice of triclinic cell parameter c and corresponding "
                       "angles, resulting in c*c <= cx*cx + cy*cy")

   return np.array([[a, bx, cx],
                    [0, by, cy],
                    [0,  0, cz]])

# ------------------------------------
# Monoclinic
# ------------------------------------
def simple_monoclinic(a,b,c, beta):
    """
     Simple monoclinic lattice

     Given cell constants (a,b,c) and angle (beta), define lattice vectors
     for the simple monoclinic lattice, stored columnwise.

     Parameters
     ----------
     a,b,c : double
             lattice constants
     beta  : double
             cell angle, in radians
     Returns
     -------
     lattice : ndarray
               Matrix containing lattice vectors .shape(3,3)

     Notes
     -----
      http://aflowlib.org/CrystalDatabase/monoclinic_lattice.html
     """
    cx = c * np.cos(beta)
    cz = c * np.sin(beta)
    cell = np.array( [[a, 0, cx],
                      [0, b,  0],
                      [0, 0, cz]])
    return cell


def base_centred_monoclinic(a, b, c, beta):
    cx = c * np.cos(beta)
    cz = c * np.sin(beta)
    cell = np.array([[ 0.5 * a,  0.5 * a, cx],
                     [-0.5 * b,  0.5 * b,  0],
                     [ 0,        0,       cz]])
    return cell



# ----------------------------------------------
# Orthorhombic
# ----------------------------------------------
def simple_orthorhombic(a,b,c):
    """
     Simple Orthorhombic lattice

     Given cell constants (a,b,c) define lattice vectors
     for the simple orthorhombic  lattice, stored columnwise.

     Parameters
     ----------
     a,b,c : double
             lattice constants

     Returns
     -------
     lattice : ndarray
               Matrix containing lattice vectors .shape(3,3)

     Notes
     -----
       http://aflowlib.org/CrystalDatabase/orthorhombic_lattice.html
     """
    assert a < b
    assert b < c
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

def base_centred_orthorhombic(spacegroup: str, a,b,c):
    if spacegroup[0].lower() =='a':
        return base_centred_orthorhombic_A(a,b,c)
    elif spacegroup[0].lower() =='c':
        return base_centred_orthorhombic_C(a,b,c)
    else:
        raise Exception('Space group must start with A or C')

def body_centred_orthorhombic(a,b,c):
    cell = 0.5 * np.array([[-a,  a,  a],
                           [ b, -b,  b],
                           [ c,  c, -c]])
    return cell

def face_centred_orthorhombic(a,b,c):
    cell = 0.5 * np.array([[ 0,  a,  a],
                           [ b,  0,  b],
                           [ c,  c,  0]])
    return cell


# -----------------------
# Tetragonal
# -----------------------
def simple_tetragonal(a, c):
    """
     Simple Tetragonal lattice

     Given cell constants (a, c) define lattice vectors
     for the simple tetragonal lattice, stored columnwise.

     Parameters
     ----------
     a, c : double
            lattice constants

     Returns
     -------
     lattice : ndarray
               Matrix containing lattice vectors .shape(3,3)

     Notes
     -----
       http://aflowlib.org/CrystalDatabase/tetragonal_lattice.html
     """
    cell = np.array([[a, 0, 0],
                     [0, a, 0],
                     [0, 0, c]])
    return cell

def body_centred_tetragonal(a, c):
    cell = 0.5 * a * np.array([[-1,     1,    1],
                               [ 1,    -1,    1],
                               [ c/a, c/a, -c/a]])
    return cell

# -----------------------
# Trigonal/Hexagonal
# -----------------------

# Adapted from Electronic structure (2004). Martin pg
# 2 triangular lattices
def graphene(a, c):
    sqrt3 = np.sqrt(3.)
    cell = a * np.array([[1, 0.5,          0],
                         [0, 0.5 * sqrt3,  0],
                         [0, 0,          c/a]])
    return cell


def hexagonal(a, c):
    """
     Simple trigonal or hexagonal lattice

     Given cell constants (a, c) define lattice vectors
     for the hexagonal lattice, stored columnwise.

     Parameters
     ----------
     a, c : double
            lattice constants

     Returns
     -------
     lattice : ndarray
               Matrix containing lattice vectors .shape(3,3)

     Notes
     -----
      It shares the same primitive vectors, but not point operations, as the hexagonal crystal system
      NOTE: I find this statement to be confusing
      http://aflowlib.org/CrystalDatabase/trigonal_lattice.html
     """
    cell = a * np.array([[ 0.5,              0.5,                0],
                         [-0.5 * np.sqrt(3), 0.5 * np.sqrt(3),   0],
                         [ 0,                0,                c/a]])
    return cell


def simple_trigonal(a, c):
    return hexagonal(a, c)


# Simple rhombohedral
# http://aflowlib.org/CrystalDatabase/trigonal_lattice.html
def rhombohedral_hex_setting(a, c):
    """
     Rhombohedral lattice in the hexgonal setting

     Given cell constants (a, c) define lattice vectors
     for the rhombohedral lattice in the hexgaonal setting, stored columnwise.

     Parameters
     ----------
     a, c : double
            lattice constants

     Returns
     -------
     lattice : ndarray
               Matrix containing lattice vectors .shape(3,3)

     Notes
     -----
      The International Tables address this ambiguity by listing atomic positions
      for the rhombohedral lattice in a hexagonal setting, where all coordinates are
      referenced to the conventional cell, and in a rhombohedral setting, where the
      coordinates are referenced to rhombohedral lattice
      http://aflowlib.org/CrystalDatabase/trigonal_lattice.html
     """
    sqrt3 = np.sqrt(3.)
    cell = np.array( [[ 0.5*a,       0,       -0.5*a       ],
                      [-0.5*a/sqrt3, a/sqrt3, -0.5*a/sqrt3 ],
                      [ c/3,         c/sqrt3,      c/sqrt3 ]])
    return cell

def rhombohedral_rhom_setting(a, alpha):
    """
     Rhombohedral lattice in the rhombohedral setting

     Given cell constant (a) and cell angle (alpha) define lattice vectors
     for the rhombohedral lattice in the rhombohedral setting, stored columnwise.

     Parameters
     ----------
     a, c : double
            lattice constants

     Returns
     -------
     lattice : ndarray
               Matrix containing lattice vectors .shape(3,3)

     Notes
     -----
      The International Tables address this ambiguity by listing atomic positions
      for the rhombohedral lattice in a hexagonal setting, where all coordinates are
      referenced to the conventional cell, and in a rhombohedral setting, where the
      coordinates are referenced to rhombohedral lattice.

      Defined according to A.11 in <a href="https://doi.org/10.1016/j.commatsci.2010.05.010">
      High-throughput electronic band structure calculations: Challenges and tools</a>, which
      differs from the expression in
      <a href="http://aflowlib.org/CrystalDatabase/trigonal_lattice.html">AFLOW trigonal lattice</a>,
      only by the orientation of the vectors relative to the Cartesian axes
     """
    ax = a * np.cos(0.5 * alpha);
    ay = a * np.sin(0.5 * alpha);
    cx = a * np.cos(alpha) / np.cos(0.5 * alpha)
    cz = a * np.sqrt(1 - (np.cos(alpha) * np.cos(alpha)) / (np.cos(0.5 * alpha) * np.cos(0.5 * alpha)));
    cell = np.array([[ ax,  ax,  cx],
                     [-ay,  ay,   0],
                     [ 0,    0,  cz]]);
    return cell

def rhombohedral(a, second_param, setting):
    if setting == 'hex':
        c = second_param
        return rhombohedral_hex_setting(a, c)
    elif setting == 'rhom':
        alpha = second_param
        return rhombohedral_rhom_setting(a, alpha)
    else:
        raise Exception("setting = 'hex' or 'rhom'")


# -----------------------
# Cubic
# -----------------------
def simple_cubic(a):
    """
     Simple cubic lattice

     Given cell constant (a) define lattice vectors
     for the simple cubic lattice, stored columnwise.

     Parameters
     ----------
     a : double
         lattice constants

     Returns
     -------
     lattice : ndarray
               Matrix containing lattice vectors .shape(3,3)

     Notes
     -----
       http://aflowlib.org/CrystalDatabase/cubic_lattice.html
     """
    cell = a * np.array( [[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
    return cell

def body_centred_cubic(a):
    cell = 0.5 * a * np.array([[-1,  1,  1],
                               [ 1, -1,  1],
                               [ 1,  1, -1]])
    return cell

def face_centred_cubic(a):
    cell = 0.5 * a * np.array([[0, 1, 1],
                               [1, 0, 1],
                               [1, 1, 0]])
    return cell
