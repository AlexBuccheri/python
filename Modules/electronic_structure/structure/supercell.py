import numpy as np 

from . import atoms

class CellMinMax:
    def __init__(self, min, max):
        self.min = min
        self.max = max


def supercell_limits(n, centred_on_zero = False):
    limits = []
    if centred_on_zero:
        for i in range(0,3):
            n_mid = int(0.5 * n[i])
            limits.append(CellMinMax(-n_mid, n_mid+1))

    if not centred_on_zero:
        for i in range(0,3):
            limits.append(CellMinMax(0, n[i]))

    return limits


# For a number of translations n = [nx,ny,nz], find the flattened index associated
# with a cell defined by indices m = [mx,my,mz], where m_i \in n_i
#
# m: Indices of cell to flatten
# n: Total number of translations in each dimension
# centred_on_zero: Is the supercell centred on the translation (and hence indices) (0,0,0) ?
#
#  TODO(Alex) Write a smarter way of doing this with analytic expressions
#  Can compare to this as a point of reference
#
def flatten_supercell_limits(m, n, centred_on_zero = False):
    assert(len(m) == 3)
    assert(len(n) == 3)

    index_map = np.zeros(shape=(n[0],n[1],n[2]), dtype=int)
    flattened_index = 0
    for iz in range(0, n[2]):
        for iy in range(0, n[1]):
            for ix in range(0, n[0]):
                index_map[ix,iy,iz] = flattened_index
                flattened_index += 1

    ix = m[0]
    iy = m[1]
    iz = m[2]

    if centred_on_zero:
        limits = supercell_limits(n, centred_on_zero=centred_on_zero)
        ix += -limits[0].min
        iy += -limits[1].min
        iz += -limits[2].min

    return index_map[ix,iy,iz]

# Lattice vectors stored columnwise
def translation_vectors(lattice, n, centred_on_zero = False):

    limits = supercell_limits(n, centred_on_zero = centred_on_zero)
    translations = []

    for iz in range(limits[2].min, limits[2].max):
        for iy in range(limits[1].min, limits[1].max):
            for ix in range(limits[0].min, limits[0].max):
                translations.append(np.matmul(lattice, np.array([ix, iy, iz])))
    return translations


# unit_cell = atomic positions in fractional or angstrom
# translations = list of translation vectors in same units
def build_supercell(unit_cell, translations):
    supercell = []
    for translation in translations:
        for atom in unit_cell:
            position = translation + atom.position
            supercell.append(atoms.Atom(atom.species, position))
    return supercell


# TODO(Alex) Give this funciton a better name
# cells = [cell_1, cell_2, ..., cell_n] where cell_i contains global index for
# each basis atom in the cell  and n =  n1*n2*n3 is the total number of
# unit cells in the supercell
def list_global_atom_indices_per_cells(unit_cell, translations):
    cells = []
    iatom = 0
    for translation in translations:
        cell = []
        for atom in unit_cell:
            cell.append(iatom)
            iatom += 1
        cells.append(cell)

    return cells

# Lattice vectors define the dimensions & shape of the box
def molecule_in_supercell(lattice):
    supercell = []
    return supercell 
