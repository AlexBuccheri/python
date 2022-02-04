"""
Generate consistent band paths
"""
import ase
import numpy as np

def k_path_with_ase(lattice_vectors):
    """

    lattice_vectors in angstrom, for example:
    [[0., 2.7145, 2.7145], [2.7145, 0., 2.7145], [2.7145, 2.7145, 0.]]
    https://wiki.fysik.dtu.dk/ase/gettingstarted/tut04_bulk/bulk.html
    :return:
    """
    cell = ase.atoms.Cell(lattice_vectors)
    path = cell.bandpath()
    # special_points = path.special_points
    # k_path = path.path
    # print(special_points)
    # print(k_path)
    # print(path.kpts)
    return path


#   KPointsAndWeights = Klines {
#      1   0.375  0.375  0.75   # K
#     20   0.0    0.0    0.0    # Gamma
#     20   0.5    0.5    0.5    # L
#   }

def dftb_kpoints_str(ase_band_path) -> str:
    """ Convert explicit k-points for bath through BZ, to a DFTB+ input

    :param ase_band_path:
    :return: DFTB+ KPointsAndWeights string
    """
    string = "KPointsAndWeights = {\n"
    n_k_points = ase_band_path.kpts.shape[0]
    weight = "1.0"

    for ik in range(0, n_k_points):
        k_point = ase_band_path.kpts[ik]
        k_str = np.array2string(k_point, precision=8, separator=' ', suppress_small=False)[1:-1]
        string += "  " + k_str + " " + weight + "\n"

    string += "}\n"

    return string

# Should confirm this works for the band structure - put in manually, then plot
path = k_path_with_ase([[0., 2.7145, 2.7145], [2.7145, 0., 2.7145], [2.7145, 2.7145, 0.]])
print(dftb_kpoints_str(path))