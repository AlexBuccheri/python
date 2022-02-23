"""Generate Bulk Band Structures using DFTB+ (where possible) and TB Lite (xTB1, xTB2).

* Create directories
* Copy geometry and charges
* Define band structure input, with consistent band path
* Run
* Convert file format and parse result
"""
import os.path
import numpy as np
from typing import List, Union, Optional
from pathlib import Path
import shutil
import json
import matplotlib.pyplot as plt

import ase
from ase.io.gen import read_gen
from ase.atoms import Atoms

from tb_lite.src.dftb_input import generate_band_structure_input
from tb_lite.src.parsers import parse_dftb_output, parse_dftb_bands, shift_bands

# Path type
path_type = Union[Path, str]


def generate_band_structure_inputs(scc_to_band_directory: dict):
    """ Given a list of directories with converged calculations, generate band structure inputs.
    """
    for scc_directory, bs_directory in scc_to_band_directory.items():
        # Make directory if it does not exist
        Path(bs_directory).mkdir(parents=True, exist_ok=True)

        # Copy the charges and structure
        for file in ['charges.bin', 'geometry.gen']:
            shutil.copyfile(os.path.join(scc_directory, file), os.path.join(bs_directory, file))

        # Generate a new input file, for band structure
        atoms: Atoms = read_gen(os.path.join(bs_directory, 'geometry.gen'))
        input_xml_str = generate_band_structure_input(atoms.get_cell(), 'GFN1-xTB')

        with open(os.path.join(bs_directory, 'dftb_in.hsd'), 'w') as fid:
            fid.write(input_xml_str)


def get_standardised_band_path(lattice_vectors):
    """
    FOR RUI

    ASE API for getting a standard band path and a fixed k-grid
    sampling the path.

    Some more methods to try:  band_path.plot()


    :param lattice_vectors:
    :return:
    """
    cell = ase.atoms.Cell(lattice_vectors)
    band_path: ase.dft.kpoints.BandPath = cell.bandpath()
    k_points = band_path.kpts
    high_sym_points = band_path.special_points

    return k_points, high_sym_points



def flatten_k_path(k_points):
    """

    :param k_points:
    :return:
    """

    n_k_points = k_points.shape[0]
    flat_path = np.empty(shape=n_k_points)
    flat_path[0] = 0.

    for ik in range(1, n_k_points):
        dk = k_points[ik, :] - k_points[ik - 1, :]
        flat_path[ik] = flat_path[ik-1] + np.linalg.norm(dk)

    return flat_path



def plot_band_structure(k_points, high_sym_points, bands: np.ndarray, n_occupied: Optional[int]=None):
    """
        c) Band structure, and return in a plottable format
        d) Label path using ASE output
    """
    fig, ax = plt.subplots(figsize=(6, 9))

    n_bands = bands.shape[1]
    n_k_points = bands.shape[0]

    # Assign hugh symmetry symbols to k indices
    symbols = []
    all_indices = []
    for symbol, point in high_sym_points.items():
        mask = (k_points == point).all(-1)
        indices = np.where(mask == True)[0].tolist()
        all_indices += indices
        symbols += [symbol] * len(indices)

    all_indices = np.asarray(all_indices)
    symbols = np.asarray(symbols)

    sort = np.argsort(all_indices)
    sorted_indices = all_indices[sort]

    # Looks like one might be able to gey away with np.arange(0, n_k_points)
    k = flatten_k_path(k_points)

    ax.set_xticks(k[sorted_indices])
    ax.set_xticklabels(symbols[sort])

    for i in range(0, n_bands):
        plt.plot(k, bands[:, i])
    plt.show()




def number_of_occupied_bands(detailed_out_str: str) -> float:
    """Parse the number of electrons, and if spin-polarised

    :return: Number of electrons
    """
    spin_polarised = False
    electrons_per_band = 1 if spin_polarised else 2
    n_electrons = parse_dftb_output(detailed_out_str)['n_electrons_up']
    n_occupied_bands = float(n_electrons) / float(electrons_per_band)
    return n_occupied_bands


# TODO Test me
class BandGap:
    """ Compute band edges and gap.
    """

    def __init__(self, n_occupied_bands: Union[int, float], band_energies: np.ndarray):
        """ Initialise Bandgap class.

        :param n_occupied_bands: Number of occupied bands. Should be an integer
        One could replace with the index of the highest-occupied band (either way, it becomes
        ambiguous at finite temperatures).
        :param band_energies: Numpy array of shape(n_k_points, n_bands).
        """
        int_n_occupied_bands = int(n_occupied_bands)
        self.n_occupied_bands = int_n_occupied_bands
        self.band_energies = band_energies
        if int_n_occupied_bands != n_occupied_bands:
            raise ValueError('Number of occupied bands is not an integer. \n'
                             'Number of electrons is likely not an integer value.')

    def band_edge_energies(self) -> tuple:
        # Python indexing at 0
        E_vb_max = np.max(self.band_energies[:, self.n_occupied_bands - 1])
        E_cb_min = np.min(self.band_energies[:, self.n_occupied_bands])
        return E_vb_max, E_cb_min

    def band_gap(self) -> float:
        E_vb_max, E_cb_min = self.band_edge_energies()
        return E_cb_min - E_vb_max

    def band_edge_k_points(self, k_points: np.ndarray) -> tuple:
        """ Get the k-points associated with the band edges.
        TODO Test me
        :param k_points: np array of k-points.
        :return: Tuple of k-points at the VBM and CBm, respectively.
        """
        assert k_points.shape[1] == 3, "Expect k_points.shape == (n_k_points, 3)"
        i_vb_max = np.argmax(self.band_energies[:, self.n_occupied_bands - 1])
        i_cb_min = np.argmin(self.band_energies[:, self.n_occupied_bands])
        return k_points[i_vb_max, :], k_points[i_cb_min, :]


# Try these to ensure all the post-processing works
# bs_root = '/Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/band_structures/'
# converged_dirs = [bs_root + x for x in ['diamond', 'silicon', 'germanium']]
# generate_band_structure_inputs(converged_dirs)

with open('/Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/band_structures/diamond_bands/detailed.out',
          "r") as fid:
    detailed_str = fid.read()
data = parse_dftb_output(detailed_str)
n = number_of_occupied_bands(detailed_str)
print(n)
band_energies = parse_dftb_bands(
    '/Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/band_structures/diamond_bands/')
# print(band_energies)
band_details = BandGap(n, band_energies)
vmax, cmin = band_details.band_edge_energies()
print(vmax, cmin, band_details.band_gap())
shifted_bands = shift_bands(band_energies, vmax)

lattice_vectors_diamond = np.array([[0.000000000000000, 1.783500000000000, 1.783500000000000],
                                    [1.783500000000000, 0.000000000000000, 1.783500000000000],
                                    [1.783500000000000, 1.783500000000000, 0.000000000000000]])

k_points, high_sym_points = get_standardised_band_path(lattice_vectors_diamond)

plot_band_structure(k_points, high_sym_points, shifted_bands, [])


# if __name__ == "__main__":
#     # All operations required to produce band structure plots.
#
#     with open('converged_energies.json', 'r', encoding='utf-8') as fid:
#         scc_data: dict = json.load(fid)
#
#     # Directory pairs relating SCC to band structure
#     bands_root = "/home/alex/tblite/bulk/bands"
#     scc_to_band_directory = {}
#     for scc in scc_data.values():
#         scc_dir = scc['directory']
#         material_name = os.path.basename(scc_dir)
#         scc_to_band_directory[scc_dir] = os.path.join(bands_root, material_name)
#
#     # Generate inputs
#     generate_band_structure_inputs(scc_to_band_directory)
#
#     # Run inputs
#
#     # Convert outputs to a useful format
#
#     # Print band gap information
#
#
#     # Parse the useful BS format
#     # Plot band structure
#     # save to file
