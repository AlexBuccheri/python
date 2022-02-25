"""Generate Bulk Band Structures using DFTB+ (where possible) and TB Lite (xTB1, xTB2).

* Create directories
* Copy geometry and charges
* Define band structure input, with consistent band path
* Run
* Convert file format and parse result
"""
import os.path
import pathlib

import numpy as np
from typing import List, Union, Optional, Tuple
from pathlib import Path
import shutil
import json
import matplotlib.pyplot as plt

import ase
from ase.io.gen import read_gen
from ase.atoms import Atoms

from tb_lite.src.dftb_input import generate_band_structure_input
from tb_lite.src.parsers import parse_dftb_output, parse_dftb_bands, parse_number_of_occupied_bands, parse_geometry_gen

# Path type
path_type = Union[Path, str]


# --------------------------------------------
# Generate band structure files for TB Lite
# --------------------------------------------

def generate_band_structure_inputs(scc_to_band_directory: dict):
    """ Given a list of directories with converged calculations, generate band structure inputs
    for DFTB+ TB Lite.
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


# --------------------------------------------------------------------------------------
# Generic
# TODO(Alex) Give to Rui
# TODO NOTE! ASE expects lattice vectors row-wise. Confirm a) this is correct
#  and b)I've done this in all instances (see data/bulk_crystals)
# --------------------------------------------------------------------------------------

def get_standardised_band_path(lattice_vectors) -> Tuple[np.ndarray, dict]:
    """ ASE standardised band path and a fixed k-grid sampling the path.

    Notes:
      Check out band_path.plot()

    :param lattice_vectors: Lattice vectors stored row wise np array or as [a, b, c]
    :return: Tuple of k-points (n_k_points, 3) and high symmetry points {symbol: k_point/s}
    """
    cell = ase.atoms.Cell(lattice_vectors)
    band_path: ase.dft.kpoints.BandPath = cell.bandpath()
    return band_path.kpts, band_path.special_points


class BandGap:
    """ Compute band edges and gap.
    """

    def __init__(self, n_occupied_bands: Union[int, float], band_energies: np.ndarray):
        """ Initialise BandGap class.

        :param n_occupied_bands: Number of occupied bands. Should be an integer.
        One could replace with the index of the highest-occupied band (either way, it becomes
        ambiguous at finite temperatures).
        :param band_energies: Numpy array of shape(n_k_points, n_bands).
        """
        int_n_occupied_bands = int(n_occupied_bands)
        self.n_occupied_bands = int_n_occupied_bands
        # Python indexing at 0
        self.ivm = int_n_occupied_bands - 1
        self.band_energies = band_energies
        if int_n_occupied_bands != n_occupied_bands:
            raise ValueError('Number of occupied bands is not an integer. \n'
                             'Number of electrons is likely not an integer value.')

    def band_edge_energies(self) -> tuple:
        """TODO(Alex) Document
        :return:
        """
        E_vb_max = np.max(self.band_energies[:, self.ivm])
        E_cb_min = np.min(self.band_energies[:, self.ivm + 1])
        return E_vb_max, E_cb_min

    def band_gap(self) -> float:
        """ Band gap, defined as E_CBm - E_VBM
        :return: Band gap
        """
        E_vb_max, E_cb_min = self.band_edge_energies()
        return E_cb_min - E_vb_max

    def band_edge_k_points(self, k_points: np.ndarray) -> tuple:
        """ Get the k-points associated with the band edges.
        TODO Test me
        :param k_points: np array of k-points.
        :return: Tuple of k-points at the VBM and CBm, respectively.
        """
        assert k_points.shape[1] == 3, "Expect k_points.shape == (n_k_points, 3)"
        ik_vb_max = np.argmax(self.band_energies[:, self.ivm])
        ik_cb_min = np.argmin(self.band_energies[:, self.ivm + 1])
        return k_points[ik_vb_max, :], k_points[ik_cb_min, :]


def set_bands_zeropoint(bands: np.ndarray, zero_point: float) -> np.ndarray:
    """ Shift the bands to a new zero point.

    :param bands: Bands with shape(n_k_points, n_bands).
    :param zero_point: New zero point of the energy.
    :return: shifted_bands
    """
    shifted_bands = bands.copy()
    for ik in range(bands.shape[0]):
        shifted_bands[ik, :] -= zero_point
    return shifted_bands


def flatten_k_path(k_points: np.ndarray) -> np.ndarray:
    """ Given a set of k-points, return a vector of accumulated
    displacement norms.

    For use with band structure plotting.

    TODO Improve description.

    :param k_points: An array of k-points.
    :return: Flat path.
    """
    assert k_points.shape[1] == 3, "Expect k_points.shape = (n_k_points, 3)"

    n_k_points = k_points.shape[0]
    flat_path = np.empty(shape=n_k_points)

    flat_path[0] = 0.
    for ik in range(1, n_k_points):
        dk = k_points[ik, :] - k_points[ik - 1, :]
        flat_path[ik] = flat_path[ik - 1] + np.linalg.norm(dk)

    return flat_path


def xticks_and_xticklabels(k_points: np.ndarray, high_sym_points: dict) -> tuple:
    """
    Establish which k-points correspond to high symmetry points.
    For use with plotting.

    TODO(Alex) Tidy this up and rename

    :param k_points: k-points grid.
    :param high_sym_points: Dict of form {'G': np.array([0., 0., 0.]), ...}
    :return:
    """
    # Assign high symmetry symbols to k indices
    symbols = []
    all_indices = []
    for symbol, point in high_sym_points.items():
        mask = (k_points == point).all(-1)
        indices = np.where(mask == True)[0].tolist()
        all_indices += indices
        symbols += [symbol] * len(indices)

    all_indices = np.asarray(all_indices)
    symbols = np.asarray(symbols)

    in_order = np.argsort(all_indices)
    # Indices of k-points which correspond to high symmetry points, in order
    # Symbols associated with the high symmetry points, in consistent order
    return all_indices[in_order], symbols[in_order]


def plot_band_structure(k_points: np.ndarray, high_sym_points: dict, bands: np.ndarray):
    """ Plot a band structure.
    """
    fig, ax = plt.subplots(figsize=(6, 9))

    n_bands = bands.shape[1]
    high_sym_indices, high_sym_symbols = xticks_and_xticklabels(k_points, high_sym_points)

    # Looks like one might be able to gey away with np.arange(0, n_k_points)
    k = flatten_k_path(k_points)

    ax.set_xticks(k[high_sym_indices])
    ax.set_xticklabels(high_sym_symbols)
    plt.ylabel('Energy (ev)')

    for i in range(0, n_bands):
        plt.plot(k, bands[:, i])
    return fig, ax


def process_band_structure(directory: Union[str, pathlib.Path]):
    """ Take a DFTB+ band structure output and return a band structure plot.

    Parsing is DFTB+ specific, band structure generation is generic.
    """
    pj = os.path.join

    # DFTB+ main output
    with open(pj(directory, 'detailed.out'), "r") as fid:
        detailed_str = fid.read()

    # Geometry
    geo_data = parse_geometry_gen(pj(directory, 'geometry.gen'))

    # Band structure
    band_energies = parse_dftb_bands(directory)

    # Band gap details. eV as DFTB+ band structure parsed in eV
    n_occ = parse_number_of_occupied_bands(detailed_str)
    band_details = BandGap(n_occ, band_energies)
    vmax, cmin = band_details.band_edge_energies()
    print("VBM, CBm, and band gap (eV):", vmax, cmin, band_details.band_gap())

    # Band structure plot
    shifted_bands = set_bands_zeropoint(band_energies, vmax)
    k_points, high_sym_points = get_standardised_band_path(geo_data['lattice'])
    fig, ax = plot_band_structure(k_points, high_sym_points, shifted_bands)
    plt.show()


# TODO(Alex) Get thi working and clean
def some_example_of_generating_bands():
    with open('converged_energies.json', 'r', encoding='utf-8') as fid:
        scc_data: dict = json.load(fid)

    # Directory pairs relating SCC to band structure
    bands_root = "/home/alex/tblite/bulk/bands"
    scc_to_band_directory = {}
    for scc in scc_data.values():
        scc_dir = scc['directory']
        material_name = os.path.basename(scc_dir)
        scc_to_band_directory[scc_dir] = os.path.join(bands_root, material_name)

    # Generate inputs
    generate_band_structure_inputs(scc_to_band_directory)

    # Run inputs


process_band_structure('/Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/band_structures/diamond_bands/')