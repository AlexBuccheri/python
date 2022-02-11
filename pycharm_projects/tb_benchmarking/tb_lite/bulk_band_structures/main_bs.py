"""
Generate Bulk Band Structures using DFTB+ (where possible) and TB Lite (xTB1, xTB2)
for these materials:
* IV: Si, Ge, Diamond - GOT
* TODO  GET Wide-gap: ZnO, TiO2 (both phases), ZrO2, WO3
* TODO  GET BN (cubic and hexagonal)
* TODO  GET III-V: GaN, GaP, GaAs, InN, InP, InAs
* TODO  GET Narrow band-gap II-VI: PbS, PbSe, PbTe
* TODO  GET MoS2, WS2
"""
import os.path
import numpy as np
from typing import List
from pathlib import Path
import shutil

from ase.io.gen import read_gen
from ase.atoms import Atoms

from tb_lite.src.dftb_input import generate_band_structure_input


def generate_band_structure_inputs(converged_directories: List[str]):
    """ Given a list of directories with converged calculations,
        Generate band structure inputs.
    :return:
    """
    for converged_directory in converged_directories:
        # Create new directory - assumes a naming convention: some/directory/material_name
        head, material_name = os.path.split(converged_directory)
        bands_directory = os.path.join(head, material_name + "_bands")
        Path(bands_directory).mkdir(parents=True, exist_ok=True)

        # Copy the charges and structure
        for file in ['charges.bin', 'geometry.gen']:
            shutil.copyfile(os.path.join(converged_directory, file), os.path.join(bands_directory, file))

        # Generate a new input file, for band structure
        atoms: Atoms = read_gen(os.path.join(bands_directory, 'geometry.gen'))
        input_xml_str = generate_band_structure_input(atoms.get_cell(), 'GFN1-xTB')
        with open(os.path.join(bands_directory, 'dftb_in.hsd'), 'w') as fid:
            fid.write(input_xml_str)


def plot_band_structure():
    """
    Parse:
        a) Direct gap
        b) CBm - VBM
        c) Band structure, and return in a plottable format
        d) Label path using ASE output
    :return:
    """
    return None


def number_of_occupied_bands() -> float:
    """
    Parse the number of electrons, and if spin-polarised

    :return:
    """
    spin_polarised = False
    electrons_per_band = 2 if spin_polarised else 1
    n_electrons = 8  # TODO PARSE ME
    n_occupied_bands = float(n_electrons) / float(electrons_per_band)
    return n_occupied_bands


# TODO Test me
class BandGap:
    def __init__(self, n_occupied_bands, band_energies: np.ndarray):
        self.n_occupied_bands = n_occupied_bands
        self.band_energies = band_energies

    def band_edge_energies(self):
        E_vb_max = np.max(self.band_energies[:, self.n_occupied_bands])
        E_cb_min = np.min(self.band_energies[:, self.n_occupied_bands + 1])
        return E_vb_max, E_cb_min

    def band_gap(self):
        E_vb_max, E_cb_min = self.band_edge_energies()
        return E_cb_min - E_vb_max

    def band_edge_k_points(self, k_points):
        assert k_points.shape[1] == 3, "Expect k_points.shape == (n_k_points, 3)"
        # TODO Check the syntax
        i_vb_max = np.amax(self.band_energies[:, self.n_occupied_bands])
        i_cb_min = np.amin(self.band_energies[:, self.n_occupied_bands + 1])
        return k_points[i_vb_max, :], k_points[i_cb_min, :]


# bs_root = '/Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/band_structures/'
# converged_dirs = [bs_root + x for x in ['diamond', 'silicon', 'germanium']]
# generate_band_structure_inputs(converged_dirs)
