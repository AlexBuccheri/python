"""
Parsers
"""
import re
import os
import subprocess
import numpy as np

import ase
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure


def parse_qcore_structure(file_name: str) -> dict:
    """ Parse the atomic structure and lattice vectors from QCore.
    **Worked specifically in the context of the molecular crystal input files**

    TODO(Alex) Would be easier to get the positions from cif, using ASE or pymatgen
      Whole routine is some HACKY GARBAGE SHIT

    Each file has the same structure:
    mgo_volume := xtb(
      structure(
          fractional= [['H', 0.1988, 0.26105, 0.14368],
                       ['H', 0.30118, 0.76105, 0.35632],
                        .
                        .
                        .
                       ['O', 0.96018, 0.57252, 0.74478]]

          lattice(
                   lattice_vectors = [[3.5263673508626314, 0.0, -0.6223031990366413], [1.4474643463714529e-15, 9.0009574, 5.511496834585842e-16], [0.0, 0.0, 6.952313179919385]]
                   lattice_vectors_unit = angstrom
                  )

    :param file_name: File name
    :return: dict containing species of all atoms, atomic positions in fractional coordinates, lattice vectors and their
    units.
    """
    with open(file_name) as fid:
        qcore_input_lines = fid.readlines()

    # Find the structure block
    start = 0
    end = 0
    for i, line in enumerate(qcore_input_lines):
        line_str = line.strip()
        match_fractional = re.match("fractional", line_str)
        match_lattice = re.match("lattice", line_str)

        if match_fractional:
            start = i
        if match_lattice:
            end = i
            break

    # Remove fractional from first line
    # fractional= [['H', 0.1988, 0.26105, 0.14368],
    symbol_index = qcore_input_lines[start].find('=')
    qcore_input_lines[start] = qcore_input_lines[start][symbol_index + 3:]

    species = []
    fractional_positions = []
    for i in range(start, end - 1):
        # ['O', 0.96018, 0.57252, 0.74478]
        species_str, x, y, z = qcore_input_lines[i].strip()[1:-2].split(',')
        species.append(species_str.strip("\'"))
        fractional_positions.append([float(r) for r in [x, y, z]])

    # The rest can be extracted with re
    symbol_index = qcore_input_lines[end+1].find('=')
    qcore_input_lines[end + 1] = qcore_input_lines[end+1][symbol_index+1:]

    def remove_brackets(stuff):
        stuff = [r.replace('[', '') for r in stuff]
        stuff = [r.replace(']', '') for r in stuff]
        return stuff

    vectors = qcore_input_lines[end + 1].split(',')
    vectors = remove_brackets(vectors)
    a = [float(v) for v in vectors[0:3]]
    b = [float(v) for v in vectors[3:6]]
    c = [float(v) for v in vectors[6:9]]

    unit = qcore_input_lines[end + 2].split()[-1].strip()

    return {'species': species, 'fractional_positions': fractional_positions,
            'lattice': [a, b, c], 'lattice_vectors_unit': unit}


def parse_qcore_settings(file_name: str) -> dict:
    """ Parse some select calculation settings from Qcore input

    :param file_name: File name
    :return: dictionary of settings
    """
    with open(file_name) as fid:
        qcore_input_lines = fid.readlines()

    k_points = None
    electronic_temperature = None
    for line in qcore_input_lines:
        line_str = line.strip()
        match_mp = re.match("monkhorst_pack", line_str)
        match_etmp = re.match("temperature", line_str)
        if match_mp:
            k_str = line_str.split('=')[-1].strip()[1:-1]
            k_points = [float(k) for k in k_str.split(',')]
        if match_etmp:
            electronic_temperature = float(line_str.split()[-2])

    data = {'k_points': k_points, 'electronic_temperature': electronic_temperature}

    return data


def parse_dftb_output(input: str) -> dict:
    """
    Parse selected DFTB+ TB lite outputs (from detailed.out) in eV.

    End of file has the typical structure:
    '''
       Fermi level:                        -0.4495587982 H          -12.2331 eV
       Band energy:                       -21.9227015154 H         -596.5471 eV
       TS:                                  0.0000949682 H            0.0026 eV
       Band free energy (E-TS):           -21.9227964836 H         -596.5496 eV
       Extrapolated E(0K):                -21.9227489995 H         -596.5484 eV
       Input / Output electrons (q):     44.0000000000     44.0000000000

       Energy H0:                         -18.2466149624 H         -496.5157 eV
       Energy SCC:                         -0.0415167324 H           -1.1297 eV
       Total Electronic energy:           -18.2881316948 H         -497.6454 eV
       Repulsive energy:                    0.0000000000 H            0.0000 eV
       Total energy:                      -18.2881316948 H         -497.6454 eV
       Extrapolated to 0:                 -18.2881791788 H         -497.6467 eV
       Total Mermin free energy:          -18.2882266629 H         -497.6480 eV
       Force related energy:              -18.2882266629 H         -497.6480 eV

      SCC converged

      Full geometry written in geo_end.{xyz|gen}

      Geometry converged
    '''

    :param input: File string contents
    :return: results dictionary
    """
    energy_h0 = re.findall(r'^Energy H0: .*$', input, flags=re.MULTILINE)[0].split()[-2]
    energy_scc = re.findall(r'^Energy SCC: .*$', input, flags=re.MULTILINE)[0].split()[-2]
    total_electronic_energy = re.findall(r'^Total Electronic energy: .*$', input, flags=re.MULTILINE)[0].split()[-2]
    repulsive_energy = re.findall(r'^Repulsive energy: .*$', input, flags=re.MULTILINE)[0].split()[-2]
    total_energy = re.findall(r'^Total energy: .*$', input, flags=re.MULTILINE)[0].split()[-2]

    results = {'energy_h0': energy_h0, 'energy_scc': energy_scc, 'total_electronic_energy': total_electronic_energy,
               'repulsive_energy': repulsive_energy, 'total_energy': total_energy}

    # Convert all strings to floats
    for key, value in results.items():
        results[key] = float(value)

    return results


def parse_dftb_bands(directory) -> np.ndarray:
    """
    Return band energies from a DFTB+ TB Lite Output

    Need DFTB+ output to be in a sensible format, so apply this script: dp_bands band.out band
    Went to dftb-19/tools folder on my and installed dp_tools via the terminal, to this virtual env:
    `python setup.py install --prefix /Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/venv/`

    bands.shape = (n_k_points, n_bands),
    where n_k_points is the total number of k-points used across all bands paths

    :param directory: Directory containing band structure file
    :return: band energies
    """
    processed_file = os.path.join(directory, 'band_tot.dat')

    if not os.path.isfile(processed_file):
        result = subprocess.run(['dp_bands', 'band.out', 'band'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                cwd=directory)
        if result.returncode != 0:
            raise subprocess.SubprocessError(f'Not able to convert band.out: {result.stderr}')

    with open(processed_file, 'r') as fid:
        lines = fid.readlines()

    n_k_points = len(lines)
    n_bands = len(lines[0].split()[1:])

    # File signature:  ith_kpoint eigenvalue1, eigenvalue2, ..., eigenvalue_ith_band\
    bands = np.empty(shape=(n_k_points, n_bands))

    for ik, line in enumerate(lines):
        band_energies = [float(x) for x in line.split()[1:]]
        bands[ik, :] = band_energies

    return bands


def cif_to_ase_atoms(file: str) -> ase.atoms.Atoms:
    """Convert CIF to ASE Atoms.

    ASE's read_cif() doesn't seem robust, or I'm using it wrong, so use pymatgen instead.
    :param file:
    :return:
    """
    structure = CifParser(file).get_structures()[0]
    atoms = ase.atoms.Atoms(numbers=structure.atomic_numbers,
                            cell=structure.lattice,
                            scaled_positions=structure.frac_coords)
    return atoms
