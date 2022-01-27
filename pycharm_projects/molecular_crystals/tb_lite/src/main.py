"""
For each molecular crystal, vary the lattice constant from 0.8 to 1.2 * equilibrium
- Have positions in fraction
Write the TB lite input files
Run the calculation
 - note any failing
Parse the results
Plot the results
https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.cell
"""
import os
from pathlib import Path

import ase
from ase.io.dftb import write_dftb

from tb_lite.src.dftb_input import Driver, Hamiltonian, Options, generate_dftb_hsd
from tb_lite.src.tb import parse_qcore_structure, parse_qcore_settings

input_directory = 'data/entos'
materials = ['acetic.cif.in']

# TODO(Alex) Test this out for one case
# Copy the runner over


def generate_inputs(input_directory, output_directory, material):
    """

    :return:
    """
    qcore_file = os.path.join(input_directory, material)
    structure = parse_qcore_structure(qcore_file)

    driver = Driver(type='ConjugateGradient', lattice_option='No')
    ham = Hamiltonian(method='GFN1-xTB', temperature=0, scc_tolerance=1.e-6, k_points=[4, 4, 4])
    options = Options()
    dftb_input_str = generate_dftb_hsd(driver, ham, options)

    Path(output_directory).mkdir(parents=True, exist_ok=True)

    for multiplier in [0.8, 1.0]:
        # Run/output directory
        dir = output_directory + '/' + str(multiplier)
        Path(dir).mkdir(parents=True, exist_ok=True)

        # DFTB input file
        with open("dftb_in.hsd", "w") as fid:
            fid.write(dftb_input_str)

        # Output the structure:  https://wiki.fysik.dtu.dk/ase/_modules/ase/io/dftb.html#write_dftb
        # TODO Check units for Atoms positions.
        cell = []
        for vector in structure['lattice']:
            cell.append([multiplier * r for r in vector])

        atoms = ase.Atoms(symbols=structure['species'],
                          positions=structure['fractional_positions'],
                          cell=cell,
                          pbc=True)

        write_dftb("geometry.gen", atoms)




def run_inputs():
    # Run files and log run-time errors
    return None


def parse_inputs():
    # Parse files and plot
    return None




