"""
TODOs
Run the calculation
 - note any failing
Parse the results
Plot the results
"""
import os.path
from typing import List
import numpy as np

from tb_lite.src.dftb_input import Driver, Hamiltonian, Options, DftbInput
from tb_lite.src.operations import generate_inputs
from tb_lite.src.runner import BinaryRunner, SubprocessRunResults

input_directory = 'data/entos'

all_materials = [
    'urea_g.cif.in',
    'urea.cif.in',
    'uracil_g.cif.in',
    'uracil.cif.in',
    'trioxane_g.cif.in',
    'trioxane.cif.in',
    'triazine_g.cif.in',
    'triazine.cif.in',
    'succinic_g.cif.in',
    'succinic.cif.in',
    'pyrazole_g.cif.in',
    'pyrazole.cif.in',
    'pyrazine_g.cif.in',
    'pyrazine.cif.in',
    'oxacb_g.cif.in',
    'oxacb.cif.in',
    'oxaca_g.cif.in',
    'oxaca.cif.in',
    'naph_g.cif.in',
    'naph.cif.in',
    'imdazole_g.cif.in',
    'imdazole.cif.in',
    'hexdio_g.cif.in',
    'hexdio.cif.in',
    'hexamine_g.cif.in',
    'hexamine.cif.in',
    'formamide_g.cif.in',
    'formamide.cif.in',
    'ethcar_g.cif.in',
    'ethcar.cif.in',
    'cytosine_g.cif.in',
    'cytosine.cif.in',
    'cyanamide_g.cif.in',
    'cyanamide.cif.in',
    'benzene_g.cif.in',
    'benzene.cif.in',
    'anthracene_g.cif.in',
    'anthracene.cif.in',
    'ammonia_g.cif.in',
    'ammonia.cif.in',
    'adaman_g.cif.in',
    'adaman.cif.in',
    'acetic_g.cif.in',
    'acetic.cif.in',
    'CO2_g.cif.in',
    'CO2.cif.in']


def inputs(materials: List[str]):
    for material in materials:
        dftb_input = DftbInput(
            hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0, scc_tolerance=1.e-6, k_points=[4, 4, 4]),
            driver=Driver(type='ConjugateGradient', lattice_option='No'), options=Options())
        generate_inputs(input_directory, 'outputs', material, dftb_input)


def run_jobs(root):

    materials = ['CO2.cif.in']
    lattice_multipliers = [0.8, 1.0]
    run_command = ['./']
    omp = 2
    time_out = 300

    for material in materials:
        for multiplier in lattice_multipliers:
            run_dir = os.path.join(root, 'outputs', material, str(multiplier))
            print('Executing {run_dir}')
            runner = BinaryRunner('dftb+', run_command, omp, time_out, run_dir)
            process_results = runner.run()
            print(process_results.success)

run_jobs()