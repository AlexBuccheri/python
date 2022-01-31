"""
"""
import os.path
from typing import List
import numpy as np
import matplotlib.pyplot as plt

from tb_lite.src.dftb_input import Hamiltonian, Options, DftbInput
from tb_lite.src.operations import generate_inputs, directory_name
from tb_lite.src.runner import BinaryRunner
from tb_lite.src.tb_parsing import parse_tb_output


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
    """ Generate file inputs

    :param materials:
    :return:
    """
    input_directory = 'data/entos'

    for material in materials:
        dftb_input = DftbInput(
            hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0, scc_tolerance=1.e-6, k_points=[4, 4, 4]),
            options=Options())
        generate_inputs(input_directory, 'outputs', material, dftb_input)


def run_jobs(root):
    """ Run jobs

    :param root:
    :return:
    """
    materials = all_materials
    lattice_multipliers = np.arange(0.8, 1.2, step=0.025)
    omp = 2
    time_out = 120

    for material in materials:
        for multiplier in lattice_multipliers:
            run_dir = directory_name(root + '/outputs', material, multiplier)
            print(f'Executing {run_dir}')
            runner = BinaryRunner("dftb+", ['./'], omp, time_out, run_dir)
            process_results = runner.run()
            print(run_dir, process_results.success)


def get_results(root):
    """ Get results and plot energy vs volume

    :param root:
    :return:
    """
    materials = ['trioxane_g']
    lattice_multipliers = np.arange(0.8, 1.2, step=0.025)
    output_directory = os.path.join(root, 'outputs')

    # Collated results
    total_energy = np.empty(shape=(len(lattice_multipliers)))
    total_energies = {}

    for material in materials:

        for i, multiplier in enumerate(lattice_multipliers):
            dir = directory_name(output_directory, material, multiplier)
            with open(os.path.join(dir, "detailed.out")) as fid:
                file_str = fid.read()
            results = parse_tb_output(file_str)
            total_energy[i] = results['total_energy']

        total_energies[material] = total_energy

    print(total_energies)
    plt.plot(lattice_multipliers, total_energies['trioxane_g'], 'go--', linewidth=2, markersize=12)
    plt.show()



inputs(all_materials)
#run_jobs('/home/alex/tblite/molecular_crystals')
