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
from typing import List

from tb_lite.src.dftb_input import Driver, Hamiltonian, Options, DftbInput
from tb_lite.src.operations import generate_inputs

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



