"""
Module for converging the SCC of Bulk Systems of Interest

Approach:
* Generate an initial input for each system.
* Manually converge w.r.t. k-sampling.
* Increase things like temperature, if required.
* Tabulate the converged input here in the routine below.
"""
# TODO(Alex) Should just pick some converged value like 16, 16, 16 for every material
# and have the code also run 14, 14,14 - confirm converged.
# Can then fully automate
import os.path
from typing import Tuple
from pathlib import Path

import ase
from ase.io.dftb import write_dftb


from tb_lite.crystal_references.cubic import silicon, germanium, diamond
from tb_lite.src.dftb_input import DftbInput, Hamiltonian
from tb_lite.crystal_references import cubic, hexagonal
from tb_lite.src.parsers import cif_to_ase_atoms, parse_dftb_output
from tb_lite.src.runner import BinaryRunner


def generate_inputs():
    """ Generate some xTB inputs
    """
    root = '/Users/alexanderbuccheri/Python/pycharm_projects/tb_benchmarking/check'
    for material in ['silicon', 'germanium', 'diamond']:
        directory = os.path.join(root, material)
        Path(directory).mkdir(parents=True, exist_ok=True)
        write_an_xtb1_input(directory, material)


def write_an_xtb1_input(directory: str, material: str):
    """Given a material, write xTB input files to a directory.
    """
    atoms, input = get_material_xtb1(material)
    write_dftb(directory + "/geometry.gen", atoms)
    with open(directory + "/dftb_in.hsd", "w") as fid:
        fid.write(input.generate_dftb_hsd())


def converge_density(cif_file_name: str, calculation_dir, clear_dir=True):
    """ Converge the density for a material w.r.t. k-points
    Start with something highly-converged,
    :return:
    """
    pj = os.path.join
    results = []
    tol = 1.e-5  # eV

    if clear_dir:
        pass
        # TODO Clear directory to start

    for i, k in enumerate([12, 14, 16, 18]):
        # Geometry and input
        atoms = cif_to_ase_atoms(cif_file_name)
        input = DftbInput(
            hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[k, k, k]))

        # Write inputs
        write_dftb(pj(calculation_dir, "geometry.gen"), atoms)
        with open(pj(calculation_dir, "dftb_in.hsd"), "w") as fid:
            fid.write(input.generate_dftb_hsd())

        # Run DFTB+/TB Lite
        runner = BinaryRunner(binary='dftb+', run_cmd=[''], omp_num_threads=1, directory=calculation_dir, time_out=240)
        job_result = runner.run()

        # Parse result
        with open(pj(calculation_dir, "detailed.out"), "r") as fid:
            output_str = fid.read()
        output = parse_dftb_output(output_str)
        result = {'k_sampling': [k, k, k], 'total_energy': output['total_energy']}
        results.append(result)

        # TODO Clear directory

        # Check change in total energy
        if k > 0:
            diff = results[i]['total_energy'] - results[i-1]['total_energy']
            if abs(diff) <= tol:
                return results

        return results





def get_material_xtb1(material_name: str) -> Tuple[ase.atoms.Atoms, DftbInput]:
    """ TB lite Inputs for bulk crystals of Interest

    Manually-converged inputs:
     * Si
     * Ge
     * Diamond
    Inputs Requiring Convergence:
     * Graphite
     * Graphene
     * ZrO2
     * ZnO
     * WS2
     * GaAs
     * InAs
     * PbS
     * BN -cubic, hexagonal and wurzite
     * MoS2
     * WS2
     * TiO2 rutile
     * TiO2 anatase

    :param material_name: Material key
    :return: Atoms object and DFTB+ Input object.
    """
    # Group IV elemental crystals
    if material_name == 'silicon':
        input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[8, 8, 8]))
        return silicon(), input

    elif material_name == 'germanium':
        input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[12, 12, 12]))
        return germanium(), input

    elif material_name == 'diamond':
        input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[8, 8, 8]))
        return diamond(), input

    elif material_name == 'zinc_oxide':
        file_path = hexagonal.hexagonal_cifs.get(material_name).file
        # Should read with pymatgen and convert to ASE atoms... Doesn't seem particularly robust
        # atoms = read_cif(file_path, index=0, primitive_cell=True)
        # print(list(atoms))
        # input = DftbInput(hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[]))
        return [], []

    else:
        print(f'material_name is not valid: {material_name}')

    #elif material_name == '':
        #
    # Graphite
    # Graphene
    # ZrO2
    #
    # WS2
    # GaAs
    # InAs
    # PbS
    # BN -cubic, hexagonal and wurzite
    # MoS2
    # WS2
    # TiO2 rutile
    # TiO2 anatase
    # Consider high throughput on the set form the Sotti paper ~ 400 crystals.




file = hexagonal.hexagonal_cifs.get('zinc_oxide').file
material, atoms = cif_to_ase_atoms(file)

print(material)

