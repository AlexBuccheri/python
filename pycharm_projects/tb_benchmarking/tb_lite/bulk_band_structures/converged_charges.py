"""Module for converging the SCC of Bulk Systems of Interest
"""
import sys
import os.path
from typing import List, Union
from pathlib import Path
import json

import ase
from ase.io.dftb import write_dftb

from tb_lite.crystal_references import cubic, hexagonal, tetragonal
from tb_lite.src.dftb_input import DftbInput, Hamiltonian
from tb_lite.src.runner import BinaryRunner
from tb_lite.src.parsers import cif_to_ase_atoms, parse_dftb_output, clear_directory


def converge_densities(calculations: dict) -> dict:
    """ Automatically converge SCC for a set of calculations.

    :return: Dictionary of all results
    """
    all_results = calculations.copy()

    for name, settings in calculations.items():
        print(f"Running calculation {name}")

        if set(settings) != {'atoms', 'directory'}:
            raise KeyError(f"Calculation is missing all required settings.\n"
                           f"Requires (atoms, directory). Contains {set(settings)}")

        all_results[name]['energy_vs_k'] = converge_density(settings['cif'],
                                                            settings['directory'],
                                                            clear_dir=True)
        print(all_results[name]['energy_vs_k'])

    return all_results


def converge_density(atoms: ase.atoms.Atoms, calculation_dir, clear_dir=True) -> List[dict]:
    """ Converge the density for a material w.r.t. k-points.

    Note, some materials will converge very quickly w.r.t. k. Others may require
    the full range of sampling to be tested.

    :return: List of k_sampling and total_energy, per calculation.
    """
    pj = os.path.join
    results = []
    # Tolerance for convergence in total energy, in eV
    tol = 1.e-5

    if clear_dir:
        clear_directory(calculation_dir)

    for i, k in enumerate([4, 6, 8, 10, 12, 14, 16, 18, 20]):
        # Input file
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
        if job_result.success:
            with open(pj(calculation_dir, "detailed.out"), "r") as fid:
                output_str = fid.read()
            output = parse_dftb_output(output_str)
            result = {'k_sampling': [k, k, k], 'total_energy': output['total_energy']}
        else:
            result = {'k_sampling': [k, k, k], 'total_energy': 'NULL'}

        results.append(result)
        clear_directory(calculation_dir)

        # Check change in total energy
        if k > 0:
            diff = results[i]['total_energy'] - results[i - 1]['total_energy']
            if abs(diff) <= tol:
                return results

        return results


def define_directory(root: Path, name: Union[str, Path]) -> Path:
    """ Define directories for calculation to run in.

    :return: Dictionary of directories settings
    """
    return root / Path(name)


if __name__ == "__main__":
    # Converge SCC calculations using DFTB+ TB Lite
    # Materials dict could be simplified if every system was defined with a CIF file

    # To Add: WO3, GaAs, InP  PbTe, CdSe, graphene
    # To Add: GaN, GaP, InN, InAs, PbSe
    # Added: Si, Ge, Diamond, MoS2, WS2, ZnO, BN hexagonal, BN cubic, MgO, copper, sodium_chloride,
    # ZrO2, PbS, TiO2 anatase, TiO2 rutile

    major, minor = sys.version_info[0:2]
    if major < 3 or minor < 6:
        sys.exit("Script relies on consistent dictionary ordering")

    # Materials {key: value} = ' ase Atoms object
    materials = {'silicon':    cubic.silicon(),
                 'germanium':  cubic.germanium(),
                 'diamond':    cubic.diamond(),
                 'zinc_oxide': cif_to_ase_atoms(hexagonal.hexagonal_cifs.get('zinc_oxide').file),
                 'mos2':       cif_to_ase_atoms(hexagonal.hexagonal_cifs.get('molybdenum_disulfide').file),
                 'ws2':        cif_to_ase_atoms(hexagonal.hexagonal_cifs.get('tungsten_disulfide').file),
                 'bn_hex':     hexagonal.boron_nitride_hexagonal(),
                 'bn_cubic':   cif_to_ase_atoms(cubic.fcc_cifs.get('boron_nitride').file),
                 'mgo':        cif_to_ase_atoms(cubic.fcc_cifs.get('magnesium_oxide').file),
                 'copper':     cif_to_ase_atoms(cubic.fcc_cifs.get('copper').file),
                 'nacl':       cif_to_ase_atoms(cubic.fcc_cifs.get('sodium_chloride').file),
                 'zro2':       cif_to_ase_atoms(cubic.fcc_cifs.get('zirconium_dioxide').file),
                 'pbs':        cif_to_ase_atoms(cubic.fcc_cifs.get('pbs').file),
                 'tio2_rutile': tetragonal.tio2_rutile(),
                 'tio2_ana':    tetragonal.tio2_anatase()
                 }

    # Calculations
    root = Path("/home/alex/tblite/bulk/scc")
    calculations = {}
    for name, atoms in materials.items():
        calculations = {name: {'atoms': atoms, 'directory': define_directory(root, name)}}

    all_results = converge_densities(calculations)

    with open('data.json', 'w', encoding='utf-8') as fid:
        json.dump(all_results, fid, ensure_ascii=False, indent=4)
