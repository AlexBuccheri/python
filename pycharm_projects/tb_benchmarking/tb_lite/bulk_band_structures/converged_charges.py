"""Module for converging the SCC of Bulk Systems of Interest
"""
import os.path
from typing import List, Union
from pathlib import Path
import json

import ase
from ase.io.dftb import write_dftb

from tb_lite.crystal_references import cubic, hexagonal, tetragonal, monoclinic
from tb_lite.src.dftb_input import DftbInput, Hamiltonian
from tb_lite.src.runner import BinaryRunner
from tb_lite.src.parsers import cif_to_ase_atoms, parse_dftb_output, clear_directory


def converge_densities(calculations: dict) -> dict:
    """ Automatically converge SCC for a set of calculations.

    :return: Dictionary of all results
    """
    all_results = {}

    for name, settings in calculations.items():
        print(f"Running calculation {name}")

        if set(settings) != {'atoms', 'directory'}:
            raise KeyError(f"Calculation is missing all required settings.\n"
                           f"Requires (atoms, directory). Contains {set(settings)}")

        all_results[name] = {'directory': settings['directory'].as_posix(),
                             'energy_vs_k': converge_density(settings['atoms'],
                                                             settings['directory'],
                                                             clear_dir=True)
                             }
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

    k_grids = [4, 6, 8, 10, 12, 14, 16, 18, 20]

    for i, k in enumerate(k_grids):
        # Clear any existing files
        clear_directory(calculation_dir)

        # Input file
        input = DftbInput(
            hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=[k, k, k]))

        # Write inputs
        write_dftb(pj(calculation_dir, "geometry.gen"), atoms)
        with open(pj(calculation_dir, "dftb_in.hsd"), "w") as fid:
            fid.write(input.generate_dftb_hsd())

        # Run DFTB+/TB Lite
        runner = BinaryRunner(binary='dftb+', run_cmd=['./'], omp_num_threads=1, directory=calculation_dir, time_out=600)
        job_result = runner.run()

        # Parse result
        if job_result.success:
            with open(pj(calculation_dir, "detailed.out"), "r") as fid:
                output_str = fid.read()
            output = parse_dftb_output(output_str)
            result = {'k_sampling': [k, k, k], 'total_energy': output['total_energy']}
        else:
            print('Calculation failed')
            result = {'k_sampling': [k, k, k], 'total_energy': 'NULL'}

        results.append(result)

        # Check change in total energy
        if i > 0:
            #print(i, [k, k, k], results[i]['total_energy'], results[i-1]['total_energy'])
            diff = results[i]['total_energy'] - results[i - 1]['total_energy']

            if abs(diff) <= tol:
                # Re-run the prior calculation so charges are available
                converged_kgrid = [k_grids[i-1]] * 3
                input = DftbInput(
                    hamiltonian=Hamiltonian(method='GFN1-xTB', temperature=0.0, scc_tolerance=1.e-6, k_grid=converged_kgrid))
                write_dftb(pj(calculation_dir, "geometry.gen"), atoms)
                with open(pj(calculation_dir, "dftb_in.hsd"), "w") as fid:
                    fid.write(input.generate_dftb_hsd())
                runner = BinaryRunner(binary='dftb+', run_cmd=['./'], omp_num_threads=1, directory=calculation_dir, time_out=600)
                job_result = runner.run()
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
    #
    # To Add:  GaP, InN, InAs, PbSe, graphene
    # Added: Si, Ge, Diamond, MoS2, WS2, ZnO, BN hexagonal, BN cubic, MgO, copper, sodium_chloride,
    # ZrO2, PbS, TiO2 anatase, TiO2 rutile
    #
    # Note, would be easier to just query the Materials Project API for a list of bulk materials
    # with the correct symmetry, and take the first hit in each case.
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
                 'pbs':        cif_to_ase_atoms(cubic.fcc_cifs.get('lead_sulfide').file),
                 'tio2_rutile': tetragonal.tio2_rutile(),
                 'tio2_ana':    tetragonal.tio2_anatase(),
                 'cdse':       cif_to_ase_atoms(hexagonal.hexagonal_cifs.get('cadmium_selenide').file),
                 'gan':        cif_to_ase_atoms(cubic.fcc_cifs.get('gallium_nitride').file),
                 'graphite':   cif_to_ase_atoms(hexagonal.hexagonal_cifs.get('graphite').file),
                 'gaas':       cif_to_ase_atoms(cubic.fcc_cifs.get('gallium_arsenide').file),
                 'wo3_monoclinic': cif_to_ase_atoms(monoclinic.simple_monoclinic_cifs.get("tungsten_oxide").file),
      #           'inp':        cif_to_ase_atoms(cubic.fcc_cifs.get("indium_phosphide").file),
                 'pbte':       cif_to_ase_atoms(cubic.fcc_cifs.get("lead_telluride").file)
                 }

    # Calculations
    root = Path("/home/alex/tblite/bulk/scc")
    calculations = {}
    for name, atoms in materials.items():
        calculations[name]= {'atoms': atoms, 'directory': define_directory(root, name)}

    for name, settings in calculations.items():
        Path.mkdir(settings['directory'], exist_ok=True)

    all_results = converge_densities(calculations)

    print("Dump all results\n")
    print(all_results)

    with open('data.json', 'w', encoding='utf-8') as fid:
        json.dump(all_results, fid, ensure_ascii=False, indent=4)

