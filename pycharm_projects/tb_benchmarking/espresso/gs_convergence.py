""" Converge GS of bulk materials with QE.
"""
from espresso.inputs import set_espresso_input
from espresso.scf import run_scf_calculation
from tb_lite.crystal_references.crystal_systems import bulk_materials
from tb_lite.src.parse_convergence import parse_converged_kgrid

if __name__ == "__main__":

    # Settings from TB Lite
    json_file = "tb_lite/bulk_band_structures/converged_energies.json"
    # Review this choice
    scf_tolerance = 1.e-6
    converged_k_grids, unconverged_k_grids = parse_converged_kgrid(json_file, scf_tolerance)

    # Espresso env Settings
    env = {'ASE_ESPRESSO_COMMAND': "/users/sol/abuccheri/packages/qe-7.1/build/bin/pw.x -in PREFIX.pwi > PREFIX.pwo",
           'ESPRESSO_PSEUDO': ''}
    pseudo_dir = '/users/sol/abuccheri/rutgers_pseudos/pbesol'
    run_dir = '/users/sol/abuccheri/espresso_tests/workflow_test'

    material = 'silicon'
    atoms = bulk_materials[material]
    k_grid = converged_k_grids[material]

    for ecut in [50.]: #, 60., 70., 80., 90., 100.]:
        qe_input = set_espresso_input(material, ecutwfc=ecut)
        atoms, calculator = run_scf_calculation(run_dir, qe_input, atoms, k_grid)
        print(calculator.get_fermi_level())
