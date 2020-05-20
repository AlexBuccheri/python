import numpy as np
import json
import subprocess
import matplotlib.pyplot as plt

from options import *


# Set up an entos xTB input file for primitive silicon
def primitive_silicon(named_result, kgrid_options):
    atoms = [['Si', 0, 0, 0],
             ['Si', 0.25, 0.25, 0.25]]

    lattice_options = LatticeOpt(a=5.431, length_unit='angstrom')
    structure_string = generate_structure_string(unit='fractional', lattice=lattice_options, atoms=atoms, bravais='fcc')
    cutoff_options = TranslationCutoffOptions(h0=40, overlap=40, repulsive=40)
    ewald_options = EwaldOptions(real_cutoff=10, recip_cutoff=2, alpha=0.5)
    temperature_option = SingleOption(InputEntry(command="temperature", value=0, unit="kelvin"))

    options_string = options_to_string(cutoff_options, ewald_options, kgrid_options, temperature_option)
    xtb_string = generate_xtb_string(named_result, structure_string, options_string)
    xtb_string.replace('\n', ' ')

    return xtb_string


# Run a job
print_level = 1
entos_exe = '/Users/alexanderbuccheri/Codes/entos/cmake-build-debug/entos'
k_grids = [[1,1,1], [2,2,2], [3,3,3], [4,4,4], [5,5,5], [6,6,6], [7,7,7], [8,8,8]]

# Converge total energy w.r.t. MP grid density
n_k_points = []
n_irreducible_k_points = []
total_energy = []

for k, k_grid in enumerate(k_grids):
    n_k_points.append(np.prod(k_grid))
    named_result = "si_kgrid_" + str(n_k_points[k])
    kgrid_options = KGridOptions(grid_integers=k_grid, symmetry_reduction=True)
    silicon_input_string = primitive_silicon(named_result, kgrid_options)

    entos_command = [entos_exe, '--format', 'json', '-s', silicon_input_string]
    entos_json_result = subprocess.check_output(entos_command)
    result = json.loads(entos_json_result)

    if print_level > 0:
        print("xTB input script: ", silicon_input_string)
        print("Result keys: ", result[named_result].keys())

    total_energy.append(result[named_result]['energy'])
    n_irreducible_k_points.append(result[named_result]['n_k_points'])


# Want to time it - debug verses release

print("# N grid points, N irreducible grid points, total energy (H)")
for k in range(0, len(k_grids)):
    print(n_k_points[k], n_irreducible_k_points[k], total_energy[k])


