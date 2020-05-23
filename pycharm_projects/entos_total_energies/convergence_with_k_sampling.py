# Converge total energy w.r.t. number of k-points
# These are important checks to do but not likely to be app tests

import numpy as np
import json
import subprocess
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import time
import typing

# My modules
from options import *


entos_exe = '/Users/alexanderbuccheri/Codes/entos/cmake-build-debug/entos'


# ------------------------
# Run and plot functions
# ------------------------

def plot_convergence(n_k_points:typing.List[int],
                     n_irreducible_k_points:typing.List[int],
                     total_energy:typing.List[float],
                     title:str=None):
    """Plot convergence in total energy w.r.t. MP grid density
    """
    fig, ax = plt.subplots()
    ax.plot(n_k_points, total_energy, color='blue', marker='o', linestyle='dashed', linewidth=2, markersize=12)

    # Add number of irreducible k-points as numbers next to each point
    for i, n_irreducible_k_point in enumerate(n_irreducible_k_points):
        ax.annotate(n_irreducible_k_point, (n_k_points[i] + 0.02, total_energy[i] + 0.02), size=14, color='red')

    # See: https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/custom_legends.html
    red_line = mlines.Line2D([0], [0], color='w', markeredgecolor='red', marker='$1,2,...$', markersize=20,
                              label='Num of irreducible k-points')
    plt.legend(handles=[red_line])

    if title != None:
        plt.title(title, fontsize=16)
    plt.xlabel("Number of Monkhorst-Pack k-points", fontsize=16)
    plt.ylabel("Total Energy (H)", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    return fig, ax


def run_convergence(named_result:str,
                    k_grids:typing.List,
                    xtb_string_function:typing.Callable,
                    symmetry_reduction=False):
    """Converge total energy w.r.t. MP grid density
    """
    n_k_points = []
    n_irreducible_k_points = []
    total_energy = []
    timing = []

    for k, k_grid in enumerate(k_grids):
        n_k_points.append(np.prod(k_grid))
        named_result += named_result + str(n_k_points[k])

        kgrid_options = KGridOptions(grid_integers=k_grid, symmetry_reduction=symmetry_reduction)
        input_string = xtb_string_function(named_result, kgrid_options)

        entos_command = [entos_exe, '--format', 'json', '-s', input_string]
        t1 = time.perf_counter()
        entos_json_result = subprocess.check_output(entos_command)
        t2 = time.perf_counter()
        timing.append(t2 - t1)
        result = json.loads(entos_json_result)

        if print_level > 0:
            print("xTB input script: ", input_string)
            print("Result keys: ", result[named_result].keys())

        total_energy.append(result[named_result]['energy'])
        n_irreducible_k_points.append(result[named_result]['n_k_points'])

    print("# N grid points, N irreducible grid points, total energy (H), time (s)")
    for k in range(0, len(k_grids)):
        print(n_k_points[k], n_irreducible_k_points[k], total_energy[k], timing[k])

    return (n_k_points, n_irreducible_k_points, total_energy)


# -----------------------------
# Material settings
# -----------------------------

def primitive_neutral_fcc_input_string(named_result:str, kgrid_options, element, a):

    assert isinstance(kgrid_options, KGridOptions)
    atoms = [["'" + element + "'", 0, 0, 0],
             ["'" + element + "'", 0.25, 0.25, 0.25]]

    lattice_options = LatticeOpt(a=a, length_unit='angstrom')
    structure_string = generate_structure_string(unit='fractional', lattice=lattice_options, atoms=atoms, bravais='fcc')
    cutoff_options = TranslationCutoffOptions(h0=40, overlap=40, repulsive=40)
    ewald_options = EwaldOptions(real_cutoff=10, recip_cutoff=2, alpha=0.5)
    temperature_option = SingleOption(InputEntry(command="temperature", value=0, unit="kelvin"))

    options_string = options_to_string(cutoff_options, ewald_options, kgrid_options, temperature_option)
    xtb_string = generate_xtb_string(named_result, structure_string, options_string)
    xtb_string.replace('\n', ' ')

    return xtb_string


def primitive_silicon_input_string(named_result:str, kgrid_options):
    """ Set up an entos xTB input file for primitive silicon
        Ref for lattice constant: Probably wiki.
    """
    element = 'Si'
    a = 5.431
    return primitive_neutral_fcc_input_string(named_result, kgrid_options, element=element, a=a)


def primitive_diamond_input_string(named_result:str, kgrid_options):
    """ Set up an entos xTB input file for primitive diamond.
        Ref for experimental lattice constant: https://doi.org/10.1103/PhysRevB.24.6121
    """
    element = 'C'
    a = 3.567
    return primitive_neutral_fcc_input_string(named_result, kgrid_options, element=element, a=a)


def primitive_germanium_input_string(named_result:str, kgrid_options):
    """ Set up an entos xTB input file for primitive diamond.
        Ref for experimental lattice constant: https://doi.org/10.1103/PhysRevB.24.6121
    """
    assert isinstance(kgrid_options, KGridOptions)
    atoms = [['Ge', 0, 0, 0],
             ['Ge', 0.25, 0.25, 0.25]]

    lattice_options = LatticeOpt(a=5.652, length_unit='angstrom')
    structure_string = generate_structure_string(unit='fractional', lattice=lattice_options, atoms=atoms, bravais='fcc')
    cutoff_options = TranslationCutoffOptions(h0=40, overlap=40, repulsive=40)
    ewald_options = EwaldOptions(real_cutoff=10, recip_cutoff=2, alpha=0.5)
    temperature_option = SingleOption(InputEntry(command="temperature", value=0, unit="kelvin"))
    solver_option = SingleOption(InputEntry(command="solver", value='non_iterative'))

    options_string = options_to_string(cutoff_options, ewald_options, kgrid_options, temperature_option, solver_option)
    xtb_string = generate_xtb_string(named_result, structure_string, options_string)
    xtb_string.replace('\n', ' ')

    return xtb_string


# --------------------------------
# Specific Material Calculations
# --------------------------------

# Could generalise this function but it doesn't save much typing
def run_silicon_convergence(show_plot=True):
    k_grids = [[1, 1, 1], [1, 2, 1], [1, 2, 2], [2,2,2], [3,3,3], [4,4,4], [5,5,5], [6,6,6], [7,7,7], [8,8,8]]
    (n_k_points, n_irreducible_k_points, total_energy) = \
        run_convergence('silicon', k_grids, primitive_silicon_input_string, symmetry_reduction=True)
    if show_plot:
        plot_title = "Convergence with k-sampling. Bulk silicon"
        fig, ax = plot_convergence(n_k_points, n_irreducible_k_points, total_energy, title=plot_title)
        plt.show()
    return

def run_diamond_convergence(show_plot=True):
    k_grids = [[1, 1, 1], [1, 2, 1], [1, 2, 2], [2,2,2], [3,3,3], [4,4,4], [5,5,5], [6,6,6]]
    (n_k_points, n_irreducible_k_points, total_energy) = \
        run_convergence('diamond', k_grids, primitive_diamond_input_string, symmetry_reduction=True)
    if show_plot:
        plot_title = "Convergence with k-sampling. Bulk diamond"
        fig, ax = plot_convergence(n_k_points, n_irreducible_k_points, total_energy, title=plot_title)
        plt.show()
    return

def run_germanium_convergence(show_plot=True):
    k_grids = [[1, 1, 1], [1, 2, 1], [1, 2, 2], [2,2,2], [3,3,3], [4,4,4], [5,5,5], [6,6,6]]
    (n_k_points, n_irreducible_k_points, total_energy) = \
        run_convergence('germanium', k_grids, primitive_germanium_input_string, symmetry_reduction=True)

    print("\n Script Note: Germanium fails for self-consistent calculations. "
          "Same behaviour witnessed in CP2K (August 2019). As such, run non-SCC")

    if show_plot:
        plot_title = "Convergence with k-sampling. Bulk Ge"
        fig, ax = plot_convergence(n_k_points, n_irreducible_k_points, total_energy, title=plot_title)
        plt.show()
    return


# --------------------------------
# Run specific Material Calculations
# --------------------------------
print_level = 1

#run_silicon_convergence()
#run_diamond_convergence()
#run_germanium_convergence()






