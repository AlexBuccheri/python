"""
Given some converged ground state calculation set up GW calculations
with optimised bases, for a range of energy parameter cutoffs.
"""
from pathlib import Path
import shutil
from collections import OrderedDict
from distutils.dir_util import copy_tree

from parse.lorecommendations_parser import parse_lorecommendations
from parse.parse_linengy import parse_lo_linear_energies
from parse.parse_basis_xml import parse_basis_as_string
from parse.set_gw_input import GWInput, set_gw_input_string
from process.optimised_basis import DefaultLOs, filter_lo_functions, generate_optimised_basis_string
from gw_benchmark_inputs.set3.A1_groundstate import converged_ground_state_input as A1_gs_input
from job_schedulers import slurm


def write_file(file_name, string):
    fid = open(file_name, "w")
    fid.write(string)
    fid.close()
    return

def write_input_file_with_gw_settings(root_path:str, gs_input:str, gw_input:GWInput) -> str:
    """
    Generate GW path, create a GW input file and write it

    :param root_path:
    :param gs_input:
    :param gw_input:
    :return:
    """
    q_str = ''
    for q in gw_input.ngridq:
        q_str += str(q)
    gw_root = root_path + "/gw_q" + q_str + "_omeg" + str(gw_input.nomeg) + "_nempty"+str(gw_input.nempty)

    Path(gw_root).mkdir(parents=True, exist_ok=True)

    gw_input_string = set_gw_input_string(A1_gs_input,  gw_input)
    write_file(gw_root + "/input.xml", gw_input_string)
    return gw_root

def write_optimised_lo_basis(species:str, l_max:int, energy_cutoff:dict,
                             lorecommendations, default_basis_string, default_los, job_dir):
    """
    To use the same cut-off for both species and all l-channels
      optimised_lo_cutoff = [energy_cutoff] * (l_max + 1)

    """
    optimised_lo = filter_lo_functions(lorecommendations, default_los, energy_cutoff)
    basis_string = generate_optimised_basis_string(default_basis_string, optimised_lo, max_matching_order=1)
    basis_file = species.capitalize() + '.xml'
    write_file(job_dir + '/' + basis_file, basis_string)
    return

def restructure_energy_cutoffs(n_energies:int, energy_cutoffs: dict) -> list:
    """
    Get in a more useful structure to iterate over

    TODO Automate to find n_energies

    """
    restructured_energies = []

    for inum in range(0, n_energies):
        data = {}
        for species, l_channels in  energy_cutoffs.items():
            data[species] = {l:energies[inum] for l, energies in l_channels.items()}
        restructured_energies.append(data)

    return restructured_energies


def set_up_g0w0(root_path:str):

    # Material
    species = ['zr', 'o']
    l_max = {'zr': 4, 'o': 3}

    # GW root and exciting input file
    # nempty needs to be > 1000 to account for total number of empty states
    # as a consequence of these large LO basis sets.
    # Note, i0 and i1 were ok
    gw_root = write_input_file_with_gw_settings(root_path,
                                               A1_gs_input,
                                               GWInput(taskname="g0w0", nempty=2000, ngridq=[2, 2, 2], skipgnd=False, n_omega=32)
                                               )

    # Default basis settings
    default_linear_energies = parse_lo_linear_energies(root_path + "/groundstate")
    default_los = {'zr': DefaultLOs(default_linear_energies['zr'], energy_tol=0.8),
                   'o': DefaultLOs(default_linear_energies['o'],  energy_tol=0.8)}

    # Default basis strings with .format tags
    default_basis_string = {'zr': parse_basis_as_string(root_path + "/groundstate/Zr.xml"),
                            'o': parse_basis_as_string(root_path + "/groundstate/O.xml")}

    # LO recommendation energies
    lorecommendations = parse_lorecommendations(root_path + '/lorecommendations.dat', species)

    # Optimised LO energy cutoffs. For first three calculations, keep Zr l=2 fixed and increase all other channels ALOT.
    # See what effect that has. Then increase zr l=2 more to check that it is also converged.
    # Final calculation. Keep all channels the same as the run, but increase lmax=2 to check it's converged
    n_energies_per_channel = 8
    #                             i0    i1   i2   i3   i4   i5   i6   i7      Struggling to run i3, but i1 looks converged, hence why I've added i4
    energy_cutoffs =  {'zr': {0: [150, 200, 250, 250, 200, 120, 120, 120],
                              1: [150, 200, 250, 250, 200, 120, 120, 120],
                              2: [300, 300, 300, 350, 350, 350, 400, 450],
                              3: [150, 200, 250, 250, 200, 120, 120, 120],
                              4: [150, 200, 250, 250, 200, 120, 120, 120]},

                       'o':  {0: [150, 200, 250, 250, 200, 120, 120, 120],
                              1: [150, 200, 250, 250, 200, 120, 120, 120],
                              2: [150, 200, 250, 250, 200, 120, 120, 120],
                              3: [150, 200, 250, 250, 200, 120, 120, 120]}
                       }

    # Note, i6 misses a function at 401 Ha and whilst i7 runs, it returns a metallic solution
    # Re-running i7 with

    # Slurm script settings
    env_vars = OrderedDict([('EXE', '/users/sol/abuccheri/exciting/bin/excitingmpismp'),
                            ('OUT', 'terminal.out')
                            ])
    module_envs = ['intel/2019']
    slurm_directives = slurm.set_slurm_directives(time=[0, 24, 0, 0],
                                                  partition='all',
                                                  exclusive=True,
                                                  nodes=4,
                                                  ntasks_per_node=2,
                                                  cpus_per_task=18,
                                                  hint='nomultithread')

    species_basis_string = ''
    for s in species:
        species_basis_string += s.capitalize() + str(l_max[s]) + '_'

    print('here')

    for ie, energy_cutoff in enumerate(restructure_energy_cutoffs(n_energies_per_channel, energy_cutoffs)):
        # Copy groundstate directory to GW directory
        # Use an index not max energy, as the max energy does not change in 3/4 runs
        job_dir = gw_root + '/max_energy_i' + str(ie)
        print('Creating directory, with input.xml, run.sh and optimised basis:', job_dir)
        print(root_path)
        copy_tree(root_path + '/groundstate', job_dir)

        # Copy input.xml with GW settings
        shutil.copy(gw_root + "/input.xml", job_dir + "/input.xml")

        # New Slurm script
        slurm_directives['job-name'] = "gw_A1_lmax_" + species_basis_string + str(ie) + 'loEcutoff'
        write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))

        # Write optimised basis
        write_optimised_lo_basis('zr', l_max['zr'], energy_cutoff['zr'], lorecommendations['zr'],
                                 default_basis_string['zr'], default_los['zr'], job_dir)
        write_optimised_lo_basis('o', l_max['o'], energy_cutoff['o'], lorecommendations['o'],
                                 default_basis_string['o'], default_los['o'], job_dir)

    return


set_up_g0w0("/users/sol/abuccheri/gw_benchmarks/A1_set3/zr_lmax4_o_lmax3_rgkmax7")
