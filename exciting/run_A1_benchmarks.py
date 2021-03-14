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
from gw_benchmark_inputs.A1_groundstate import converged_ground_state_input as A1_gs_input
from job_schedulers import slurm


def write_file(file_name, string):
    fid = open(file_name, "w")
    fid.write(string)
    fid.close()
    return

def write_input_file_with_gw_settings(gw_root, gs_input, gw_input):
    gw_input_string = set_gw_input_string(A1_gs_input,  gw_input)
    write_file(gw_root + "/input.xml", gw_input_string)
    return

def write_optimised_lo_basis(species:str, l_max:int, energy_cutoff:float,
                             lorecommendations, default_basis_string, default_los):
    # Cut-off same for both species and all l-channels
    optimised_lo_cutoff = [energy_cutoff] * (l_max + 1)
    optimised_lo = filter_lo_functions(lorecommendations, default_los, optimised_lo_cutoff)
    basis_string = generate_optimised_basis_string(default_basis_string, optimised_lo, max_matching_order=1)
    basis_file = species.capitalize() + '.xml'
    write_file(job_dir + '/' + basis_file, basis_string)
    return


# Produce one of these files per top-level directory, such that all the settings are preserved
def set_up_g0w0(root_path:str):

    # Material
    species = ['zr', 'o']
    l_max = {'zr': 4, 'o': 3}

    # GW root
    gw_root = root_path + "/gw_q222_omeg32_nempty500"
    Path(gw_root).mkdir(parents=True, exist_ok=True)

    # exciting input file
    write_input_file_with_gw_settings(gw_root, A1_gs_input,
                                     GWInput(taskname="g0w0", nempty=500, ngridq=[2, 2, 2], skipgnd=False))

    # Default basis settings
    default_linear_energies = parse_lo_linear_energies(root_path + "/groundstate")
    default_los = {'zr': DefaultLOs(default_linear_energies['zr'], energy_tol=0.1),
                   'o': DefaultLOs(default_linear_energies['o'],  energy_tol=0.1)}

    # Default basis strings with tags
    default_basis_string = {'zr': parse_basis_as_string(root_path + "/groundstate/Zr.xml"),
                            'o': parse_basis_as_string(root_path + "/groundstate/O.xml")}

    # LO recommendation energies
    lorecommendations = parse_lorecommendations(root_path + '/lorecommendations.dat', species)

    # Optimised LO energy cutoffs
    energy_cutoffs = [100]   #[20, 40, 60, 80, 100, 120, 140, 160]


    # Slurm script settings
    env_vars = OrderedDict([('EXE', '/users/sol/abuccheri/exciting/bin/excitingmpismp'),
                            ('OUT', 'terminal.out')
                            ])
    module_envs = ['intel/2019']
    slurm_directives = slurm.set_slurm_directives(time=[0, 6, 0, 0],
                                                  partition='all',
                                                  exclusive=True,
                                                  nodes=1,
                                                  ntasks_per_node=2,
                                                  cpus_per_task=18,
                                                  hint='nomultithread')

    species_basis_string = [s.capitalize() + str(l_max[s]) + '_' for s in species]

    for energy_cutoff in energy_cutoffs:
        # Copy groundstate directory to GW directory
        job_dir = gw_root + '/max_energy_' + str(energy_cutoff)
        print('Creating directory, with input.xml, run.sh and optimised basis:', job_dir)
        copy_tree(root_path +'/groundstate', job_dir)

        # Copy input.xml with GW settings
        shutil.copy(gw_root + "/input.xml", job_dir + "/input.xml")

        # New Slurm script
        slurm_directives['job-name'] = "gw_A1_lmax_" + species_basis_string + str(energy_cutoff) +'loEcutoff'
        write_file(job_dir + '/run.sh', slurm.set_slurm_script(slurm_directives, env_vars, module_envs))

        # Write optimised basis
        write_optimised_lo_basis('zr', l_max['zr'], energy_cutoff, lorecommendations['zr'], default_basis_string['zr'], default_los['zr'])
        write_optimised_lo_basis('o', l_max['o'], energy_cutoff, lorecommendations['o'], default_basis_string['o'], default_los['o'])

    return


set_up_g0w0("/users/sol/abuccheri/second_test")