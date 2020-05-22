# Translate the atoms in the primitive cell
# by a fixed amount and confirm that the total energy is unchanged
# Note, arbitrary shift should be small enough that all atoms stay within the unit cell

import json
import subprocess
import collections


entos_exe = '/Users/alexanderbuccheri/Codes/entos/cmake-build-debug/entos'

# Run entos energy differences from python
run_energy_differences = False

# Generate equivalent app test inputs from python
generate_app_test_inputs = True


# -----------------------------------
# Functions
# -----------------------------------

def list_to_string(mylist: list, joiner: str) -> str:
    return joiner.join([str(i) for i in mylist])

def get_energies(named_result, inputs):
    energy = []
    for input in inputs:
        entos_command = [entos_exe, '--format', 'json', '-s', input.replace('\n', ' ')]
        entos_output = subprocess.check_output(entos_command)
        result = json.loads(entos_output)
        energy.append(result[named_result]['energy'])
    return energy

def get_energy_differences(named_result, strings_function, arbitrary_shift):
    input_strings = strings_function(named_result, arbitrary_shift)
    energies = get_energies(named_result, input_strings)

    energy_diffs = []
    for i in range(1, len(energies)):
        energy_diffs.append(abs(energies[i] - energies[i-1]))

    return energy_diffs

def evaluate_energy_differences(energy_difference, tol=1.e-6):
    all_pass = True
    for key, value in energy_difference.items():
        if value > tol:
            print(key + ' primitive unit cell')
            print("Energy should be invariant w.r.t. translation of basis positions,"
                  "however we find |delta E| (Ha) = ", value)
            all_pass = False
    return all_pass

# ----------------------------------------
# entos string functions for translations
# ----------------------------------------

def silicon_translational_invariance_strings(named_result, arbitrary_shift):
    input_string = named_result + """ := xtb(
        structure( fractional=[[Si,{pos0:s}],
                               [Si,{pos1:s}]] 
           lattice(a = 5.431 angstrom
                   bravais = fcc
                   )
                 )
          h0_cutoff = 40 bohr
          overlap_cutoff = 40 bohr
          repulsive_cutoff = 40 bohr
          ewald_real_cutoff = 10 bohr
          ewald_reciprocal_cutoff = 2 
          ewald_alpha = 0.5 
          monkhorst_pack = [2,2,2] 
          symmetry_reduction = true 
          temperature = 0 kelvin
        )
        """

    fractional_0 = [0.00, 0.00, 0.00]
    fractional_1 = [0.25, 0.25, 0.25]

    fractional_0_shifted = [round(i + arbitrary_shift, 3) for i in fractional_0]
    fractional_1_shifted = [round(i + arbitrary_shift, 3) for i in fractional_1]

    fractional_0_str = list_to_string(fractional_0, ',')
    fractional_1_str = list_to_string(fractional_1, ',')
    input_noshift = input_string.format(pos0=fractional_0_str, pos1=fractional_1_str)

    fractional_0_shifted_str = list_to_string(fractional_0_shifted, ',')
    fractional_1_shifted_str = list_to_string(fractional_1_shifted, ',')
    input_shift = input_string.format(pos0=fractional_0_shifted_str, pos1=fractional_1_shifted_str)
    return (input_noshift, input_shift)


def rutile_translational_invariance_strings(named_result, arb_shift):
    # Reference: https://materialsproject.org/materials/mp-2657/
    input_string = named_result + """ := xtb(
                       structure( fractional=[[Ti, {pos0:s}],
                                              [Ti, {pos1:s}],
                                              [ O, {pos2:s}],
                                              [ O, {pos3:s}],
                                              [ O, {pos4:s}],
                                              [ O, {pos5:s}]]
                              lattice( a = 4.6068 angstrom
               	                    c = 2.9916 angstrom
                                       bravais = tetragonal
                                      )
                                )

                        repulsive_cutoff = 40.0 bohr
                        overlap_cutoff = 40.0 bohr
                        h0_cutoff = 40.0 bohr
                        temperature = 0 kelvin
                        ewald_real_cutoff = 10 bohr
                        ewald_reciprocal_cutoff = 2
                        ewald_alpha = 0.5
                        monkhorst_pack = [2, 2, 2]
                        symmetry_reduction = false
                       ) """

    positions = [[0.000000, 0.000000, 0.000000],
                 [0.500000, 0.500000, 0.500000],
                 [0.695526, 0.695526, 0.000000],
                 [0.304474, 0.304474, 0.000000],
                 [0.195526, 0.804474, 0.500000],
                 [0.804474, 0.195526, 0.500000]]

    fractional_strings = {}
    fractional_shift_strings = {}
    for i, position in enumerate(positions):
        fractional_strings['pos' + str(i)] = list_to_string(position, ',')
        shifted_position = [xyz + arb_shift for xyz in position]
        fractional_shift_strings['pos' + str(i)] = list_to_string(shifted_position, ',')

    # ** unpacks the dictionary
    input_noshift = input_string.format(**fractional_strings)
    input_shift = input_string.format(**fractional_shift_strings)
    return (input_noshift, input_shift)


def boron_nitride_hex_translational_invariance_strings(named_result, arb_shift):
    # Reference: https://materialsproject.org/materials/mp-984/
    input_string = named_result + """ := xtb(
         structure( fractional=[[B, {pos0:s}],
                                [B, {pos1:s}],
                                [N, {pos2:s}],
                                [N, {pos3:s}]]
                lattice( a = 2.51242804 angstrom
                         c = 7.70726501 angstrom
                         bravais = hexagonal
                        )
                  )
      repulsive_cutoff = 40.0 bohr
      overlap_cutoff = 40.0 bohr
      h0_cutoff = 40.0 bohr
      temperature = 0 kelvin
      ewald_real_cutoff = 10 bohr
      ewald_reciprocal_cutoff = 2
      ewald_alpha = 0.5
      monkhorst_pack = [2, 2, 2]
      symmetry_reduction = false
     ) """

    positions = [[0.333333, 0.666667, 0.250000],
                 [0.666667, 0.333333, 0.750000],
                 [0.333333, 0.666667, 0.750000],
                 [0.666667, 0.333333, 0.250000]]

    fractional_strings = {}
    fractional_shift_strings = {}
    for i, position in enumerate(positions):
        fractional_strings['pos'+str(i)] = list_to_string(position, ',')
        shifted_position = [xyz + arb_shift for xyz in position]
        fractional_shift_strings['pos'+str(i)] = list_to_string(shifted_position, ',')

    input_noshift = input_string.format(**fractional_strings)
    input_shift = input_string.format(**fractional_shift_strings)

    return (input_noshift, input_shift)


# Fill in skeleton routines
#
# def diamond_translational_invariance():
#     return abs(energy[1] - energy[0])
#
# def sodium_chloride():
#     return abs(energy[1] - energy[0])

# def anatase():
#     return abs(energy[1] - energy[0])
#
# def calcium_carbonate():
#     return abs(energy[1] - energy[0])
#
# def lead_sulphide():
#     return abs(energy[1] - energy[0])
#
# graphite:
# Looks like CIF contains asymmetric cell and vesta converts to conventional
# Need primitive cell


# -----------------------------------------------
# Generate entos input files
# -----------------------------------------------

def boron_nitride_input_str():
    named_result = 'bn_hex'
    (input1, input2) = boron_nitride_hex_translational_invariance_strings(named_result, arb_shift=0.13)
    assert_string = " assert(load = " + named_result + " variable = n_iterations value = 6)\n"
    assert_string += " assert(load = " + named_result + " variable = energy value =  -9.431841) \n\n"
    print(input1 + '\n' + assert_string + input2 + assert_string)

if generate_app_test_inputs:
    boron_nitride_input_str()


# -----------------------------------------------
# Compute energy differences
# -----------------------------------------------

if run_energy_differences:
    energy_difference = collections.OrderedDict()

    energy_difference['silicon'] = get_energy_differences("silicon",
                                                          silicon_translational_invariance_strings,
                                                          arbitrary_shift=0.43)[0]

    energy_difference['rutile'] = get_energy_differences("rutile",
                                                         rutile_translational_invariance_strings,
                                                         arbitrary_shift=0.134)[0]

    energy_difference['bn_hex'] = get_energy_differences("bn_hex",
                                                         boron_nitride_hex_translational_invariance_strings,
                                                         arbitrary_shift=0.13)[0]

    all_pass = evaluate_energy_differences(energy_difference)
    assert all_pass








