# Translate the atoms in the primitive cell
# by a fixed amount and confirm that the total energy is unchanged

import json
import subprocess
import collections

def list_to_string(mylist: list, joiner: str) -> str:
    return joiner.join([str(i) for i in mylist])

entos_exe = '/Users/alexanderbuccheri/Codes/entos/cmake-build-debug/entos'

def silicon_translational_invariance():
    input_string = """bulk_si := xtb(
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

    arbitrary_shift = 0.43
    fractional_0_shifted = [round(i + arbitrary_shift, 3) for i in fractional_0]
    fractional_1_shifted = [round(i + arbitrary_shift, 3) for i in fractional_1]

    fractional_0_str = list_to_string(fractional_0, ',')
    fractional_1_str = list_to_string(fractional_1, ',')
    input_noshift = input_string.format(pos0=fractional_0_str, pos1=fractional_1_str).replace('\n', ' ')

    fractional_0_shifted_str = list_to_string(fractional_0_shifted, ',')
    fractional_1_shifted_str = list_to_string(fractional_1_shifted, ',')
    input_shift = input_string.format(pos0=fractional_0_shifted_str, pos1=fractional_1_shifted_str).replace('\n', ' ')

    energy = []
    for input in [input_noshift, input_shift]:
        entos_command = [entos_exe, '--format', 'json', '-s', input]
        entos_output = subprocess.check_output(entos_command)
        result = json.loads(entos_output)
        energy.append(result['bulk_si']['energy'])

    return abs(energy[1] - energy[0])


def rutile_translational_invariance():
    # Reference: https: // materialsproject.org / materials / mp - 2657 /
    input_string = """bulk_rutile := xtb( 
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

    arb_shift = 0.134

    fractional_strings = {}
    fractional_shift_strings = {}
    for i, position in enumerate(positions):
        fractional_strings['pos'+str(i)] = list_to_string(position, ',')
        shifted_position = [xyz + arb_shift for xyz in position]
        fractional_shift_strings['pos'+str(i)] = list_to_string(shifted_position, ',')

    # ** unpacks the dictionary
    input_noshift = input_string.format(**fractional_strings).replace('\n', ' ')
    input_shift = input_string.format(**fractional_shift_strings).replace('\n', ' ')

    energy = []
    for input in [input_noshift, input_shift]:
        entos_command = [entos_exe, '--format', 'json', '-s', input]
        entos_output = subprocess.check_output(entos_command)
        result = json.loads(entos_output)
        energy.append(result['bulk_rutile']['energy'])
        print(energy)

    return abs(energy[1] - energy[0])


#TODO(Alex) Fill in skeleton routines
# def diamond_translational_invariance():
#     return abs(energy[1] - energy[0])
#
# def boron_nitride_hex():
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


# Compute energy differences
energy_difference = collections.OrderedDict()

energy_difference['silicon'] = silicon_translational_invariance()

energy_difference['rutile'] = rutile_translational_invariance()

# Evaluate energy differences
for key, value in energy_difference:
    if energy_difference > 1.e-6:
        print(key + ' primitive unit cell')
        print("Energy should be invariant w.r.t. translation of basis positions,"
              "however we find |delta E| (Ha) = ", value)




