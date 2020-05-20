# Translate the atoms in the primitive cell
# by a fixed amount and confirm that the total energy is unchanged

import json
import subprocess

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


    #TODO(Alex) Fill in skeleton routines
    def diamond_translational_invariance():
        return abs(energy[1] - energy[0])

    def boron_nitride_hex():
        return abs(energy[1] - energy[0])

    def sodium_chloride():
        return abs(energy[1] - energy[0])

    def anatase():
        return abs(energy[1] - energy[0])

    def calcium_carbonate():
        return abs(energy[1] - energy[0])

# Main Routine
energy_difference = silicon_translational_invariance()

if energy_difference > 1.e-6:
    print('Silicon primitive unit cell')
    print("Energy should be invariant w.r.t. translation of basis positions,"
          "however we find |delta E| (Ha) = ", energy_difference)

