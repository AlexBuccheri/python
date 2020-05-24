# Translate the atoms in the primitive cell
# by a fixed amount and confirm that the total energy is unchanged
# Note, arbitrary shift should be small enough that all atoms stay within the unit cell

import numpy as np
import json
import subprocess
import collections
import typing


entos_root = '/Users/alexanderbuccheri/Codes/entos/'
entos_exe = entos_root + 'cmake-build-debug/entos'

# Run entos energy differences from python
run_energy_differences = False

# Generate equivalent app test inputs from python
generate_app_test_inputs = False


# -----------------------------------
# Functions
# -----------------------------------

def list_to_string(mylist: list, joiner: str) -> str:
    return joiner.join([str(i) for i in mylist])


def get_energies(named_result: str, inputs: typing.List[str]) -> typing.List[float]:
    energy = []
    for input in inputs:
        entos_command = [entos_exe, '--format', 'json', '-s', input.replace('\n', ' ')]
        entos_output = subprocess.check_output(entos_command)
        result = json.loads(entos_output)
        energy.append(result[named_result]['energy'])
    return energy


def get_energy_differences(named_result: str,
                           strings_function: typing.Callable,
                           arbitrary_shift: float) -> typing.List[float]:

    input_strings = strings_function(named_result, arbitrary_shift)
    energies = get_energies(named_result, input_strings)

    energy_diffs = []
    for i in range(1, len(energies)):
        energy_diffs.append(abs(energies[i] - energies[i-1]))

    return energy_diffs


def evaluate_energy_differences(energy_difference: typing.Dict, tol=1.e-6) -> bool:
    all_pass = True
    for key, value in energy_difference.items():
        if value > tol:
            print(key + ' primitive unit cell')
            print("Energy should be invariant w.r.t. translation of basis positions,"
                  "however we find |delta E| (Ha) = ", value)
            all_pass = False
    return all_pass


def check_shifted_positions_witin_cell(named_result, positions: typing.List, arb_shift:float):

    all_positions_in_cell = True

    for i,position in enumerate(positions):
        xyz = np.asarray(position) + arb_shift
        below_zero = np.where((xyz < 0))[0]
        exceed_one = np.where((xyz > 1))[0]
        if len(below_zero) > 0:
            print("Component/s of atom position " + str(i) + " of " + named_result + " below 0: " + str(xyz))
            all_positions_in_cell = False
        elif len(exceed_one) > 0:
            print("Component/s of atom position " + str(i) + " of " + named_result + " exceed 1: " + str(xyz))
            all_positions_in_cell = False

    return all_positions_in_cell


def substitute_positions_in_input(input_string: str, positions: typing.List, arb_shift: float) \
        -> typing.Tuple[str]:
    """ Fills in positions using string.format

        Inputs:
          positions. List of positions. Must be in fractional coordinates

        Returns a tuple containing string with a) original positions subbed and b) shifted positions subbed
    """

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


# ----------------------------------------
# entos string functions for translations
# ----------------------------------------

def silicon_translational_invariance_strings(named_result: str, arb_shift: float) -> typing.Tuple[str]:
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

    positions = [[0.00, 0.00, 0.00],[0.25, 0.25, 0.25]]

    assert check_shifted_positions_witin_cell(named_result, positions, arb_shift)
    return substitute_positions_in_input(input_string, positions, arb_shift)


def rutile_translational_invariance_strings(named_result: str, arb_shift: float) -> typing.Tuple[str]:
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

    assert check_shifted_positions_witin_cell(named_result, positions, arb_shift)
    return substitute_positions_in_input(input_string, positions, arb_shift)


def anatase_translational_invariance_strings(named_result: str, arb_shift: float) -> typing.Tuple[str]:
    """ TiO2 Anatase. Space group:  I4_1/amd (141)
        Indirect band gap: 2.062 eV
        https://materialsproject.org/materials/mp-390/#

        SPGLib confirms this is the primitive cell and is the space group 141
        However, it looks body-centred cubic rather than body-centred tetragonal to me
        from inspection of the lattice vectors and parameters   """

    input_string = named_result + """ := xtb(
                       structure( fractional=[[Ti, {pos0:s}],
                                              [Ti, {pos1:s}],
                                              [ O, {pos2:s}],
                                              [ O, {pos3:s}],
                                              [ O, {pos4:s}],
                                              [ O, {pos5:s}]]
                              lattice( a = 5.55734663 angstrom
                                       c = 5.55734663 angstrom
                                       bravais = 'body centred tetragonal'
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

    positions = [[0.500000, 0.500000, 0.000000],
                 [0.250000, 0.750000, 0.500000],
                 [0.456413, 0.956413, 0.500000],
                 [0.706413, 0.706413, 0.000000],
                 [0.043587, 0.543587, 0.500000],
                 [0.293587, 0.293587, 0.000000]]

    assert check_shifted_positions_witin_cell(named_result, positions, arb_shift)
    return substitute_positions_in_input(input_string, positions, arb_shift)



def boron_nitride_hex_translational_invariance_strings(named_result: str, arb_shift: float) -> typing.Tuple[str]:
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

    assert check_shifted_positions_witin_cell(named_result, positions, arb_shift)
    return substitute_positions_in_input(input_string, positions, arb_shift)


def calcium_titanate_translational_invariance_strings(named_result: str, arb_shift: float) -> typing.Tuple[str]:
    # Ref: https://materialsproject.org/materials/mp-4019/
    input_string = named_result + """ := xtb(
         structure( fractional=[[Ca, {pos0:s}],  
                                [Ca, {pos1:s}], 
                                [Ca, {pos2:s}], 
                                [Ca, {pos3:s}], 
                                [Ti, {pos4:s}], 
                                [Ti, {pos5:s}], 
                                [Ti, {pos6:s}], 
                                [Ti, {pos7:s}], 
                                [O , {pos8:s}],
                                [O , {pos9:s}],
                                [O , {pos10:s}],
                                [O , {pos11:s}],
                                [O , {pos12:s}],
                                [O , {pos13:s}],
                                [O , {pos14:s}],
                                [O , {pos15:s}],
                                [O , {pos16:s}],
                                [O , {pos17:s}],
                                [O , {pos18:s}],
                                [O , {pos19:s}]] 
                lattice(  a=5.40444906 angstrom
                          b=5.51303112 angstrom
                          c=7.69713264 angstrom
                          bravais = orthorhombic
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

    positions = [[0.991521, 0.044799, 0.750000],
                 [0.491521, 0.455201, 0.250000],
                 [0.508479, 0.544799, 0.750000],
                 [0.008479, 0.955201, 0.250000],
                 [0.500000, 0.000000, 0.500000],
                 [0.000000, 0.500000, 0.500000],
                 [0.000000, 0.500000, 0.000000],
                 [0.500000, 0.000000, 0.000000],
                 [0.921935, 0.520580, 0.250000],
                 [0.421935, 0.979420, 0.750000],
                 [0.578065, 0.020580, 0.250000],
                 [0.078065, 0.479420, 0.750000],
                 [0.707456, 0.291917, 0.959281],
                 [0.207456, 0.208083, 0.040719],
                 [0.792544, 0.791917, 0.540719],
                 [0.292544, 0.708083, 0.459281],
                 [0.707456, 0.291917, 0.540719],
                 [0.207456, 0.208083, 0.459281],
                 [0.292544, 0.708083, 0.040719],
                 [0.792544, 0.791917, 0.959281]]

    assert check_shifted_positions_witin_cell(named_result, positions, arb_shift)
    return substitute_positions_in_input(input_string, positions, arb_shift)



def aluminum_hexathiohypodiphosphate_trans_invariance_strings(named_result: str, arb_shift: float) \
        -> typing.Tuple[str]:
    """   AlPS4. Simple orthorhombic. Space group: P222 [16]
          Indirect band gap: 2.634 eV
          https://materialsproject.org/materials/mp-27462/"""

    input_string = named_result + """ := xtb(
           structure( fractional=[[Al, {pos0:s}],  
                                  [Al, {pos1:s}], 
                                  [P,  {pos2:s}], 
                                  [P,  {pos3:s}], 
                                  [S,  {pos4:s}], 
                                  [S,  {pos5:s}], 
                                  [S,  {pos6:s}], 
                                  [S,  {pos7:s}], 
                                  [S,  {pos8:s}],
                                  [S,  {pos9:s}],
                                  [S,  {pos10:s}],
                                  [S,  {pos11:s}]] 
                  lattice(  a=5.71230345 angstrom
                            b=5.71644625 angstrom
                            c=11.46678755 angstrom
                            bravais = orthorhombic
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
                 [0.500000, 0.000000, 0.500000],
                 [0.000000, 0.500000, 0.000000],
                 [0.000000, 0.000000, 0.500000],
                 [0.197847, 0.276435, 0.101916],
                 [0.197847, 0.723565, 0.898084],
                 [0.802153, 0.276435, 0.898084],
                 [0.802153, 0.723565, 0.101916],
                 [0.776404, 0.800507, 0.601208],
                 [0.776404, 0.199493, 0.398792],
                 [0.223596, 0.800507, 0.398792],
                 [0.223596, 0.199493, 0.601208]]

    assert check_shifted_positions_witin_cell(named_result, positions, arb_shift)

    return substitute_positions_in_input(input_string, positions, arb_shift)


def copper_trans_invariance_strings(named_result: str, arb_shift: float) \
        -> typing.Tuple[str]:
    """   Cu FCC. Space group: FM-3m [225]
	      https://materialsproject.org/materials/mp-30/
	"""

    input_string = named_result + """ := xtb(
           structure( fractional=[[Cu, {pos0:s}]] 
                  lattice(  a=2.56061937 angstrom
                            bravais = fcc
                          )
                  )
        repulsive_cutoff = 40.0 bohr
        overlap_cutoff = 40.0 bohr
        h0_cutoff = 40.0 bohr
        ! Will also work at 0 k
        temperature = 100 kelvin
        ewald_real_cutoff = 10 bohr
        ewald_reciprocal_cutoff = 2
        ewald_alpha = 0.5
        monkhorst_pack = [2, 2, 2]
        symmetry_reduction = true
       ) """

    positions = [[0.000000, 0.000000, 0.000000]]
    assert check_shifted_positions_witin_cell(named_result, positions, arb_shift)
    return substitute_positions_in_input(input_string, positions, arb_shift)




# -----------------------------------------------
# Generate entos input files
# -----------------------------------------------

def get_assert_string(named_result: str, variables: typing.Dict) -> str:
    assert_string = ""
    for variable_name, variable in variables.items():
        assert_string += " assert(load = " + named_result + \
        " variable = " + variable_name + " value = " + str(variable) + "\n"
    return assert_string + "\n"


def silicon_input_str() -> str:
    named_result = 'silicon'
    assert_variables = {"n_iter":8, "energy": -3.669686}
    (input1, input2) = silicon_translational_invariance_strings(named_result, arb_shift=0.43)
    assert_string = get_assert_string(named_result, assert_variables)
    return input1 + '\n' + assert_string + input2 + '\n' + assert_string


def rutile_input_str() -> str:
    named_result = 'rutile'
    assert_variables = {"n_iter":35, "energy": -150.135037}
    (input1, input2) = rutile_translational_invariance_strings(named_result, arb_shift=0.134)
    assert_string = get_assert_string(named_result, assert_variables)
    return input1 + '\n' + assert_string + input2 + '\n' + assert_string


def anatase_input_str() -> str:
    named_result = 'anatase'
    assert_variables = {"n_iter":11, "energy": -19.142514}
    (input1, input2) = anatase_translational_invariance_strings(named_result, arb_shift=0.04)
    assert_string = get_assert_string(named_result, assert_variables)
    return input1 + '\n' + assert_string + input2 + '\n' + assert_string


def boron_nitride_hex_input_str() -> str:
    named_result = 'bn_hex'
    comments = "! Periodic D3 should really be included for an accurate total energy \n"
    (input1, input2) = boron_nitride_hex_translational_invariance_strings(named_result, arb_shift=0.13)
    assert_variables = {"n_iter":6, "energy": -9.431841}
    assert_string = get_assert_string(named_result, assert_variables)
    return comments + input1 + '\n' + assert_string + input2 + '\n' + assert_string


def calcium_titanate_input_str() -> str:
    named_result = 'CaTiO3'
    (input1, input2) = calcium_titanate_translational_invariance_strings(named_result, arb_shift=0.04)
    assert_variables = {"n_iter":0, "energy":0}
    assert_string = get_assert_string(named_result, assert_variables)
    return input1 + '\n' + assert_string + input2 + '\n' + assert_string


def aluminum_hexathiohypodiphosphate_input_str() -> str:
    named_result = 'AlPS4'
    (input1, input2) = aluminum_hexathiohypodiphosphate_trans_invariance_strings(named_result, arb_shift=0.1)
    assert_variables = {"n_iter":0, "energy":0}
    assert_string = get_assert_string(named_result, assert_variables)
    return input1 + '\n' + assert_string + input2 + '\n' + assert_string


def copper_input_str() -> str:
    named_result = 'Copper'
    (input1, input2) = copper_trans_invariance_strings(named_result, arb_shift=0.53)
    assert_variables = {"n_iter":3, "energy":-3.827222, "entropy":-0.0001098}
    assert_string = get_assert_string(named_result, assert_variables)
    return input1 + '\n' + assert_string + input2 + '\n' + assert_string


if generate_app_test_inputs:
    string_functions = [silicon_input_str(),
                        rutile_input_str(),
                        boron_nitride_hex_input_str(),
                        copper_input_str(),
                        anatase_input_str()]
    input_string = "! Translation invariance of the Total Energy. \n" \
                   "! For each material, do two calculations and assert that the energies are the same\n\n" \
                   "! Autogenerated by translation_invariance.py. Material refs within"

    for material_input_str in string_functions:
        input_string += material_input_str + "\n\n"
    fid = open(entos_root + 'test/xtb_periodic_trans_inv.in', 'w')
    fid.write(input_string)
    fid.close()
    print("Exceptions:" )
    print("aluminum_hexathiohypodiphosphate gets stuck in a local minimum and does not converge")
    print("calcium_titanate matrix loses symmetry after ~ 7 iterations. Need to look into this")



# ----------------------------------------------------
# Compute energy differences via python calling entos
# ----------------------------------------------------

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








