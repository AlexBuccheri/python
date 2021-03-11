import numpy as np

prefix = """
<?xml version="1.0" encoding="UTF-8"?>
<spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <sp chemicalSymbol="O" name="oxygen" z="-8.00000" mass="29165.12203">
    <muffinTin rmin="0.100000E-05" radius="1.4500" rinf="17.0873" radialmeshPoints="600"/>
    <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="2" l="1" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="2" l="1" kappa="2" occ="2.00000" core="false"/>
    <basis>
      <default type="lapw" trialEnergy="0.1500" searchE="false"/>

"""

end = """
    </basis>
  </sp>
</spdb>

"""

lo_0 = """
      <custom l="0" type="lapw" trialEnergy="0.1500" searchE="true"/>
      <lo l="0">
        <wf matchingOrder="0" trialEnergy="-0.5" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/>
      </lo>
      <lo l="0">
        <wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="-0.5" searchE="true"/>
      </lo>
"""

lo_1 = """
      <custom l="1" type="lapw" trialEnergy="0.1" searchE="true"/>
      <lo l="1">
        <wf matchingOrder="0" trialEnergy="0.1" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.1" searchE="true"/>
      </lo>
      <lo l="1">
        <wf matchingOrder="1" trialEnergy="0.1" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.1" searchE="true"/>
      </lo>
"""

lo_2 = """
      <custom l="2" type="lapw" trialEnergy="1.00" searchE="true"/>
      <lo l="2">
        <wf matchingOrder="0" trialEnergy="1.00" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="1.00" searchE="true"/>
      </lo>

      <lo l="2">
        <wf matchingOrder="1" trialEnergy="1.00" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="1.00" searchE="true"/>
      </lo>
"""

lo_3 = """
      <custom l="3" type="lapw" trialEnergy="3.00" searchE="true"/>
      <lo l="3">
        <wf matchingOrder="0" trialEnergy="3.00" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="3.00" searchE="true"/>
      </lo>

      <lo l="3">
        <wf matchingOrder="1" trialEnergy="3.00" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="3.00" searchE="true"/>
      </lo>
"""
#Replace? the 1.00 trial energy with 6.30665393833073,  # This was the lowest energy => I wonder if the current lo=3 is wrong. Swith that out with this
# I'll just turn the guess up to 3

o_l0_trial_energies = np.array([
5.80938705814048,
16.7061666034956,
31.9255790520402,
51.2472557436106,
74.5831424365372,
101.882980904486,
133.112048374788
])

o_l1_trial_energies = np.array([
5.29847506754117,
14.5811519729642,
27.9916336056773,
45.4458909783032,
66.8836089941066,
92.2610168978132,
121.558932063140
])

o_l2_trial_energies = np.array([
10.8316585465262,
22.2926692634601,
37.7901975158305,
57.2635537887480,
80.6928856385164,
108.058232962616,
139.335880747533
])

o_l3_trial_energies = np.array([
16.0308558920391,
29.6509047060728,
47.2372775015344,
68.7552370920934,
94.1889070233316,
123.535884607105
])


def generate_lo_l_string(l, energies, max_matching_order):
    template = """
       <lo l="{l}">
        <wf matchingOrder="{mo1}" trialEnergy="{te}" searchE="false"/>
        <wf matchingOrder="{mo2}" trialEnergy="{te}" searchE="false"/>
       </lo>
        """
    string = ''
    matching_orders = [(i, i + 1) for i in range(0, max_matching_order)]
    for trial_energy in energies:
        for mo1, mo2 in matching_orders:
            string += template.format(l=l, te=round(trial_energy, 2), mo1=mo1, mo2=mo2)

    return string


def filter_trial_energies(trial_energies: list, energy_cutoff: float):
    """
    Expects a list of np arrays - each containing a set of trial energies
    w.r.t. a given l-value

    :param trial_energies:
    :return:
    """
    filtered_trial_energies = []
    for trial_energys_array in trial_energies:
        filtered_trial_energies.append(trial_energys_array[trial_energys_array < energy_cutoff])

    return filtered_trial_energies


def o_basis(energy_cutoff: float, max_matching_order: int):

    filtered_trial_energies = filter_trial_energies([o_l0_trial_energies,
                                                   o_l1_trial_energies,
                                                   o_l2_trial_energies,
                                                   o_l3_trial_energies],  energy_cutoff)

    lo_strings = []
    for l_value, trial_energyies_array in enumerate(filtered_trial_energies):
        lo_strings.append(generate_lo_l_string(l_value, trial_energyies_array, max_matching_order))

    return prefix + lo_0 + lo_strings[0] \
                 + lo_1 + lo_strings[1] \
                 + lo_2 + lo_strings[2] \
                 + lo_3 + lo_strings[3] \
                 + end
