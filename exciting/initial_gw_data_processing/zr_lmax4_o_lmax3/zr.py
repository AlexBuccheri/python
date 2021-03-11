import numpy as np

prefix = """
<?xml version="1.0" encoding="UTF-8"?>
<spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <sp chemicalSymbol="Zr" name="zirconium" z="-40.0000" mass="166291.1791">
    <muffinTin rmin="0.100000E-05" radius="2.0000" rinf="26.3269" radialmeshPoints="600"/>
    <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="3" l="1" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="3" l="1" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="2" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="2" kappa="3" occ="6.00000" core="true"/>
    <atomicState n="4" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="4" l="1" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="4" l="1" kappa="2" occ="4.00000" core="false"/>
    <atomicState n="4" l="2" kappa="2" occ="2.00000" core="false"/>
    <atomicState n="5" l="0" kappa="1" occ="2.00000" core="false"/>
    <basis>

      <default type="lapw" trialEnergy="0.1500" searchE="false"/>
"""

end = """
    </basis>
  </sp>
</spdb>

"""

lo_0 = """
   <custom l="0" type="lapw" trialEnergy="-2.0" searchE="true"/>  
      <lo l="0"> 
	<wf matchingOrder="0" trialEnergy="-2.0" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="-2.0" searchE="true"/> 
      </lo> 
      <lo l="0"> 
	<wf matchingOrder="1" trialEnergy="-2.0" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="-2.0" searchE="true"/> 
      </lo> 

"""

lo_1 = """
    <custom l="1" type="lapw" trialEnergy="-0.5" searchE="true"/>  
      <lo l="1"> 
	<wf matchingOrder="0" trialEnergy="-0.5" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/> 
      </lo> 
      <lo l="1"> 
	<wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="-0.5" searchE="true"/> 
      </lo> 
"""

lo_2 = """
    <custom l="2" type="lapw" trialEnergy="0.33" searchE="true"/>  
      <lo l="2"> 
	<wf matchingOrder="0" trialEnergy="0.33" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> 
      </lo>
      <lo l="2"> 
	<wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="0.33" searchE="true"/> 
      </lo> 

"""

lo_3 = """
    <custom l="3" type="lapw" trialEnergy="1.00" searchE="true"/>  
      <lo l="3"> 
	<wf matchingOrder="0" trialEnergy="1.00" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="1.00" searchE="true"/> 
      </lo>
      <lo l="3"> 
	<wf matchingOrder="1" trialEnergy="1.00" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="1.00" searchE="true"/> 
      </lo> 
"""

lo_4 = """
      <custom l="4" type="lapw" trialEnergy="1.00" searchE="true"/>
      <lo l="4">
        <wf matchingOrder="0" trialEnergy="1.00" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="1.00" searchE="true"/>
      </lo>
      <lo l="4">
        <wf matchingOrder="1" trialEnergy="1.00" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="1.00" searchE="true"/>
      </lo>
"""


zr_l0_trial_energies = np.array([
3.41674973498303,
11.8576335922056,
23.5553488286426,
38.2225139492090,
55.7183934928437,
75.9631434089484,
98.9062300331434,
124.506902544624
])

zr_l1_trial_energies = np.array([
4.25309289713687,
12.6683041144246,
24.1783693122255,
38.5543138792516,
55.6730176346131,
75.4659088967451,
97.8905521589860,
122.912605026687
])

zr_l2_trial_energies = np.array([
5.63201190641682,
13.6742717606980,
24.5438162077787,
38.1197697668750,
54.3281060520844,
73.1368245214146,
94.5250199257823,
118.473406611847
])

zr_l3_trial_energies = np.array([
6.64673661942383,
13.6764927798014,
23.4106472587631,
35.8037851146779,
50.8128677385146,
68.4234854530309,
88.6093920574303,
111.352631210723
])

zr_l4_trial_energies = np.array([
11.1770562687516,
20.0456009666540,
31.4363977820287,
45.3837373317974,
61.8787737989140,
80.9149255382461,
102.493764412311,
126.606730752395
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


def zr_basis(energy_cutoff: float, max_matching_order: int):

    filtered_trial_energies = filter_trial_energies([zr_l0_trial_energies,
                                                   zr_l1_trial_energies,
                                                   zr_l2_trial_energies,
                                                   zr_l3_trial_energies,
                                                   zr_l4_trial_energies],  energy_cutoff)

    lo_strings = []
    for l_value, trial_energyies_array in enumerate(filtered_trial_energies):
        lo_strings.append(generate_lo_l_string(l_value, trial_energyies_array, max_matching_order))

    return prefix + lo_0 + lo_strings[0] \
                 + lo_1 + lo_strings[1] \
                 + lo_2 + lo_strings[2] \
                 + lo_3 + lo_strings[3] \
                 + lo_4 + lo_strings[4] \
                 + end
