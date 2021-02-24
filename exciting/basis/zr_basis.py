import numpy as np


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


def zr_basis(energy_cutoff: float, max_matching_order: int):

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
      <custom l="1" type="lapw" trialEnergy="-1.2" searchE="true"/>
      
      <lo l="1">
	     <wf matchingOrder="0" trialEnergy="-1.2" searchE="true"/>
	     <wf matchingOrder="1" trialEnergy="-1.2" searchE="true"/>
      </lo>
      <lo l="1">
	     <wf matchingOrder="1" trialEnergy="-1.2" searchE="true"/>
	     <wf matchingOrder="2" trialEnergy="-1.2" searchE="true"/>
      </lo>

    """

    lo_2 = """
      <custom l="2" type="lapw" trialEnergy="0.15" searchE="true"/>
      
      <lo l="2">
	     <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/>
	     <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
      </lo>
      <lo l="2">
	     <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
	     <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/>
      </lo>

    """

    zr_l0_trial_energies = np.array([
        3.41664190279722,
        11.8575202888064,
        23.5552337618056,
        38.2223979680527,
        55.7182768821595,
        75.9630266420627,
        98.9061129488830
    ])

    zr_l0_trial_energies = zr_l0_trial_energies[zr_l0_trial_energies< energy_cutoff]

    zr_l1_trial_energies = np.array([
        4.25298375036404,
        12.6681901859002,
        24.1782538890203,
        38.5541975851575,
        55.6729010214330,
        75.4657918893579,
        97.8904348668412
    ])

    zr_l1_trial_energies = zr_l1_trial_energies[zr_l1_trial_energies< energy_cutoff]

    zr_l2_trial_energies = np.array([
        5.63189923196280,
        13.6741564955162,
        24.5436998681673,
        38.1196527377421,
        54.3279885923796,
        73.1367070205889,
        94.5249022361277
    ])

    zr_l2_trial_energies = zr_l2_trial_energies[zr_l2_trial_energies< energy_cutoff]

    l0_extras = generate_lo_l_string(0, zr_l0_trial_energies, max_matching_order)
    l1_extras = generate_lo_l_string(1, zr_l1_trial_energies, max_matching_order)
    l2_extras = generate_lo_l_string(2, zr_l2_trial_energies, max_matching_order)


    return prefix + lo_0 + l0_extras + lo_1 + l1_extras + lo_2 + l2_extras + end
