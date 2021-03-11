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


def o_basis(energy_cutoff: float, max_matching_order: int):

    prefix = """
    <?xml version="1.0" encoding="UTF-8"?>
    <spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
      <sp chemicalSymbol="O" name="oxygen" z="-8.00000" mass="29165.12203">
        <muffinTin rmin="0.100000E-05" radius="1.4500" rinf="17.0873" radialmeshPoints="400"/>
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

    o_l0_trial_energies = np.array([
        5.80927914521301,
        16.7060568693384,
        31.9254693308077,
        51.2471465637003,
        74.5830334054306,
        101.882872176025
    ])

    o_l0_trial_energies = o_l0_trial_energies[o_l0_trial_energies < energy_cutoff]

    o_l1_trial_energies = np.array([
        5.29836476840879,
        14.5810416642800,
        27.9915249042968,
        45.4457829750705,
        66.8835003001113,
        92.2609075875359,
        121.558822945226
    ])

    o_l1_trial_energies = o_l1_trial_energies[o_l1_trial_energies < energy_cutoff]

    l0_extras = generate_lo_l_string(0, o_l0_trial_energies, 1)
    l1_extras = generate_lo_l_string(1, o_l0_trial_energies, 1)

    return prefix + lo_0 + l0_extras + lo_1 + l1_extras + end

