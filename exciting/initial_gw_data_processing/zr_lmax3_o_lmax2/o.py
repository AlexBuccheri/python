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
      <custom l="2" type="lapw" trialEnergy="0.5" searchE="true"/>
      <lo l="2">
        <wf matchingOrder="0" trialEnergy="0.5" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.5" searchE="true"/>
      </lo>

      <lo l="2">
        <wf matchingOrder="1" trialEnergy="0.5" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.5" searchE="true"/>
      </lo>
"""

o_l0_trial_energies = np.array([
5.80933432136286,
16.7061133526094,
31.9255256532359,
51.2472022710829,
74.5830889201517,
101.882927329434,
133.111994727049,
168.251921800440,
207.293972838563,
250.229957037615,
297.051825519260,
347.754361276701,
402.334374197427,
460.788883319951,
523.115086098318,
589.310835158727,
659.374457069366,
733.304461959578,
811.099536279284
])

o_l1_trial_energies = np.array([
5.29842200275321,
14.5810985652875,
27.9915801318094,
45.4458377595188,
66.8835554524094,
92.2609632916323,
121.558878410012,
154.769211953336,
191.880976504408,
232.883458378917,
277.771146978118,
326.540812714950,
379.188563170813,
435.710774565492,
496.104999662198,
560.369343348915,
628.501988265783,
700.501345254497,
776.366124504543,
856.095212327207
])

o_l2_trial_energies = np.array([
10.8316055215777,
22.2926163420635,
37.7901440867745,
57.2635002866576,
80.6928321035865,
108.058179381480,
139.335827097075,
174.512012732428,
213.582148658821,
256.541208582489,
303.382581540232,
354.101780800580,
408.696381894463,
467.164047229256,
529.502352105417,
595.709396621513,
665.783726742382,
739.724028987369,
817.529116132766,
899.197977366254
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
                                                   o_l2_trial_energies],  energy_cutoff)

    lo_strings = []
    for l_value, trial_energyies_array in enumerate(filtered_trial_energies):
        lo_strings.append(generate_lo_l_string(l_value, trial_energyies_array, max_matching_order))

    return prefix + lo_0 + lo_strings[0] \
                 + lo_1 + lo_strings[1] \
                 + lo_2 + lo_strings[2] \
                 + end

