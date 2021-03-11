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

lo_3 = """
     <custom l="3" type="lapw" trialEnergy="0.15" searchE="true"/>  
     <lo l="3"> 
	   <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/> 
	   <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/> 
     </lo>
     <lo l="3"> 
	   <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/> 
	   <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/> 
     </lo> 

"""

zr_l0_trial_energies = np.array([
3.41670315145522,
11.8575860474206,
23.5553007908067,
38.2224658045548,
55.7183451003712,
75.9630948999106,
98.9061811382636,
124.506853824274,
152.736272343468,
183.572220398408,
216.996148182596,
252.993582051900,
291.552301603424,
332.661903448781,
376.313508144067,
422.499286698384,
471.212374359221
])

zr_l1_trial_energies = np.array([
4.25304622502754,
12.6682565365447,
24.1783211999211,
38.5542651940606,
55.6729691523214,
75.4658603020841,
97.8905035375948,
122.912555855434,
150.509342851795,
180.663562164709,
213.361244689579,
248.592017973277,
286.347108457152,
326.619269248177,
369.402427158881,
414.691370914542,
462.481754830305,
512.769949416961
])

zr_l2_trial_energies = np.array([
5.63196429907165,
13.6742236645990,
24.5437678240075,
38.1197212744329,
54.3280573522054,
73.1367757409677,
94.5249711511503,
118.473357332827,
144.970307877388,
174.004453170129,
205.566662504965,
239.649613505682,
276.246713317352,
315.352778807108,
356.963486764696,
401.075309655956,
447.685328477261,
496.790960066153,
548.389835251048
])

zr_l3_trial_energies = np.array([
6.64668830805706,
13.6764440964337,
23.4105983964913,
35.8037366438940,
50.8128187612058,
68.4234367351880,
88.6093431245865,
111.352582177430,
136.640619737900,
164.462044651354,
194.809848623883,
227.677564312143,
263.059978533674,
300.952641608507,
341.351303087630,
384.252140702652,
429.651589870071,
477.546482836021,
527.934112857646,
580.812229375338
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
                                                   zr_l3_trial_energies],  energy_cutoff)

    lo_strings = []
    for l_value, trial_energyies_array in enumerate(filtered_trial_energies):
        lo_strings.append(generate_lo_l_string(l_value, trial_energyies_array, max_matching_order))

    return prefix + lo_0 + lo_strings[0] \
                 + lo_1 + lo_strings[1] \
                 + lo_2 + lo_strings[2] \
                 + lo_3 + lo_strings[3] \
                 + end

