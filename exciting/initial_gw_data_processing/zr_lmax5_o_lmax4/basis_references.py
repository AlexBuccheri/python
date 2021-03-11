import numpy as np

# Copy-paste the ground state basis here, then add tags and replace with
# o_basis.format(** {custom_lo0: lo0_string, custom_lo1: lo1_string})

oxygen_groundstate_basis = """
<spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <sp chemicalSymbol="O" name="oxygen" z="-8.00000" mass="29165.12203">
    <muffinTin rmin="0.100000E-05" radius="1.4500" rinf="17.0873" radialmeshPoints="600"/>
    <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="2" l="1" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="2" l="1" kappa="2" occ="2.00000" core="false"/>
    <basis>
      <default type="lapw" trialEnergy="0.1500" searchE="false"/>

      <custom l="0" type="lapw" trialEnergy="0.1500" searchE="true"/>
      <lo l="0">
        <wf matchingOrder="0" trialEnergy="-0.5" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/>
      </lo>
      <lo l="0">
        <wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="-0.5" searchE="true"/>
      </lo>

      {custom_l0}

      <custom l="1" type="lapw" trialEnergy="0.1" searchE="true"/>
      <lo l="1">
        <wf matchingOrder="0" trialEnergy="0.1" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.1" searchE="true"/>
      </lo>
      <lo l="1">
        <wf matchingOrder="1" trialEnergy="0.1" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.1" searchE="true"/>
      </lo>

      {custom_l1}

      <custom l="2" type="lapw" trialEnergy="0.15" searchE="true"/>
      <lo l="2">
        <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
      </lo>

      <lo l="2">
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/>
      </lo>

      {custom_l2}

      <custom l="3" type="lapw" trialEnergy="0.15" searchE="true"/>
      <lo l="3">
        <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
      </lo>

      <lo l="3">
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/>
      </lo>

      {custom_l3}

      <custom l="4" type="lapw" trialEnergy="0.15" searchE="true"/>
      <lo l="4">
        <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
      </lo>

      <lo l="4">
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/>
      </lo>

      {custom_l4}

    </basis>
  </sp>
</spdb>

"""

zr_groundstate_basis = """
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
      <custom l="0" type="lapw" trialEnergy="-2.0" searchE="true"/>  
      <lo l="0"> 
	<wf matchingOrder="0" trialEnergy="-2.0" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="-2.0" searchE="true"/> 
      </lo> 
      <lo l="0"> 
	<wf matchingOrder="1" trialEnergy="-2.0" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="-2.0" searchE="true"/> 
      </lo> 

     {custom_l0} 

      <custom l="1" type="lapw" trialEnergy="-0.5" searchE="true"/>  
      <lo l="1"> 
	<wf matchingOrder="0" trialEnergy="-0.5" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/> 
      </lo> 
      <lo l="1"> 
	<wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="-0.5" searchE="true"/> 
      </lo> 

    {custom_l1} 

      <custom l="2" type="lapw" trialEnergy="0.33" searchE="true"/>  
      <lo l="2"> 
	<wf matchingOrder="0" trialEnergy="0.33" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> 
      </lo>
      <lo l="2"> 
	<wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="0.33" searchE="true"/> 
      </lo> 

    {custom_l2}

      <custom l="3" type="lapw" trialEnergy="0.15" searchE="true"/>  
      <lo l="3"> 
	<wf matchingOrder="0" trialEnergy="0.15" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="0.15" searchE="true"/> 
      </lo>
      <lo l="3"> 
	<wf matchingOrder="1" trialEnergy="0.15" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="0.15" searchE="true"/> 
      </lo> 

    {custom_l3}

      <custom l="4" type="lapw" trialEnergy="0.15" searchE="true"/>
      <lo l="4">
        <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
      </lo>
      <lo l="4">
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/>
      </lo>

    {custom_l4}

      <custom l="5" type="lapw" trialEnergy="0.15" searchE="true"/>
      <lo l="5">
        <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
      </lo>
      <lo l="5">
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/>
      </lo>

    {custom_l5}

    </basis>
  </sp>
</spdb>
"""

# LINENGY.OUT. A reference to one of 4 functions per l-value
# Trial energy the same in all cases
#
# Zr
#  local-orbital functions :
#   l.o. =  1, l =  0, order =  1 :   -1.390000000
#   l.o. =  3, l =  1, order =  1 :  -0.5100000000
#   l.o. =  5, l =  2, order =  1 :   0.3300000000
#   l.o. =  7, l =  3, order =  1 :    1.000000000
#   l.o. =  9, l =  4, order =  1 :    1.000000000
#   l.o. = 11, l =  5, order =  1 :    1.000000000

# Use all the lorecommendations for both species, for any trial energy that exceeds the ones in the groundstate basis
# and go up to something reasonable like 200 Ha (not testing above that)

zr_l0_param_energies = np.array([
 3.41680898437054,
 11.8576938160972,
 23.5554092179115,
 38.2225743653328,
 55.7184541414946,
 75.9632042132661,
 98.9062907682391,
 124.506963631815,
 152.736382049958,
 183.572330514498,
 216.996258309127
])

zr_l1_param_energies = np.array([
4.25315291890153,
12.6683646478874,
24.1784297417876,
38.5543743382504,
55.6730783216737,
75.4659697341828,
97.8906129248615,
122.912665751060,
150.509452639976,
180.663672344259,
213.361354872195
])

zr_l2_param_energies = np.array([
5.63207202422344,
13.6743324027838,
24.5438767537143,
38.1198303318921,
54.3281665696292,
73.1368854195481,
94.5250807411363,
118.473467377649,
144.970417790885,
174.004563436752,
205.566772757493,
239.649724026780
])

zr_l3_param_energies = np.array([
2.35656385982863,
6.64679733266455,
13.6765531715818,
23.4107076815681,
35.8038458268313,
50.8129281422818,
68.4235464148143,
88.6094529201154,
111.352692361235,
136.640729756977,
164.462154984826,
194.809958925107,
227.677675107155
])

zr_l4_param_energies = np.array([
4.77455121020908,
11.1771167783142,
20.0456615354422,
31.4364582614499,
45.3837978604668,
61.8788345591218,
80.9149864427151,
102.493825240835,
126.606791532102,
153.247775037703,
182.409893765899
])

zr_l5_param_energies = np.array([
7.18997802942021,
15.3532461395206,
25.8441635227260,
38.7938609014216,
54.2379030848716,
72.1950076024806,
92.6583889975120,
115.630977693237,
141.114673179722,
169.107320134930,
199.607716989205,
232.612589814082
])

zr_trial_energies = [zr_l0_param_energies,
                     zr_l1_param_energies,
                     zr_l2_param_energies,
                     zr_l3_param_energies,
                     zr_l4_param_energies,
                     zr_l5_param_energies]

# LINENGY.OUT. A reference to one of 4 functions per l-value
# Trial energy the same in all functions with a given l-value
#
# O
#  local-orbital functions :
#   l.o. =  1, l =  0, order =  1 :  -0.4125000000E-01
#   l.o. =  3, l =  1, order =  1 :   0.1000000000
#   l.o. =  5, l =  2, order =  1 :    1.000000000
#   l.o. =  7, l =  3, order =  1 :    1.000000000
#   l.o. =  9, l =  4, order =  1 :    1.000000000

o_l0_param_energies = np.array([
5.80945052840072,
16.7062300478316,
31.9256424986285,
51.2473193415797,
74.5832061082049,
101.883044651259,
133.112112213231,
168.252039431831,
207.294090299430
])

o_l1_param_energies = np.array([
0.533864174665581,
5.29853864332659,
14.5812155082034,
27.9916971277413,
45.4459545917901,
66.8836726867726,
92.2610806819811,
121.558995913797,
154.769329550561,
191.881094229837,
232.883576256490
])

o_l2_param_energies = np.array([
3.53299155030141,
10.8317220427379,
22.2927328327085,
37.7902610059448,
57.2636174281591,
80.6929493429828,
108.058296728187,
139.335944599908,
174.512130382165,
213.582266432544
])

o_l3_param_energies = np.array([
6.30671604138950,
16.0309187705152,
29.6509677063678,
47.2373408134373,
68.7553006190284,
94.1889706920050,
123.535948355816,
156.789575073662,
193.937150415444,
234.969897767763
])

o_l4_param_energies = np.array([
9.40057261707854,
21.4575643161227,
37.1973890071254,
56.8411287418714,
80.3902745607293,
107.832442393812,
139.164029808692,
174.387320333020,
213.500065365429,
256.495012493561
])

o_trial_energies = [o_l0_param_energies,
                    o_l1_param_energies,
                    o_l2_param_energies,
                    o_l3_param_energies,
                    o_l4_param_energies]