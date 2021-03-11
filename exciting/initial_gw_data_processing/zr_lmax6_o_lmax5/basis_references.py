"""
Groundstate basis files for Zr l-max = 6 and O l-max = 5
with custom tags added
"""


oxygen_groundstate_basis = """
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

      <custom l="6" type="lapw" trialEnergy="0.15" searchE="true"/>
      <lo l="6">
        <wf matchingOrder="0" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
      </lo>
      <lo l="6">
        <wf matchingOrder="1" trialEnergy="0.15" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.15" searchE="true"/>
      </lo>
      
      {custom_l6}

    </basis>
  </sp>
</spdb>
"""
