"""
Inputs for GW Scaling calculations

q = [2,2,2] and img frequency = 45 such that n_q * n_f = 360 (cores over 10 nodes).

nempty=100, such that the calculations are fast
60 Ha cut-off in LOs

Ref directory for basis functions:
/users/sol/abuccheri/gw_benchmarks/A1/zr_lmax3_o_lmax2_rgkmax7/gw_q222_omeg32_nempty800/max_energy_60
"""

input_xml ="""
<input>

   <title>A1-ZrO2-Cubic-Primitive-PBEsol</title>

   <structure speciespath=".">
      <crystal  scale="1.8897259886">
         <basevect>0.00000000   2.53574055   2.53574055</basevect>
         <basevect>2.53574055   0.00000000   2.53574055</basevect>
         <basevect>2.53574055   2.53574055   0.00000000</basevect>
      </crystal>

      <species speciesfile="Zr.xml" rmt="2.0">
         <atom coord="0.00 0.00 0.00"></atom>
      </species>

      <species speciesfile="O.xml" rmt="1.6">
         <atom coord="0.25 0.25 0.25"></atom>
         <atom coord="0.75 0.75 0.75"></atom>
      </species>
   </structure>

   <groundstate
      do="fromscratch"
      rgkmax="7.0"
      ngridk="6 6 6"
      xctype="GGA_PBE_SOL"
      epsengy="1.e-8"
      gmaxvr="20.0"
      >
   </groundstate>

   <gw
    taskname="g0w0"
    nempty="100"
    ngridq="2 2 2"
    skipgnd="false"
    >
  
    <mixbasis
      lmaxmb="4"
      epsmb="1.d-3"
      gmb="1.d0"
    ></mixbasis>
  
    <freqgrid
      nomeg="45"
      freqmax="1.0"
    ></freqgrid>
  
    <barecoul
      pwm="2.0"
      stctol="1.d-16"
      barcevtol="0.1"
    ></barecoul>
  
    <selfenergy
      actype="pade"
      singularity="mpb"
    ></selfenergy>
  
   </gw>

</input>
"""


# 60 Ha cut-off in LOs
zr_basis_xml = """
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

      
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="3.42" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="3.42" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="11.86" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="11.86" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="23.56" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="23.56" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="38.22" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="38.22" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="55.72" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="55.72" searchE="false"/>
       </lo>
         

      <custom l="1" type="lapw" trialEnergy="-0.5" searchE="true"/>  
      <lo l="1"> 
	<wf matchingOrder="0" trialEnergy="-0.5" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/> 
      </lo> 
      <lo l="1"> 
	<wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="-0.5" searchE="true"/> 
      </lo> 

      
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="4.25" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="4.25" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="12.67" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="12.67" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="24.18" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="24.18" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="38.55" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="38.55" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="55.67" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="55.67" searchE="false"/>
       </lo>
         

      <custom l="2" type="lapw" trialEnergy="0.33" searchE="true"/>  
      <lo l="2"> 
	<wf matchingOrder="0" trialEnergy="0.33" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> 
      </lo>
      <lo l="2"> 
	<wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="0.33" searchE="true"/> 
      </lo> 

      
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="5.63" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="5.63" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="13.67" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="13.67" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="24.54" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="24.54" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="38.12" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="38.12" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="54.33" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="54.33" searchE="false"/>
       </lo>
         

      <custom l="3" type="lapw" trialEnergy="1.00" searchE="true"/>  
      <lo l="3"> 
	<wf matchingOrder="0" trialEnergy="1.00" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="1.00" searchE="true"/> 
      </lo>
      <lo l="3"> 
	<wf matchingOrder="1" trialEnergy="1.00" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="1.00" searchE="true"/> 
      </lo> 

      
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="6.65" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="6.65" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="13.68" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="13.68" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="23.41" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="23.41" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="35.8" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="35.8" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="50.81" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="50.81" searchE="false"/>
       </lo>
         

    </basis>
  </sp>
</spdb>
"""

# 60 Ha cut-off in LOs
o_basis_xml = """
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

      
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="5.81" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="5.81" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="16.71" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="16.71" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="31.93" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="31.93" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="51.25" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="51.25" searchE="false"/>
       </lo>
         

      <custom l="1" type="lapw" trialEnergy="0.1" searchE="true"/>
      <lo l="1">
        <wf matchingOrder="0" trialEnergy="0.1" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.1" searchE="true"/>
      </lo>
      <lo l="1">
        <wf matchingOrder="1" trialEnergy="0.1" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.1" searchE="true"/>
      </lo>

      
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="5.3" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="5.3" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="14.58" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="14.58" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="27.99" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="27.99" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="45.45" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="45.45" searchE="false"/>
       </lo>
         

      <custom l="2" type="lapw" trialEnergy="0.10" searchE="true"/>
      <lo l="2">
        <wf matchingOrder="0" trialEnergy="0.10" searchE="true"/>
        <wf matchingOrder="1" trialEnergy="0.10" searchE="true"/>
      </lo>

      <lo l="2">
        <wf matchingOrder="1" trialEnergy="0.10" searchE="true"/>
        <wf matchingOrder="2" trialEnergy="0.10" searchE="true"/>
      </lo>

      
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="10.83" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="10.83" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="22.29" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="22.29" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="37.79" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="37.79" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="57.26" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="57.26" searchE="false"/>
       </lo>
         

    </basis>
  </sp>
</spdb>
"""