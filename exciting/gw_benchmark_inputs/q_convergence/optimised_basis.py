"""
Optimised LO sets for ZrO2
"""

# 5s 12p 10d 8f   sourced from the respective channels in:
# 5s corresponds to i0    /users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8/s_channel/gw_q222_omeg32_nempty2000/max_energy_i0/Zr.xml
# 12p corresponds to i5   /users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8/p_channel/gw_q222_omeg32_nempty2000/max_energy_i5/Zr.xml
# 10d corresponds to i4   /users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8/d_channel/gw_q222_omeg32_nempty2000/max_energy_i4/Zr.xml
# 8f corresponds to i4    /users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8/f_channel/gw_q222_omeg32_nempty2000/max_energy_i4/Zr.xml
zr_basis = """<?xml version="1.0" encoding="UTF-8"?>
<spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <sp chemicalSymbol="Zr" name="zirconium" z="-40.0000" mass="166291.1791">
    <muffinTin rmin="0.100000E-05" radius="2.0000" rinf="26.3269" radialmeshPoints="600"/>
    <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="3" l="1" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="3" l="1" kappa="2" occ="4.00000" core="false"/>
    <atomicState n="3" l="2" kappa="2" occ="4.00000" core="false"/>
    <atomicState n="3" l="2" kappa="3" occ="6.00000" core="false"/>
    <atomicState n="4" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="4" l="1" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="4" l="1" kappa="2" occ="4.00000" core="false"/>
    <atomicState n="4" l="2" kappa="2" occ="2.00000" core="false"/>
    <atomicState n="5" l="0" kappa="1" occ="2.00000" core="false"/>
    <basis>
      <default type="lapw" trialEnergy="0.1500" searchE="false"/>

      <custom l="0" type="lapw" trialEnergy="-15.0" searchE="true"/>
      <lo l="0">
	    <wf matchingOrder="0" trialEnergy="-15.0" searchE="true"/>
	    <wf matchingOrder="1" trialEnergy="-15.0" searchE="true"/>
      </lo>
      <lo l="0">
	    <wf matchingOrder="1" trialEnergy="-15.0" searchE="true"/>
	    <wf matchingOrder="2" trialEnergy="-15.0" searchE="true"/>
      </lo>
      <lo l="0">
	    <wf matchingOrder="0" trialEnergy="-2.0" searchE="true"/>
	    <wf matchingOrder="1" trialEnergy="-2.0" searchE="true"/>
      </lo>
      <lo l="0">
	    <wf matchingOrder="0" trialEnergy="-15.0" searchE="true"/>
	    <wf matchingOrder="0" trialEnergy="-2.0" searchE="true"/>
      </lo>
      
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="7.37" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="7.37" searchE="false"/>
       </lo>
       
       
        <custom l="1" type="lapw" trialEnergy="-0.5" searchE="true"/>
        <lo l="1">
     	  <wf matchingOrder="0" trialEnergy="-11.5" searchE="true"/>
     	  <wf matchingOrder="1" trialEnergy="-11.5" searchE="true"/>
        </lo>
        <lo l="1">
          <wf matchingOrder="1" trialEnergy="-11.5" searchE="true"/>
          <wf matchingOrder="2" trialEnergy="-11.5" searchE="true"/>
        </lo>
        <lo l="1">
   	      <wf matchingOrder="0" trialEnergy="-0.5" searchE="true"/>
   	      <wf matchingOrder="1" trialEnergy="-0.5" searchE="true"/>
        </lo>
        <lo l="1">
   	      <wf matchingOrder="0" trialEnergy="-11.5" searchE="true"/>
   	      <wf matchingOrder="0" trialEnergy="-0.5" searchE="true"/>
        </lo>

      
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="8.41" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="8.41" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="22.47" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="22.47" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="41.2" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="41.2" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="64.28" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="64.28" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="91.57" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="91.57" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="122.98" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="122.98" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="158.45" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="158.45" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="197.93" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="197.93" searchE="false"/>
       </lo>
         
         
        <custom l="2" type="lapw" trialEnergy="0.33" searchE="true"/>
        <lo l="2">
          <wf matchingOrder="0" trialEnergy="-6.4" searchE="true"/>
          <wf matchingOrder="1" trialEnergy="-6.4" searchE="true"/>
        </lo>
        <lo l="2">
  	      <wf matchingOrder="1" trialEnergy="-6.4" searchE="true"/>
  	      <wf matchingOrder="2" trialEnergy="-6.4" searchE="true"/>
        </lo>
        <lo l="2">
  	      <wf matchingOrder="0" trialEnergy="0.33" searchE="true"/>
  	      <wf matchingOrder="1" trialEnergy="0.33" searchE="true"/>
        </lo>
        <lo l="2">
  	      <wf matchingOrder="0" trialEnergy="-6.4" searchE="true"/>
  	      <wf matchingOrder="0" trialEnergy="0.33" searchE="true"/>
        </lo>

       <lo l="2">
        <wf matchingOrder="0" trialEnergy="9.83" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="9.83" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="22.89" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="22.89" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="40.29" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="40.29" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="61.86" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="61.86" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="87.52" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="87.52" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="117.22" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="117.22" searchE="false"/>
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
        <wf matchingOrder="0" trialEnergy="10.18" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="10.18" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="21.43" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="21.43" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="36.93" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="36.93" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="56.58" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="56.58" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="80.34" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="80.34" searchE="false"/>
       </lo>
        
       <lo l="3">
        <wf matchingOrder="0" trialEnergy="108.14" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="108.14" searchE="false"/>
       </lo>

    </basis>
  </sp>
</spdb>
"""

# Used for all calculations in /users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8
o_basis = """<?xml version="1.0" encoding="UTF-8"?>
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
        <wf matchingOrder="0" trialEnergy="5.85" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="5.85" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="16.74" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="16.74" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="31.96" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="31.96" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="51.29" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="51.29" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="74.62" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="74.62" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="101.92" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="101.92" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="133.15" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="133.15" searchE="false"/>
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
        <wf matchingOrder="0" trialEnergy="5.34" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="5.34" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="14.62" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="14.62" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="28.03" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="28.03" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="45.48" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="45.48" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="66.92" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="66.92" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="92.3" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="92.3" searchE="false"/>
       </lo>
        
       <lo l="1">
        <wf matchingOrder="0" trialEnergy="121.6" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="121.6" searchE="false"/>
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
        <wf matchingOrder="0" trialEnergy="10.87" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="10.87" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="22.33" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="22.33" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="37.83" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="37.83" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="57.3" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="57.3" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="80.73" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="80.73" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="108.1" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="108.1" searchE="false"/>
       </lo>
        
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="139.37" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="139.37" searchE="false"/>
       </lo>
         

    </basis>
  </sp>
</spdb>
"""