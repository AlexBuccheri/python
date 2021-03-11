converged_ground_state_input = \
"""
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
      do="skip"
      rgkmax="7.0"
      ngridk="6 6 6"
      xctype="GGA_PBE_SOL"
      epsengy="1.e-8"
      gmaxvr="20.0"
      >
   </groundstate>

  {GW_INPUT}

</input>
"""
