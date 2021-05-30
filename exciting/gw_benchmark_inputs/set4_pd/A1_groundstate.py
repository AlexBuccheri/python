"""
These settings correspond to sufficiently converged ground state calculations.
What will differ between subsequent ground state runs are the basis files.
"""

# deband was required in order for the ground state calculation to converge:
# It is the initial step length used when searching for the band energy, which is used as the APW linearisation energy
# The default was larger, at 0.0025.
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

      <species speciesfile="Zr.xml" rmt="1.6">
         <atom coord="0.00 0.00 0.00"></atom>
      </species>

      <species speciesfile="O.xml" rmt="1.6">
         <atom coord="0.25 0.25 0.25"></atom>
         <atom coord="0.75 0.75 0.75"></atom>
      </species>
   </structure>

   <groundstate
      do="fromfile"
      rgkmax="8.0"
      ngridk="6 6 6"
      xctype="GGA_PBE_SOL"
      epsengy="1.e-7"
      deband="0.001"
      gmaxvr="20.0"
      >
   </groundstate>

  {GW_INPUT}

</input>
"""
