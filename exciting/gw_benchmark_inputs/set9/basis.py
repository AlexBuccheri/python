"""
LO cut-offs per channel for optimised basis, and ground state input settings
"""
import sys

from gw_benchmark_inputs.input_utils import restructure_energy_cutoffs

n_energies_per_channel = 6


def set_lo_channel_cutoffs(l_max: dict) -> list:
    """
    Define the LO energy cut-offs per LO channel
    Ranges from (6, 5) to (7, 6) for (Zr, O)

    :return: list energy_cutoffs: Energy cut-offs for each LO channel of Zr and O
    """

    # (6, 5)
    if l_max['zr'] == 6 and l_max['o'] == 5:
        energy_cutoffs = {'zr': {0: [80, 100, 120, 150, 180, 200],
                                 1: [80, 100, 120, 150, 180, 200],
                                 2: [80, 100, 120, 150, 180, 200],
                                 3: [80, 100, 120, 150, 180, 200],
                                 4: [80, 100, 120, 150, 180, 200],
                                 5: [80, 100, 120, 150, 180, 200],
                                 6: [80, 100, 120, 150, 180, 200]},

                          'o': {0: [80, 100, 120, 150, 180, 200],
                                1: [80, 100, 120, 150, 180, 200],
                                2: [80, 100, 120, 150, 180, 200],
                                3: [80, 100, 120, 150, 180, 200],
                                4: [80, 100, 120, 150, 180, 200],
                                5: [80, 100, 120, 150, 180, 200]}
                          }

    # (7, 6)
    elif l_max['zr'] == 7 and l_max['o'] == 6:
        energy_cutoffs = {'zr': {0: [80, 100, 120, 150, 180, 200],
                                 1: [80, 100, 120, 150, 180, 200],
                                 2: [80, 100, 120, 150, 180, 200],
                                 3: [80, 100, 120, 150, 180, 200],
                                 4: [80, 100, 120, 150, 180, 200],
                                 5: [80, 100, 120, 150, 180, 200],
                                 6: [80, 100, 120, 150, 180, 200],
                                 7: [80, 100, 120, 150, 180, 200]},

                          'o': {0: [80, 100, 120, 150, 180, 200],
                                1: [80, 100, 120, 150, 180, 200],
                                2: [80, 100, 120, 150, 180, 200],
                                3: [80, 100, 120, 150, 180, 200],
                                4: [80, 100, 120, 150, 180, 200],
                                5: [80, 100, 120, 150, 180, 200],
                                6: [80, 100, 120, 150, 180, 200]}
                          }

    else:
        sys.exit("L max pair for Zr and O not valid")

    return restructure_energy_cutoffs(n_energies_per_channel, energy_cutoffs)


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
      do="fromfile"
      rgkmax="8.0"
      ngridk="6 6 6"
      xctype="GGA_PBE_SOL"
      epsengy="1.e-6"
      deband="0.001"
      gmaxvr="20.0"
      >
   </groundstate>

  {GW_INPUT}

</input>
"""
