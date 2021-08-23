"""
LO cut-offs per channel for optimised basis, and ground state input settings
"""
import sys

from gw_benchmark_inputs.input_utils import restructure_energy_cutoffs


def set_lo_channel_cutoffs(l_max: dict) -> list:
    """
    Define the LO energy cut-offs per LO channel
    Ranges from (4, 3) to (7, 6) for (Zr, O)

    I assume 150 will be enough but may need to go higher

    :return: list energy_cutoffs: Energy cut-offs for each LO channel of Zr and O
    """

    n_energies_per_channel = 4

    # (4, 3)
    if l_max['zr'] == 4 and l_max['o'] == 3:
        energy_cutoffs = {'zr': {0: [80, 100, 120, 150],
                                 1: [80, 100, 120, 150],
                                 2: [80, 100, 120, 150],
                                 3: [80, 100, 120, 150],
                                 4: [80, 100, 120, 150]},

                          'o': {0: [80, 100, 120, 150],
                                1: [80, 100, 120, 150],
                                2: [80, 100, 120, 150],
                                3: [80, 100, 120, 150]}
                          }

    # (5, 4)
    elif l_max['zr'] == 5 and l_max['o'] == 4:
        energy_cutoffs = {'zr': {0: [80, 100, 120, 150],
                                 1: [80, 100, 120, 150],
                                 2: [80, 100, 120, 150],
                                 3: [80, 100, 120, 150],
                                 4: [80, 100, 120, 150],
                                 5: [80, 100, 120, 150]},

                          'o': {0: [80, 100, 120, 150],
                                1: [80, 100, 120, 150],
                                2: [80, 100, 120, 150],
                                3: [80, 100, 120, 150],
                                4: [80, 100, 120, 150]}
                          }

    # (6, 5)
    elif l_max['zr'] == 6 and l_max['o'] == 5:
        energy_cutoffs = {'zr': {0: [80, 100, 120, 150],
                                 1: [80, 100, 120, 150],
                                 2: [80, 100, 120, 150],
                                 3: [80, 100, 120, 150],
                                 4: [80, 100, 120, 150],
                                 5: [80, 100, 120, 150],
                                 6: [80, 100, 120, 150]},

                          'o': {0: [80, 100, 120, 150],
                                1: [80, 100, 120, 150],
                                2: [80, 100, 120, 150],
                                3: [80, 100, 120, 150],
                                4: [80, 100, 120, 150],
                                5: [80, 100, 120, 150]}
                          }

    # (7, 6)
    elif l_max['zr'] == 7 and l_max['o'] == 6:
        energy_cutoffs = {'zr': {0: [80, 100, 120, 150],
                                 1: [80, 100, 120, 150],
                                 2: [80, 100, 120, 150],
                                 3: [80, 100, 120, 150],
                                 4: [80, 100, 120, 150],
                                 5: [80, 100, 120, 150],
                                 6: [80, 100, 120, 150],
                                 7: [80, 100, 120, 150]},

                          'o': {0: [80, 100, 120, 150],
                                1: [80, 100, 120, 150],
                                2: [80, 100, 120, 150],
                                3: [80, 100, 120, 150],
                                4: [80, 100, 120, 150],
                                5: [80, 100, 120, 150],
                                6: [80, 100, 120, 150]}
                          }

    else:
        sys.exit("L max pair for Zr and O not valid")

    return restructure_energy_cutoffs(n_energies_per_channel, energy_cutoffs)


# Note, I reduced the MT radius of Zr from 2.0 to 1.6, but left O
# Based on the GS calculation input used for A1_more_APW/set7
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
