"""
Given a converged basis from set9/A1_zr_o.py, systematically reduce the
basis function LOs per channel, such that a minimal basis for converged
calculations can be found.
"""

# Zr MT radius = 2 and O MT radius = 1.6
#
# L_max=(6,5)
# 100 5941.336342663189

channel_cutoffs = {'zr': {'s': [],  # Removes 1, 2 and 3 highest LOs
                          'p': [],
                          'd': [],
                          'f': [],
                          'g': [],
                          '?': [],
                          '?': []},

                   'o': {'s': [],
                         'p': [],
                         'd': [],
                         'f': [],
                         'g': [],
                         '?': []}
                   }



def main():

    # Define the converged basis (6, 5) with some cut-off

    # Define cut-offs for every channel to remove 2, 4 and 6 LOs

    # Copy ground state calculation

    # Create new Zr and O species files

    # subdirectory: zr_s, zr_p, zr_..., o_s, o_p