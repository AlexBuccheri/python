from typing import List
import numpy as np

from o_basis import o_basis
from zr_basis import zr_basis

# def parse_lorecommendations(file_name):
#     fid = open(file='lorecommendations.txt', mode='r')
#     lines = fid.readlines()
#     fid.close()
#
#     # File has good structure
#     l_min = 0
#     l_max = 3
#     n_min = 0
#     n_max = 20
#     n_atoms = 3
#     energy
#
#     # Scrap header
#     lines = lines[2:]
#
#     # Skip species = 1 and l = 0
#     i = 2
#     basis_per_atom = []
#     for i_atom in range(0, 1):
#         for i_l in range(l_min, l_max + 1):
#             lo_l = []
#             for i_n in range(n_min, n_max + 1):
#                 lo_l.append({n:trial_energy})
#                 print(lines[i])
#                 i+=1
#             # Single line break + l index
#             i+=2
#         # species
#         i+=1
#         basis_per_atom.append(lo_l)
#
# parse_lorecommendations('lorecommendations.txt')

# Functions to inject



# 6 calculations
# Include Step 20 (ha?) i.e. -ve-0, >0-20, >20-40, >40-60, >60-80, >80-100

energy_cutoffs = [100.]   #[20., 40., 60., 80., 100]
max_matching_order = 1

for energy_cutoff in energy_cutoffs:
    print("ZR. Energy cutoff: " + str(energy_cutoff) + " max_matching_order:" + str(max_matching_order))
    print(zr_basis(energy_cutoff, max_matching_order))
    print("OXYGEN. Energy cutoff: " + str(energy_cutoff) + " max_matching_order:" + str(max_matching_order))
    print(o_basis(energy_cutoff, max_matching_order))






