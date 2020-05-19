import json
import subprocess
import sys

solvent = "benzene"

# Structure from Minnesota solvation database (2012)
MNSol_xyz = """0044met H4C1O1 methanol m062x_mg3s_geom

0 1
  8   -0.0462410000   -0.7520490000    0.0000000000
  6   -0.0462410000    0.6609810000    0.0000000000
  1    0.8580790000   -1.0671410000    0.0000000000
  1   -1.0866750000    0.9767460000    0.0000000000
  1    0.4379840000    1.0704500000    0.8899200000
  1    0.4379840000    1.0704500000   -0.8899200000"""

atom_symbol_dict = {8: "O", 6: "C", 1: "H"}

structure = []
for line in MNSol_xyz.split('\n')[3:]:
    atomic_number, x_coord, y_coord, z_coord = line.split()
    structure.append([atom_symbol_dict[int(atomic_number)],
                     float(x_coord),
                     float(y_coord),
                     float(z_coord)])

print(structure)


input_str = """
gas := xtb(
  structure( xyz = {structure:s} )
)
sol := xtb(
  structure( xyz = {structure:s} )
  solvation( model = gbsa solvent = {solvent:s} )
)
""".format(solvent=solvent, structure=str(structure))
input_str.replace('\n', ' ')
#entos_command = [sys.argv[1], '--mute', '--json-results', '-s', input_str]

print(input_str)