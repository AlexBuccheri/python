#!/usr/bin/env python3

from ase.io import read,write

# Use ASE to produce DLPOLY CONFIG file from CIF file
# Example: Several supercells for SiO2 (alpha-quartz)

#Silicate cell integer(x,y & z):natoms in supercell
basis = 3+6
cell_natoms={6:1944,7:3087,8:4608,9:6561}
a=read("cif_inputs/sio2.cif")

for ncell in cell_natoms:
    natoms = cell_natoms[ncell]
    print(ncell,natoms,int(basis*ncell**3.))
    write( str(ncell)+"CONFIG",a.repeat(ncell) )
    print('In FIELD file: ATOMS 9, NUMMOLS ',natoms/9)


# Note: a.repeat(6) produces a supercell with 1944 atoms
# ASE orders a unit cell '3 Si' then '6 O' = 9 atoms
# To be consistent in the FIELD file, ATOMS 9 and NUMMOLS = 1944/9 = 216

