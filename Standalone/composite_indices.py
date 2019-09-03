#!/usr/bin/env python3

import numpy as np

nx=4
ny=3
nz=2

#Converting two indices to single composite index
def ij_to_cmp(ix,iy,nx,ny, index=0):
    return iy + (ix*ny) + index

print("2 nested loops")

cnt=0
for ix in range(0,nx):
    for iy in range(0,ny):
        icell = ij_to_cmp(ix,iy,nx,ny)
        print(cnt, icell)
        cnt+=1


print("")

# Note, this indexing is stupid in python 
cnt=0
for ix in range(0,nx):
    for iy in range(0,ny):
        cnt+=1
        icell = ij_to_cmp(ix,iy,nx,ny, index=1)
        print(cnt, icell)


print("")
print("3 nested loops")

# Converting three indices to composite index
# composite index = inner inex + outer_index*Nelements_of_inner_index
def ijk_to_cmp(ix,iy,iz,nx,ny,nz):
    ixy = ij_to_cmp(ix,iy,nx,ny, index=0)
    return ixy +(ix*ny*nz) 

cnt=0
for ix in range(0,nx):
    for iy in range(0,ny):
        for iz in range(0,nz):
            icell = iz + (iy*nz) + (ix*ny*nz)
            alt_icell = ijk_to_cmp(ix,iy,iz,nx,ny,nz)
            print(cnt, icell, alt_icell)
            cnt+=1
    
