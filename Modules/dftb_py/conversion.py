#!/usr/local/bin/python3

import numpy as np
import sys
from . import classes as dftb_class

# ---------------------------------------------------
# Take a gen string and convert to DFTB Geo object
# ---------------------------------------------------
def Gen_string_to_data(string, header=False):

    #Optional header line
    hoff=0
    if header==True:
        header=string[0]
        hoff=1

    #Line 2 of string if header present, else line 1 of string
    #Number of atoms and boundary conditions
    Natoms,boundary_conditions = string[0+hoff].split()
    Natoms=int(Natoms)

    #Line 2. List of unique elements, in order
    uni_elements = string[1+hoff].split()

    #Atomic positions in supercell
    pos=np.zeros(shape=(Natoms,3))
    elements=[]

    istart=2+hoff
    ia=0
    for i in range(istart,istart+Natoms):
        idummy, el_index, x,y,z = string[i].split()
        element_string=uni_elements[int(el_index)-1]
        elements.append(element_string)
        pos[ia,:] = [float(x),float(y),float(z)]
        ia+=1

    #Origin and lattice vectors
    index=istart+Natoms
    #x,y,z=string[index].split()
    #origin = [float(x),float(y),float(z)]

    latvec=np.zeros(shape=(3,3))
    ia=0
    for i in range(index+1,index+4):
        x, y, z = string[i].split()
        latvec[ia,:]= [float(x),float(y),float(z)]
        ia+=1

    #Assign local data to DFTB geometry object (Note, position and basis essentially degenerate)
    geo_data = dftb_class.Geometry(material='dummy_material',elements=elements,
                             boundary_conditions=boundary_conditions,
                             position=pos, al=0. ,basis_vectors=pos,
                             lattice_vectors=latvec,header=header)

    return geo_data







