#!/usr/bin/env python3

#Import libraries
import math
import numpy as np
import sys


#Data found in DLPOLY CONFIG file 
class Config:
    def __init__(self,header=None,Nlines_per_record=1,boundary_index=1,Natoms=0, \
                      lattice_vector=None,atom_name=None,coord=None,velo=None,force=None,atom_index=None):
                                                  #DL POLY Variable names
       self.header=header        
       self.Nlines_per_record=Nlines_per_record
       self.boundary_index=boundary_index
       self.Natoms=Natoms

       if lattice_vector == None:
           self.lattice_vector=np.zeros(shape=(3,3))   #Stored row-wise
       else:
           self.lattice_vector=lattice_vector
       if atom_name == None:
           self.atom_name=[]
       else:
           self.atom_name=atom_name
       if coord == None:                               #xxx, yyy and zzz
           self.coord=np.zeros(shape=(3,1))
       else:
           self.coord=coord
       if velo == None:                                #vxx, vyy and vzz 
           self.velo=None
       else:
           self.velo=velo
       if force == None:                               #fxx, fyy and fzz   
           self.force=None
       else:
           self.force=force
       if atom_index == None:
           self.atom_index=np.zeros(shape=(1))
       else:
           self.atom_index=atom_index
    
    atomic_record=np.array(['coorindates', 'coordinates & velocities', 'coordinates, velocities & forces'])
    boundary_key=np.array(['finite','cubic','orthorhombic','parallel-piped','Not valid integer','Not valid integer','xy but not z'])
    # print('Boundary conditions:',boundary_key[boundary_index])
    # print('Each atomic record is comprised of ',atomic_record[Nlines_per_record-1])


    
#Take a CONFIG file parsed as a single string and assign data to DL POLY class 
def Config_string_to_data(string, data):

    data.header = string[0]
    Nlines_per_record, boundary_index, Natoms = string[1].split()
    
    data.Nlines_per_record=int(Nlines_per_record)+1
    data.boundary_index=int(boundary_index)
    data.Natoms=int(Natoms)
    #Allocate correct array size for coordinates 
    if data.coord.shape[1] < data.Natoms:
        data.coord.resize((3,data.Natoms))

    #Lattice vectors     
    for i in range(2,5):
        j=i-2
        data.lattice_vector[j,0:3] = string[i].split()  #No arg => whitespace delimiter

    #Atomic records for Natoms, each comprised of 'Nlines_per_record+1' lines
    offset=5
    k=offset
    
    if data.Nlines_per_record == 1:
        for ia in range(0,data.Natoms):
            atom_name,dummy_index = string[k].split()
            data.atom_name.append(atom_name)
            k=k+1
            data.coord[:,ia] = string[k].split()
            k=k+data.Nlines_per_record

    if data.Nlines_per_record == 2:
        print('Needs testing')
        if data.velo.shape(1) < data.Natoms: data.velo.resize((3,data.Natoms))
            
        for ia in range(0,data.Natoms):
            atom_name,dummy_index = string[k].split()
            data.atom_name.append(atom_name)
            k=k+1
            data.coord[:,ia] = string[k].split()
            data.velo[:,ia]  = string[k+1].split()
            k=k+data.Nlines_per_record

    if data.Nlines_per_record == 3:
        print('Needs testing')
        if data.velo.shape(1) < data.Natoms: data.velo.resize((3,data.Natoms))
        if data.force.shape(1) < data.Natoms: data.force.resize((3,data.Natoms))
            
        for ia in range(0,data.Natoms):
            atom_name,dummy_index = string[k].split()
            data.atom_name.append(atom_name)
            k=k+1
            data.coord[:,ia] = string[k].split()
            data.velo[:,ia]  = string[k+1].split()
            data.force[:,ia] = string[k+2].split()
            k=k+data.Nlines_per_record


#Write me 
#def config_data_to_string(config):
#    return config_str
