#!/usr/bin/env python3

#Import libraries
import math
import numpy as np
import sys



class Constants:
     boundary_key=np.array(['finite','cubic','orthorhombic','parallel-piped','Not valid integer',\
                            'Not valid integer','xy but not z'])
     atomic_record=np.array(['coorindates', 'coordinates & velocities', 'coordinates, velocities & forces'])

                            
#Data found in DLPOLY CONFIG file 
class Config:
    def __init__(self,header=None,Nlines_per_record=1,boundary_index=1,Natoms=0, \
                      lattice_vector=None,atom_name=None,coord=None,velo=None,force=None,atom_index=None):
                                                  #DL POLY Variable names
       self.header=header        
       self.Nlines_per_record=Nlines_per_record
       self.boundary_index=boundary_index
       self.Natoms=Natoms
       self.lattice_vector=lattice_vector  #Stored row-wise
       self.atom_name=atom_name
       self.coord=coord     #xxx, yyy and zzz
       self.velo=velo       #vxx, vyy and vzz 
       self.force=force     #fxx, fyy and fzz
       self.atom_index=atom_index

       #Does it make sense to initialise shapes?
       if lattice_vector is None:  self.lattice_vector=np.zeros(shape=(3,3))   
       if atom_name is None:       self.atom_name=[]         
       if coord is None:           self.coord=np.zeros(shape=(3,1))
       if atom_index is None:      self.atom_index=np.zeros(shape=(1))  
 
       # print('Boundary conditions:',Constants.boundary_key[boundary_index])
       # print('Each atomic record is comprised of ',Constants.atomic_record[Nlines_per_record-1])


    
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


# Convert DLPOLY CONFIG data into a string            
def config_data_to_string(config):
    config_string=''
    padding='          '

    if isinstance(config.header, str): config_string=config.header+'\n'

    config_string += padding+str(config.Nlines_per_record)+'    '+\
                     str(config.boundary_index)+'      '+\
                     str(config.Natoms)+'\n'

    for ia in range(0,3):
        vec = config.lattice_vector[ia, :]
        lat_str = np.array2string(vec, separator='    ', formatter={'float_kind': lambda vec: "%.8E" % vec})
        config_string += padding + lat_str[1:-1]+'\n'


    if config.Natoms==0: config.Natoms=len(config.coord)

    if config.Nlines_per_record != 1:
        print('Have not written CONFIG to string routine to output when'
              'Nlines_per_record exceeds 1')
        sys.exit('Script has stopped')

    if config.Nlines_per_record==1:
        for ia in range(0,config.Natoms):
            config_string+=config.atom_name[ia]+padding[:-len(config.atom_name[ia])-1]+str(ia+1)+'\n'
            pos = config.coord[ia,:]
            # Using this rather than 'join' to specify the formatting TURN INTO FUNCTION
            pos_str=np.array2string(pos, separator='    ', formatter={'float_kind': lambda pos: "%.8E" % pos})
            config_string += padding+pos_str[1:-1]+'\n'



    return config_string
