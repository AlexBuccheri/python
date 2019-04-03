#!/usr/bin/env python3

#Import libraries
import sys
import math
import numpy as np
import json

#----------------------------------------------------------------------------------------
# Convert DFTB+ hsd settings, structure input and force output into a JSON file
#----------------------------------------------------------------------------------------

#Headers must start with '#'
def get_header_length(fname):
    skip_rows = 0 
    fid = open(file = fname)
    lines = fid.readlines()
    if '#' not in lines[0]:
        fid.close()
        return skip_rows
    else:
        for line in lines:
            if '#' not in line: break
            skip_rows+=1
    fid.close()
    return skip_rows 

#Generate species index map of the form:
#{0:'Null, 1: 'C', 2: 'H', 3: 'A', ..., Nspecies:'X'}
def get_species_map(fname,Nheader_rows):
    fid = open(fname)
    for i, line in enumerate(fid):
        #Fixed format file 
        if i == Nheader_rows+1:
            elements = line.split()
            break
    fid.close()
    species_map = {0:'Null'}
    for i in range(0,len(elements)):
        species_map[i+1] = elements[i]
    return species_map
    
#Extract from forces header 
def get_md_step(fname):
    fid = open(fname)
    for i, line in enumerate(fid):
        md_step = line.split()[-1:]
        break
    fid.close()
    #md_step = [int(i) for i in md_step]  or list(map(int,md_step))
    return int( md_step[0] )

def get_boundary(fname,Nheader_rows):
    fid = open(fname)
    for i, line in enumerate(fid):
      #Fixed format file 
        if i == Nheader_rows:
            boundary = line.split()[1]
            break
    fid.close()
    return str(boundary)


#Returns dictionary of DFTB data with the keys:
# Keys = ['md_step', 'natoms', 'boundary', \
#         'latvec', 'hsd', 'atom', 'species_order']
def convert_data_to_dict(forces_file, gen_file, hsd_file):
    #Read forces
    Nrows = get_header_length(forces_file) 
    forces = np.loadtxt(forces_file,skiprows=Nrows)

    #Read gen file
    Nrows = get_header_length(gen_file)
    species_map = get_species_map(gen_file,Nrows)
    data = np.genfromtxt(gen_file, dtype='str', skip_header=Nrows+2, delimiter="\n")
    boundary = get_boundary(gen_file,Nrows)
    
    latvec=[]
    species_index=[]
    position=[]
    xyz=[]

    #Lattice vectors
    if(boundary.lower() =='s'):
        latvec_data = data[-3:]
        for line in latvec_data:
            xyz = [float(i) for i in line.split()]
            latvec.append( xyz )
        data = data[:-4]

    #Positions 
    for line in data:
        species_index.append( int(line.split()[1]) )
        xyz = [float(i) for i in line.split()[2:]]
        position.append( xyz )

    Natoms = len(species_index)

    #Fill DTFB dictionary
    dftb = {}
    dftb['md_step'] = get_md_step(forces_file) 
    dftb['natoms'] = int(Natoms)
    dftb['boundary'] = boundary
    if(dftb['boundary'].lower() == 's'): dftb['latvec'] = latvec
    dftb['hsd'] = np.genfromtxt(hsd_file, dtype='str', delimiter="\n").tolist()

    dftb['atom'] = []
    for ia in range(0,Natoms):
        dftb['atom'].append( {
            'species': species_map[species_index[ia]],
            'species_index':species_index[ia],
            'position':position[ia][:],
            'force':forces[ia].tolist()
            } )

    dftb['species_order'] = []
    for key,species in species_map.items():
        dftb['species_order'].append(species)

    return dftb

#Could make this more generic by looping over each key/value and printing     
def print_DFTB_dict(dftb):
    #keys=[]
    #for key in dftb_:
    #    keys.append(key)

    md_step  = dftb['md_step']
    Natoms   = dftb['natoms']
    boundary = dftb['boundary']
    latvec   = dftb['latvec']
    hsd_str  = dftb['hsd']
    atom     = dftb['atom']
    species_order = dftb['species_order']

    print('MD step: ' , md_step)
    print('Natoms: '  , Natoms)
    print('Boundary: ', boundary)
    print('Lattice vectors (ang): ',latvec)
    print('Species order: '        ,species_order)

    for ia in range(0,Natoms):
        print(ia, atom[ia]['species'], atom[ia]['species_index'])
        print('Position (ang):', atom[ia]['position'])
        print('Force (ha/bohr):', atom[ia]['force']) 
        print('')

    print('dftb_in.hsd: ')
    for line in hsd_str:
        print(line)

    print(' ')
    return 
    


#--------------------------
# Main Routine
#--------------------------
print('For a given calculation: Positions = inputs and forces = ouputs')

forces_file = 'inputs/forces_0.dat'
gen_file    = 'inputs/structure_0.gen'
hsd_file    = 'inputs/dftb_in.hsd'
json_file   = 'dftb_0_json.txt'

dftb = convert_data_to_dict(forces_file, gen_file, hsd_file)

with open(json_file, 'w') as fid:  
    json.dump(dftb, fid, indent=4)

with open(json_file) as fid:  
    dftb_0 = json.load(fid)
 
#print_DFTB_dict(dftb_0)    
