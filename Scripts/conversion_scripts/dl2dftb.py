#!/usr/local/bin/python3

################################################################
# Convert DL POLY CONFIG file to DFTB+ structure input file
# and vice-versa
################################################################

#Import libraries
import math
import numpy as np
import sys
import fileinput


#---------------------------------------------------
#User settings
#---------------------------------------------------

#Specify file names in script (1) or at command line (2)
specify_fname=1

#Use consistent boundary conditions or specify different 'C' or 'S"
dftb_boundary='C'



#---------------------------------------------------
#Specify input and output file names in script 
#---------------------------------------------------
if specify_fname==1:
    print('File names specified in script')
    input_fname='CONFIG'
    output_fname='graphene_finite_struct.gen'
    
    if input_fname.find("CONFIG")==0:
        print('Converting DL POLY ',input_fname,' to ',output_fname)
        if input_fname.find(".gen")==0:
            print('Converting DFTB+ structure file ',input_fname,' to ',output_fname)


#---------------------------------------------------
#Accept an input and output name at command line 
#---------------------------------------------------

if specify_fname==2:
    #No files are passed at command line
    if len(sys.argv)==1:
        print('Script requires a file input at command line')
        sys.exit('Script has stopped')

    #If files are passed at the command line 
    if len(sys.argv) > 1:
  
        #Only input file is passed at command line
        if sys.argv[-2] == sys.argv[0]:
            input_fname = sys.argv[-1]
            #File is DL
            if input_fname.find("CONFIG")==0:
                output_fname='dftb_structure.gen'
                print('Converting DL POLY ',input_fname,' to ',output_fname)
                #File is DFTB+
            if input_fname.find(".gen")==0:
                output_fname='CONFIG'
                print('Converting DFTB+ structure file ',input_fname,' to ',output_fname)
            else:
                print('File format is not recognised')
                sys.exit('Script has stopped')
        
        #Input and output names are passed at command line
        if sys.argv[-2]!=sys.argv[0] and sys.argv[-1]!=sys.argv[0]:
            input_fname=sys.argv[-2]
            output_fname=sys.argv[-1]
        if input_fname.find("CONFIG")==0:
            print('Converting DL POLY ',input_fname,' to ',output_fname)
        if input_fname.find(".gen")==0:
            print('Converting DFTB+ structure file ',input_fname,' to ',output_fname)
        if input_fname.find("CONFIG")!=0 and input_fname.find(".gen")!=0:
           print('Input file format is not recognised')
           sys.exit('Script has stopped')
    

        

#Based on the file extensions, detect whether to convert DL CONFIG to DFTB+ or vice versa.
#-----------------------------------------
#Convert DL CONFIG to DFTB+ .gen
#-----------------------------------------
if input_fname.find("CONFIG")==0:
    #Assumes consistent formatting for CONFIG file 

    #If the input is DL's CONFIG format 
    #Want to read this from command line 
    input_fname='CONFIG'
    data = np.genfromtxt(input_fname, dtype=str,  delimiter='\n')

    #CONFIG Expected Format
    #Line 1 = Header
    header = data[0]
    #Line 2 = Num Lines per atomic record,  Boundary conditions,  Natoms
    Nlines_per_record, boundary_index, Natoms = data[1].split()
    
    #Convert to integers
    atomic_record=np.array(['coorindates', 'coordinates & velocities', 'coordinates, velocities & forces'])
    Nlines_per_record=int(Nlines_per_record)+1
    #Boundary key
    boundary_index=int(boundary_index)
    boundary_key=np.array(['finite','cubic','orthorhombic','parallel-piped','Not valid integer','Not valid integer','xy but not z'])
    #Total number of atoms in system 
    Natoms=int(Natoms)

    #Inform user
    print('DL POLY CONFIG file contains: ')
    print( Natoms,' atoms')
    print('Boundary conditions:',boundary_key[boundary_index])
    print('Each atomic record is comprised of ',atomic_record[Nlines_per_record-1])


    #Line 3 = Lattice vectors
    print('Lattice vectors')
    a=np.zeros(shape=(3,3))
    for i in range(2,5):
        j=i-2
        a[j,0:3] = data[i].split()  #No arg => whitespace delimiter
        print(a[j,0:3])


    #Initial position, velocity and force per atom
    offset=5
    k=offset
    atom_type=[None]*Natoms
    r_v_f=np.zeros(shape=(Natoms,Nlines_per_record,3))

    #Atomic records for Natoms, each comprised of 'Nlines_per_record+1' lines 
    for ia in range(0,Natoms):
        atom_type[ia],dummy_index = data[k].split()
        #print (atom_type[ia], dummy_index)
        k=k+1
        for j in range(0,Nlines_per_record):
            #print ('record line',data[k])
            r_v_f[ia,j,:]=data[k].split()
            #print(r_v_f[ia,j,:])
            k=k+1


    #-----------------------------------------
    #Output data to dftb .gen structure format
    #-----------------------------------------
    #DFTB and DL boundaries differ
    if dftb_boundary=='S' and boundary_key[boundary_index]=='finite':
        print('Note, DFTB and DL boundary conditions differ')
    if dftb_boundary=='C' and boundary_key[boundary_index]!='finite':
        print('Note, DFTB and DL boundary conditions differ')

    #Assume the origin is (0,0,0)
    origin=np.array([0.,0.,0.])

    #Construct dictionary of elements
    Natoms=len(atom_type)    
    cnt=1
    #Initalise 1st element of each
    atom_index = {atom_type[0]: cnt}
    unq_atom_index=([atom_type[0]])

    for i in range(1,Natoms):
        exists= (atom_type[i] in atom_index)
        #If element/species is not in dictionary, add it 
        if exists == False:
            cnt=cnt+1
            atom_index.update({atom_type[i]:cnt})
            unq_atom_index.append(atom_type[i])

    #Open file to write to 
    fid=open(output_fname, mode='wt', encoding='utf-8')


    #.gen file for finite boundary conditions - Don't need unit cell vectors 
    if dftb_boundary=='C':
        print ('# Header line.',file=fid)
        print (' ',Natoms,' '+dftb_boundary+'       #',Natoms,' atoms, finite cluster',file=fid)
        print (' '.join(unq_atom_index), file=fid)   #Print list without square brackets 
        for ia in range(0,Natoms):
            print ('{0:3d} {1:2} {2:10f} {3:10f} {4:10f}'.format(ia+1, atom_index[atom_type[ia]], *r_v_f[ia,0,:]),file=fid)
        fid.close()


    #.gen file for periodic boundary conditions 
    if dftb_boundary=='S':
        print ('# Header line',file=fid)
        print (' ',Natoms,' '+dftb_boundary+'       #',Natoms,' atoms, supercell',file=fid)
        print (' '.join(unq_atom_index), file=fid)  
        #Atoms in supercell
        for ia in range(0,Natoms):
            print ('{0:3d} {1:2} {2:10f} {3:10f} {4:10f}'.format(ia+1, atom_index[atom_type[ia]], *r_v_f[ia,0,:]),file=fid)
        #Origin and cell vectors
        print ('{0:10f} {1:10f} {2:10f}'.format(*origin[:]),file=fid)
        for i in range(0,len(a)):
            print ('{0:10f} {1:10f} {2:10f}'.format(*a[i,:]),file=fid)
        fid.close()


#-----------------------------------------
#Convert DFTB+ .gen to DL CONFIG 
#-----------------------------------------
if input_fname.find(".gen")==0:
    print('.gen to DL CONFIG needs writing')
