#!/usr/local/bin/python3

################################################################
# Convert .gin data to .xyz data 
# Hard-coded for pyrope
################################################################

#Import libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import sys


#----------------------------
# User options 
#----------------------------

#Format 'xyz' for VMD or'gen' for DFTB+
fformat='gen'

#Unit output 'Ang' or 'Red'
length='Ang'

#Number of headers in input file 
NumHeaders=1

#Output file for (C)luster or periodic (S)upercell .gen file  
boundary='S'

#Hard-coded lattice constant (in Ang)
al=11.456

#Hard-coded origin and lattice vectors (in Ang)
origin=np.zeros(shape=(3))   
lvector = np.array([[al, 0., 0.], [0., al, 0.], [0., 0.,al]])

#File name 
fname='position_data/pyrope.dat'


#----------------------------
# Main Routine 
#----------------------------

print ('Format of outputted data file: ',fformat)

if boundary=='C':
    print ('Finite boundary conditions: Cluster')
if boundary=='S':
    print ('Periodic boundary conditions: Supercell')

if length=='Ang':
    print ("Atomic positions output in Angstroms")
else:
    al=1. #Going to multiply atomic positions by al, hence initialise
    print('Atomic positions output in reduced coordinates')


#Read in data 
data = np.genfromtxt(fname, dtype=str,  delimiter='\n')

N=(len(data))-NumHeaders
element=[None]*N
etype=[None]*N
pos=np.zeros(shape=(3,N))
charge=np.empty(shape=(N))
occ=np.empty(shape=(N))
rad1=np.empty(shape=(N))


#Read data, skipping header
for i in range(NumHeaders,len(data)):
    j=i-NumHeaders
    #White spaces used as delimiters to split data (hence no argument in split())
    element[j],etype[j],x,y,z,charge[j],occ[j],rad1[j] = data[i].split()
    pos[:,j]=([x,y,z])
    pos[:,j]=pos[:,j]*al


#Construct dictionary of elements
#https://stackoverflow.com/questions/1024847/add-new-keys-to-a-dictionary
Natoms=len(element)    
cnt=1
element_index = {element[0]: cnt}
unq_element=([element[0]])

for i in range(1,Natoms):
    exists= (element[i] in element_index)
    #If element is not in dictionary, add it 
    if exists == False:
        cnt=cnt+1
        element_index.update({element[i]:cnt})
        unq_element.append(element[i])


#--------------------
#Output data
#--------------------
f=open(fname[:-3]+fformat, mode='wt', encoding='utf-8')

if fformat=='xyz':
    print (Natoms,file=f)
    print ('  ',file=f)
    for i in range(0,Natoms):
        #Format specifier {column,length type}, where none type is specified for string
        #* required for tuples and lists: https://stackoverflow.com/questions/7568627/using-python-string-formatting-with-lists
        print ('{0:2} {1:10f} {2:10f} {3:10f}'.format(element[i],*pos[:,i]),file=f)
    f.close()
    print('xyz data written to '+fname[:-3]+'xyz')



if fformat=='gen':
    #finite boundary conditions 
    if boundary=='C':
        print ('# Pyrope atomic positions ('+fformat+') for cubic unit cell.',file=f)
        print (' ',Natoms,' '+boundary+'       #',Natoms,' atoms, finite cluster',file=f)
        print (' '.join(unq_element), file=f)   #Print list without square brackets 
        for i in range(0,Natoms):
            print ('{0:3d} {1:2} {2:10f} {3:10f} {4:10f}'.format(i+1, element_index[element[i]], *pos[:,i]),file=f)
        f.close()
        
    #periodic boundary conditions 
    if boundary=='S':
        print ('# Pyrope atomic positions ('+fformat+') for cubic unit cell.',file=f)
        print (' ',Natoms,' '+boundary+'       #',Natoms,' atoms, supercell',file=f)
        print (' '.join(unq_element), file=f)  
        #Atoms in supercell
        for i in range(0,Natoms):
            print ('{0:3d} {1:2} {2:10f} {3:10f} {4:10f}'.format(i+1, element_index[element[i]], *pos[:,i]),file=f)
        #Origin and cell vectors
        print ('{0:10f} {1:10f} {2:10f}'.format(*origin[:]),file=f)
        for i in range(0,len(lvector)):
            print ('{0:10f} {1:10f} {2:10f}'.format(*lvector[i,:]),file=f)
        f.close()

