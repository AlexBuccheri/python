#!/usr/local/bin/python3

#Libraries 
import sys
import numpy as np
import os

#My libraries 
sys.path.insert(0, '/Users/alexanderbuccheri/Python')
from Modules.dftb_py import classes as dftb

#---------------------------------------------------------------------------
#Test 3 cases of DFTB+ geometry class
#Note, must specify the argument names in the class if not using all of them
#---------------------------------------------------------------------------

#-------------------------
#Periodic bulk
#-------------------------
structure='cubic'
basis='prim'
elements = np.array(['S','Zn']) #Order should be consistent with anion at (0,0,0) and cation at (0.25,0.25,0.25)
material='ZnS'
boundary_conditions='S'
atomic_structure=dftb.Lattice(structure=structure,basis=basis,material=material)
periodic_bulk=dftb.Geometry(material,elements,boundary_conditions,\
                            al=atomic_structure.al,basis_vectors=atomic_structure.basis_vector,\
                            lattice_vectors=atomic_structure.lattice_vector,header=atomic_structure.header)
print('Periodic Bulk')
print(vars(periodic_bulk))
print(' ')


#-------------------------
#Periodic supercell
#-------------------------
structure='cubic'
basis='prim'
atomic_structure=dftb.Lattice(structure=structure,basis=basis,material=material)

elements = np.array(['H','C','C','C','C','C','C','C','C','C','C','C','C','C','C','H','H','H']) 
material='Graphene'
boundary_conditions='S'
position=[
[0.000000,   0.000000,   0.000000],
[4.322529,   0.000000,   0.000000],
[8.645058,   0.000000,   0.000000],
[12.967587,   0.000000,   0.000000],
[17.290116,   0.000000,   0.000000],
[0.500000,  2.881686,   0.866025],
[4.822529,   2.881686,   0.866025],
[9.145058,   2.881686,   0.866025],
[13.467587,   2.881686,   0.866025],
[17.790116,   2.881686,   0.866025],
[0.000000,   0.000000,   2.495613],
[4.322529,   0.000000,   2.495613],
[8.645058,   0.000000,   2.495613],
[12.967587,   0.000000,   2.495613],
[17.290116,   0.000000,   2.495613],
[0.500000,   2.881686,   3.361639],
[4.822529,   2.881686,   3.361639],
[9.145058,   2.881686,   3.361639]]


PS=dftb.Geometry(material,elements,boundary_conditions,position=position,\
                 al=atomic_structure.al, basis_vectors=atomic_structure.basis_vector,\
                 lattice_vectors=atomic_structure.lattice_vector,header='DFTB+ graphene supercell')

print('Periodic supercell')
print(vars(PS))
#Elements index can easily be constructed from uni_elements_index dictionary
elements_index=[]
for ele in elements:
    elements_index.append(PS.uni_elements_index[ele])
print('elements_index',elements_index)
print(' ')

#-------------------------
#Finite system
#-------------------------
elements = np.array(['H','C','C','C','C','C','C','C','C','C','C','C','C','C','C','H','H','H']) 
material='Graphene'
boundary_conditions='C'
position=[
[0.000000,   0.000000,   0.000000],
[4.322529,   0.000000,   0.000000],
[8.645058,   0.000000,   0.000000],
[12.967587,   0.000000,   0.000000],
[17.290116,   0.000000,   0.000000],
[0.500000,  2.881686,   0.866025],
[4.822529,   2.881686,   0.866025],
[9.145058,   2.881686,   0.866025],
[13.467587,   2.881686,   0.866025],
[17.790116,   2.881686,   0.866025],
[0.000000,   0.000000,   2.495613],
[4.322529,   0.000000,   2.495613],
[8.645058,   0.000000,   2.495613],
[12.967587,   0.000000,   2.495613],
[17.290116,   0.000000,   2.495613],
[0.500000,   2.881686,   3.361639],
[4.822529,   2.881686,   3.361639],
[9.145058,   2.881686,   3.361639]]
finite=dftb.Geometry(material,elements,boundary_conditions,position=position,header='DFTB+ graphene finite system')
print('Finite')
print(vars(finite))
#Elements index can easily be constructed from uni_elements_index dictionary
elements_index=[]
for ele in elements:
    elements_index.append(finite.uni_elements_index[ele])
print('elements_index',elements_index)
