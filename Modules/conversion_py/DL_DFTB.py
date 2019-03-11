#!/usr/bin/env python3

#--------------------------------------------------------------------------
# Functions for converting/copying data from DLPOLY to DFTB+, and vice versa
#--------------------------------------------------------------------------

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import math
import numpy as np

#My libraries 
sys.path.insert(0, '/Users/alexanderbuccheri/Python')
from Modules.dftb_py import classes as dftb
from Modules.dlpoly_py import dl_classes as dl

#--------------------------------------------
# Switches leading dimension of 2D array
# where leading_dim = 's'mallest or 'l'argest 
#--------------------------------------------
def reshape2Darray(array,leading_dim):
    leading_dim=str(leading_dim).lower()[0]
    #Want leading dimension = largest of the two
    if leading_dim=='l' and len(array[0])> len(array):
            array=np.transpose(array)
    #Want leading dimension = smallest of the two
    if leading_dim=='s' and len(array)> len(array[0]):
            array=np.transpose(array)           
    return array


#--------------------------------------------------------
# Takes unique elements index from DFTB+ geometry object
# and uses it to construct a corresponding index list
# for all elements in the system 
#--------------------------------------------------------
def generate_elements_index(geo):
    elements_index=[]
    for ele in geo.elements:
        elements_index.append(geo.uni_elements_index[ele])
    return elements_index


#----------------------------------------------------------------------------------------------------
# Copies geometry data from DLPOLY config object to DFTB+ geometry object
# Inputs:  Filled 'config' object 
#          override_boundary      => DFTB can do a finite simulation,
#                                   despite using periodic boundaries in DLPOLY
# Output: Filled 'geo' object
#
# Information that's not copied/transferred between objects
# Number of lines per record. Not needed by geometry class
# atom_index. Seems superfluous. geometry class can easily construct this from 'uni_elements_index'
#----------------------------------------------------------------------------------------------------

def configObject_to_geoObject(config,override_boundary=False):

    #DFTB+ class expects dimensions of position[1:Natom,1:3], so enforce
    config.coord=reshape2Darray(config.coord,leading_dim='largest')

    #DFTB+ uses finite BC despite DLPOLY option
    if override_boundary == True:
        print('DLPOLY BC is ',dl.Constants.boundary_key[config.boundary_index-1],' but DFTB will use finite.')
        boundary_conditions='c'
        
    #DLPOLY BC maps to DFTB+ BC
    elif override_boundary == False: 
        #Finite
        if config.boundary_index==1:
            boundary_conditions='c'
        #Periodic cubic.
        elif config.boundary_index==2: 
            boundary_conditions='s'
        #Other boundary conditions 
        elif config.boundary_index>2:
            boundary_conditions='s'
            print('Note, no direct mapping between DLPOLY boundary condition and DFTB+ boundary condition')
            print('DLPOLY boundary ', config.boundary_index, 'set to ',boundary_conditions,' in DFTB+ geo object')
            
    #Should not need to be unit conversion here. Both inputs in Ang (although internal DFTB+ in Bohr)
    geo=dftb.Geometry(material=None,elements=config.atom_name,boundary_conditions=boundary_conditions,\
                      position=config.coord,lattice_vectors=config.lattice_vector, header=config.header)

    #Consistency Check
    check_consistency(geo,config)

    #Additional MD initialisaiton information
    if config.Nlines_per_record==2:  print('Velocity information from config object not assigned to geo object')
    if config.Nlines_per_record==3:  print('Velocity & force information from config object not assigned to geo object')
        
    return geo

# ---------------------------------------------------
# Convert fractional coordinates to Cartesian
# ---------------------------------------------------
def convert_fractional_2_cartesian(latvec,pos_frac):
    # Set up transformation matrix using period lattice vectors
    T = np.array([latvec[0,:], latvec[1,:], latvec[2,:]])
    T = np.transpose(T)
    #Do transformation on positions
    pos_cart = np.matmul(T, np.transpose(pos_frac))
    pos_cart = np.transpose(pos_cart)
    return pos_cart


#--------------------------------------------------------------------------
# Copies geometry data from DFTB+ config object to DLPOLY geometry object
# Inputs:  Filled 'geo' object 
# Output: Filled 'config' object       
#--------------------------------------------------------------------------
def geoObject_to_configObject(geo):
    
    #DFTB geo data doesn't contain info on initial values for velocities or forices, hence
    #Nlines_per_record = 1
    #DLPOLY uses large supercell with either finite or periodic boundaries.
    #Corresponding DFTB+ atomic positions for large systems will be stored in geo.position

    #Boundary matching
    if geo.boundary_conditions.lower()=='c':
        boundary_index=1    #finite
        position=geo.position
    if geo.boundary_conditions.lower()=='s':
        boundary_index=2    #Cubic periodic
        position = geo.position
    if geo.boundary_conditions.lower()=='f':
        print('DFTB+ basis/atomic position vectors are in fractional coorinates.')
        print('Converting to same unit as the lattice vectors')
        boundary_index = 2  # Cubic periodic
        position = convert_fractional_2_cartesian(geo.lattice_vectors, geo.position)

    config=dl.Config(header=geo.header,Nlines_per_record=1,boundary_index=boundary_index,Natoms=geo.Natoms,\
                     lattice_vector=geo.lattice_vectors,atom_name=geo.elements,coord=position,\
                     atom_index=generate_elements_index(geo))

    return config


def check_consistency(geo,config):
    if geo.Natoms != config.Natoms:
        print("Number of atoms in system disagrees between DLPOLY and DFTB+")
        print("DLPOLY Natoms:",config.Natoms,' DFTB+ Natoms:',geo.Natoms)
        sys.exit('DL_DFTB conversion has stopped')
