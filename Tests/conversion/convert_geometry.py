#!/usr/local/bin/python3

###############################################################
# Alex Buccheri UoB. May-July 2018.
# Test Script 
#   Read in DFTB+ and convert to CONFIG
#   Read in CONFIG and convert to DFTB+
#   Write both out and confirm they're the same (Need to add)
###############################################################

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import numpy as np
import os

#My libraries 
sys.path.insert(0, '/Users/alexanderbuccheri/Python')
from Modules.dftb_py import classes as dftb                # DFTB Classes
from Modules.dftb_py import hsd                            # .hsd file generator 
from Modules.dftb_py import conversion as dftb_conversion  # Data <-> string routines   
from Modules.dlpoly_py import dl_classes as dl             # DLPOLY Classes
from Modules.conversion_py import DL_DFTB                  # Geometry class-mapping functions 
 

#Conversion choice 
convert_DL_to_DFTB=False
convert_DFTB_to_DL=True


#----------------------------------
#DL POLY config to DFTB+ .gen
#----------------------------------
if convert_DL_to_DFTB:
    print('Converting DLPOLY CONFIG file to DFTB+ .gen file')
    #Parses CONFIG data in string format and fills CONFIG class with the information.
    config_string = np.genfromtxt('inputs/CONFIG', dtype=str,  delimiter='\n')
    #Object
    config=dl.Config()
    #Use string to fill config
    dl.Config_string_to_data(config_string,config)
    #Use config data to fill gen data
    geo=DL_DFTB.configObject_to_geoObject(config,override_boundary=True)
    print(vars(geo))
    #Convert gen data to string Need to write this in
    gen_string = hsd.hsd_gen_string_cluster(geo)
    print(gen_string)
    #Write out string to the .gen file
    fid= open('outputs/graphene.gen', "w")
    fid.write(gen_string[:-1])#Last line is blank so don't write it
    fid.close()


# ----------------------------------
# DFTB+ .gen to DL POLY config
# ----------------------------------
if convert_DFTB_to_DL:
    print('Converting DFTB+ .gen file to DLPOLY CONFIG file')
    # Parses .gen data in string format and fills DFTB class with the information.
    gen_string = np.genfromtxt('inputs/struct.gen', dtype=str, delimiter='\n')
    #print(gen_string)
    #Fill gen object from string (DFTB+ .gen uses Ang by default)
    geo = dftb_conversion.Gen_string_to_data(gen_string)
    #print(vars(geo))
    #Convert geo object to config object (retaining boundary used in DFTB+ file)
    config= DL_DFTB.geoObject_to_configObject(geo)
    #print(vars(config))
    #Convert config data to string and write to file
    config.header='# Header: Converted from .gen file'
    config_string=dl.config_data_to_string(config)
    #print(config_string)
    fid = open('outputs/CONFIG_Ti', "w")
    fid.write(config_string[:-1])  # Last line is blank so don't write it
    fid.close()


#----------------------------------
#DFTB+ .gen to DL POLY config
#----------------------------------
#Rather than read the DFTB+ gen data from file, take it from the
#object I'v just constructed 
#config2=DL_DFTB.geoObject_to_configObject(geo)
#print(vars(config2))
#Write CONFIG data to file 
















