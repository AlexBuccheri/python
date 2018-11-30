#!/usr/bin/env python3

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import numpy as np
import os

#My libraries 
sys.path.insert(0, '/Users/alexanderbuccheri/Python')
from Modules.dlpoly_py import dl_classes as dl

#----------------------------------
#DL POLY class example.
#----------------------------------
#Parses CONFIG data in string format and fills CONFIG class with the information. 
string = np.genfromtxt('lib/CONFIG', dtype=str,  delimiter='\n')
data=dl.Config()
dl.Config_string_to_data(string,data)
