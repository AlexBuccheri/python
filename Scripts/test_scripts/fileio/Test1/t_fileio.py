#!/usr/bin/env python3

#----------------------------------------------------------------
# Alexander Buccheri. University of Bristol.
# Test 1 for fileio class. Specifically for use with my TB code
#----------------------------------------------------------------
import re
import sys
import os
import platform

#Paths to my modules (alternatively put directory in $PATH)
#Blue crystal varies according to specific log-in node
blue3='Linux-2.6.32-642.6.2.el6.x86_64-x86_64-with-redhat-6.4-Carbon'
AlexMac='Darwin-17.4.0-x86_64-i386-64bit'
Landau='Linux-2.6.32-504.16.2.el6.x86_64-x86_64-with-centos-6.6-Final'

if platform.platform()[0:12] == blue3[0:12]:
    sys.path.insert(0, '/panfs/panasas01/chem/ab17369/python_modules')
if platform.platform()[0:6] == AlexMac[0:6]:
    sys.path.insert(0, '/Users/alexanderbuccheri/Python')
if platform.platform() == Landau:
    sys.path.insert(0, '/home/buccheri/Python')

    
#My modules 
from Modules.fileio_py.fileio     import String_replacement,modify_variables,add_directory


#--------------------
#Main Routine
#--------------------
fin='input.fdf'

#Radius and dipole approximation         
Radius=[10.]
NTT=[1,2,3]

keyword =['Radius','NumTransTerms'  ]
terminating_char=['#','#']
whitespace=['     ','     ']

for rad in Radius:
    for ntt in NTT:
        fout='outputs/'+str(0.2*rad,)[0]+'nm/input'+str(ntt)+'.fdf'
        add_directory(fout)
        
        replacement= [rad , ntt]
        
        obj = String_replacement(before=keyword,after=terminating_char,\
                                replace=replacement,padding=whitespace)
        
        modify_variables(fin,fout,obj)

