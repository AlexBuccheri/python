#!/usr/bin/env python3
import sys
import numpy as np
import os

""" Collate DLPOLY output files utilised by statis.py 

Collates files used by the analysis tool statis.py into one folder
allowing for easy SCP'ing to one's local machine.
A lack of PyQT on BC3 means I can't do analysis on BC3
line.

Should be run in the job root 

"""

__author__    = "Alexander Buccheri <ab17369@bristol.ac.uk>"
__copyright__ = "Copyright 2019 A Buccheri"
__license__   = "GPL-3.0-only"


#Functions
def mkdir(dir_name, verbose=False):
    if(os.path.isdir(dir_name)==True):
        if(verbose == True): print("Directory already exists: ",dir_name)
        return False
    else:
       os.system("mkdir "+dir_name)
       if(verbose == True): print("Created directory: ",dir_name)
       return True 

def delete_dir(dir_name):
    if(os.path.isdir(dir_name)==True):
        os.system("rm -r "+dir_name)
        print("Deleted directory: ",dir_name)
    return


#Main Routine 
dT=100
dirs=np.arange(1700,1800+dT,dT)
gather=False
root='analysis'

if gather==True:
    if(os.path.isdir(root)==True):
        print('analysis directory already exists')
        sys.exit('Script has stopped')

    mkdir('analysis',True)
    print('Moving production files into analysis folder')
    
    for t in dirs:
        old_dir = str(t)+'k/prod/'
        new_dir = root+'/'+str(t)+'k/'
        mkdir(new_dir)
        print( 'mv '+old_dir+'STATIS '+new_dir+'STATIS' )
        os.system('mv '+old_dir+'STATIS '+new_dir+'STATIS') 
        os.system('mv '+old_dir+'RDFDAT '+new_dir+'RDFDAT') 
        os.system('cp statis.py '+new_dir+'statis.py') 
        
if gather==False:
    if(os.path.isdir(root)==False):
        print('No analysis directory')
        sys.exit('Code has stopped')
        
    print('Returning files back to run directory')
    for t in dirs:
        new_dir = str(t)+'k/prod/'
        old_dir = root+'/'+str(t)+'k/'
        mkdir(new_dir)
        os.system('mv '+old_dir+'STATIS '+new_dir+'STATIS') 
        os.system('mv '+old_dir+'RDFDAT '+new_dir+'RDFDAT')
    delete_dir(root)
