#!/usr/bin/env python3
#Alexander Buccheri. University of Bristol. Feb 2018

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import math
import numpy as np
import os
import matplotlib.pyplot as plt

#My module paths 
blue3='Linux-2.6.32-642.6.2.el6.x86_64-x86_64-with-redhat-6.4-Carbon'
AlexMac='Darwin-17.4.0-x86_64-i386-64bit'

if platform.platform() == blue3:
    sys.path.insert(0, '/panfs/panasas01/chem/ab17369/python_modules')
if  platform.platform() == AlexMac:
    sys.path.insert(0, '/Users/alexanderbuccheri/Python')

#My modules
from Modules.scheduler_py import pbs 





#-------------------------------------------------------------
#Job submission and processing for band structure calculation
#-------------------------------------------------------------

exe='$HOME/dftb/master_build/prog/dftb+/dftb+'
queue='head'
material='GaAs'

job_location1=material+'_bs/1.converged_charges'
job_location2=material+'_bs/2.bands'

#Charges
job1=pbs.PbsJob(exe,nodes=1,job_name=material+'_charges',queue=queue,parallel='openmp')
pbs.submit_job(job1,job_location1,queue)
os.system('cp '+job_location1+'/charges.bin '+job_location2+'/')

#Band structure
job2=pbs.PbsJob(exe,nodes=1,job_name=material+'_band',queue=queue,parallel='head')
pbs.submit_job(job2,job_location2,queue)

#Produce band structure and DOS plots 
#Need to chnage to correct directories to use
#os.system('dp_dos band.out dos_total.dat')
root = os.getcwd()
os.chdir(job_location2)
os.system('dp_bands band.out band')
os.chdir(root)
os.system('cp plot.p '+job_location2+'/')
