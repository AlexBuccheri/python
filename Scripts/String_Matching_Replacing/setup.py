#!/usr/bin/env python

#----------------------------------------------------------------
# Alexander Buccheri. University of Bristol.
# Generate input files for my TB_BSE code
#----------------------------------------------------------------
from __future__ import print_function
import re
import sys
import os
import platform
from fileio     import String_replacement,modify_variables,add_directory
from sun        import SunJob,generate_sun_script

#--------------------
#Main Routine
#--------------------
fin='input.fdf'
fin2='run.csh'


#Submission script
qname='parallel.q'
np=4
wt=[1,0,0]
exe='~/TB_BSE/test.exe'

#Radius and dipole approximation         
Radius=[25,30]
NTT=[1,2,3]

keyword =['Radius','NumTransTerms'  ]
terminating_char=['#','#']
whitespace=['     ','     ']


#set up files
for rad in Radius:
    for ntt in NTT:
	
        fout='../SZ/'+str(0.2*rad)[0]+'nm/input'+str(ntt)+'.fdf'
        fout2='../SZ/'+str(0.2*rad)[0]+'nm/run'+str(ntt)+'.csh'
	print('Writing files:',fout,fout2)
        add_directory(fout)
        
        #Input script
        replacement= [rad , ntt]        
        obj = String_replacement(before=keyword,after=terminating_char,\
                                replace=replacement,padding=whitespace)
        modify_variables(fin,fout,obj)
        
        #run script
        testname='spec_'+str(0.2*rad)[0]+'nm_ntt'+str(ntt)
        sun = SunJob(shell_type='tcsh',exe=exe,nodes=1,ppn=np,job_name=testname,queue=qname,\
                     walltime=wt,input_name='input'+str(ntt)+'.fdf')
        sun_string=generate_sun_script(sun)
        fid= open(fout2, "w+")
        fid.write(sun_string)
        fid.close()


