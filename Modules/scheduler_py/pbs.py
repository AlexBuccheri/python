#!/usr/bin/env python3

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import os
import math
import numpy as np


def submit_job(job,job_location,where):  

    #Run on head node
    if where=='head':
        bash_string = generate_head_node_string(job)
        run_script = 'run_local.csh'
        fid= open(job_location+'/'+run_script, "w")
        fid.write(bash_string)
        fid.close()
        cwd = os.getcwd()
        os.chdir(job_location)
        os.system('chmod +x '+run_script)
        os.system('./'+run_script)
        os.chdir(cwd)    
        
    #Submission to queue
    if where!='head':
        #Generate submission script file 
        pbs_string = generate_pbs_script(job)
        run_script = 'run.csh'
        fid= open(job_location+'/'+run_script, "w")
        fid.write(pbs_string)
        fid.close()
        #Change from current directory to job_location, submit, then change back
        cwd = os.getcwd()
        os.chdir(job_location)
        os.system('qsub '+run_script)
        os.chdir(cwd)

        
class PbsJob:
    #walltime=[hours,minutes,seconds] - turn this into a dictionary 
    def __init__(self,exe,nodes=1,job_name='test_run',queue='testq',walltime='None',parallel='None',output='terminal'):
        self.nodes=nodes
        self.ppn=16
        self.job_name=job_name

        PbsJob.check_queuename(queue)
        self.queue=queue
        
        if walltime=='None':
           self.walltime=self.average_walltime() 
        else:
            PbsJob.check_walltime(walltime)
            self.walltime= walltime

        self.exe=exe
        self.output=output+'.out'
        self.parallel=parallel


    def check_walltime(walltime):
        if walltime[0] > 120:
            print('Number of hours should not exceed 120')
            print('Current setting: ',walltime[0])
            sys.exit('Code has stopped')
        if walltime[1] > 59:
            print('Number of minutes should not exceed 59')
            print('Current setting: ',walltime[1])
            sys.exit('Code has stopped')
        if walltime[2] > 59:
            print('Number of seconds should not exceed 59')
            print('Current setting: ',walltime[2])
            sys.exit('Code has stopped')
            
    def average_walltime(self):
        if self.queue=='testq':
            walltime=[0,10,0] 
        elif self.queue=='head':
            #Dummy wall time (as not in queue)
            walltime=[0,1,0]      
        else:
            print('Haven''t written an average walltime for queue: ',self.queue)
            print('Therefore using one hour')
            walltime=[1,0,0]
        return walltime
                    
    def check_queuename(queue):
        names=['testq','veryshort','short','medium','long','head']
        if queue not in names:
            print('Queue name not recognised: ',queue)
            sys.exit('Code has stopped')



def generate_pbs_script(pbsjob):
    time_str = get_time_string(pbsjob.walltime)
    
    pbs_string='#!/bin/bash                \n'+ \
                '#PBS -V                   \n'+ \
                '#PBS -l nodes='+str(pbsjob.nodes)+':ppn='+str(pbsjob.ppn)+ '\n'+ \
                '#PBS -q '+pbsjob.queue+   '\n'+ \
                '#PBS -l walltime='+time_str[:-1]+'\n' \
                '#PBS -N '+pbsjob.job_name+ '\n'

    pbs_string += '#Set variables      \n'+ \
                'EXE='+pbsjob.exe    +'\n'+ \
                'OUT='+pbsjob.output +'\n'+ \
                '#Change to submission dir on compute node  \n'+ \
                'cd $PBS_O_WORKDIR     \n'

    if pbsjob.parallel!=None:
        if pbsjob.parallel == 'openmp':
            pbs_string += '# set the number of OpenMP threads to match the resource request:   \n' +\
                          'NUMPROC=`wc -l < $PBS_NODEFILE` \n' +\
                           'export OMP_NUM_THREADS=$NUMPROC \n'

    pbs_string +=  '#Run        \n'+ \
                   '$EXE > $OUT \n'

    
    return pbs_string
        


#Format time string for PBS secheduler 
def get_time_string(walltime):
    time_str=''
    for i in range(0,3):
        ws = str(walltime[i])
        if len(ws) <2:
            time_str += '0'+ws+':'
        else:
            time_str +=     ws+':'
    return time_str


#For running short serial jobs on the head node
def generate_head_node_string(pbsjob):
    head_node_string='#!/bin/bash                \n'+ \
               'EXE='+pbsjob.exe         +'\n'+ \
               'OUT='+pbsjob.output      +'\n'+ \
               'export OMP_NUM_THREADS=1   \n'+ \
               '$EXE > $OUT '
    return head_node_string
    

        




