#!/usr/bin/env python3

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import os
import math
import numpy as np


class SlurmJob:
    #walltime=[hours,minutes,seconds] - turn this into a dictionary 
    def __init__(self,exe,nodes=1,job_name='test_run',queue='test',walltime='None',parallel='None',output='terminal',\
                     cores_used=None):
        self.nodes=nodes
        self.ppn=28
        self.cpus_per_task = 1 
        self.job_name=job_name
        if cores_used == None:
            self.np = self.nodes*self.ppn* self.cpus_per_task
        else:
            self.np = self.cores_used 

        SlurmJob.check_queuename(queue)
        self.queue=queue
        
        if walltime=='None':
           self.walltime=self.average_walltime() 
        else:
            SlurmJob.check_walltime(walltime)
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
        if self.queue=='test':
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
        names=['test','cpu_test','veryshort','cpu','hmem','serial','serial_verylong','head']
        if queue not in names:
            print('Queue name not recognised: ',queue)
            sys.exit('Code has stopped')



def generate_slurm_script(slurmjob):
    time_str = get_time_string(slurmjob.walltime)
    
    slurm_string='#!/bin/bash                               \n'+ \
                '#SBATCH --job-name='+slurmjob.job_name+   '\n'+ \
                '#SBATCH --partition='+slurmjob.queue+    '\n'+ \
                '#SBATCH --nodes='+str(slurmjob.nodes)+    '\n'+ \
                '#SBATCH --ntasks-per-node='+str(slurmjob.ppn)+             '\n'+ \
                '#SBATCH --cpus-per-task='+str(slurmjob.cpus_per_task) +    '\n'+ \
                '#SBATCH --time='+time_str[:-1]+'\n' \

    slurm_string += '#Set variables      \n'+ \
                'export OMP_NUM_THREADS=1'+'\n'+ \
                'EXE='+slurmjob.exe       +'\n'+ \
                'OUT='+slurmjob.output    +'\n'+ \
                'cd $SLURM_SUBMIT_DIR       \n'

    if slurmjob.parallel == 'hybrid':
        print('Write slurm submission string generator for hybrid calculations')
        sys.exit('script has stopped')
                      
    slurm_string +=  'mpiexec -np '+str(slurmjob.np)+' $EXE > $OUT'

    return slurm_string
        


#Format time string for SLURM secheduler 
def get_time_string(walltime):
    time_str=''
    for i in range(0,3):
        ws = str(walltime[i])
        if len(ws) <2:
            time_str += '0'+ws+':'
        else:
            time_str +=     ws+':'
    return time_str
