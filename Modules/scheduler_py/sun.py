#!/usr/bin/env python3

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import os


#For sun scheduler
class SunJob:
        def __init__(self,shell_type='tcsh',exe='~/TB_BSE/test.exe',nodes=1,ppn=1,job_name='TBtest',\
                          queue='parallel.q',walltime='None',input_name='input.fdf', output_name='terminal'):
                          
                self.shell_type=shell_type
                self.exe=exe
                self.nodes=nodes
                self.ppn=SunJob.check_ppn(ppn)
                self.job_name=job_name
                SunJob.check_queuename(queue)
                self.queue=queue     
                if walltime=='None':
                        self.walltime=[1,0,0] 
                else:
                        SunJob.check_walltime(walltime)
                self.walltime= walltime
                self.input=input_name
                self.output=output_name

        def check_ppn(ppn):
            if ppn>8:
                print('Cores per node limited to 8.')
                ppn=8
            return ppn
            
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
        
        def check_queuename(queue):
               names=['parallel.q','head']
               if queue not in names:
                  print('Queue name not recognised: ',queue)
                  sys.exit('Code has stopped')
            


def generate_sun_script(sunjob):
        time_str=get_time_string(sunjob.walltime)
        
        sun_string=('#!/bin/'+sunjob.shell_type  +'\n'+ \
               '#$ -S /bin/'+sunjob.shell_type   +'\n'+ \
               '#$ -l qname='+sunjob.queue       +'\n'+ \
               '#$ -l h_rt='+time_str[:-1]       +'\n'+ \
               '#$ -N '+sunjob.job_name          +'\n'+ \
               '#$ -pe orte 8 '                  +'\n'                        
               '#$ -j y  '                       +'\n'                       
               '#$ -R y  '                       +'\n'                                                         
               '#$ -cwd  '                       +'\n'       
               '#$ -V    '                       +'\n'    
               '#$ -o  diagnostic.out  '         +'\n'
               '#$ -e  diagnostic.err  \n\n')          
                  
        sun_string += ('echo -----------START OF MY OUTPUT-------------'   +'\n'
                       'echo Job $JOB_ID started on `date`             '   +'\n'
                       'echo  \n\n' )

        sun_string += ('#Set variables'                  +'\n'
                       'set MPIRUN  = mpirun'            +'\n'
                       'set EXE ='+sunjob.exe            +'\n'+ \
                       'set INPUT ='+sunjob.input        +'\n'+ \
                       'set OUT ='+sunjob.output         +'\n\n')


        sun_string += ('#Copy files required at runtime to working dir'         +'\n'
                       'cp -f $SGE_O_WORKDIR/$INPUT $TMPDIR'                    +'\n'
                       'cp -rf /home/buccheri/TB_BSE/parameter_files/ $TMPDIR'  +'\n'
                       '#Switch to working dir'                                 +'\n'
                       'cd $TMPDIR'                                             +'\n'
                       'echo $TMPDIR'  +'\n\n')

        sun_string += '$MPIRUN -np '+str(sunjob.ppn)+' $EXE $INPUT > $OUT'  +'\n\n'

        sun_string += ('#Copy outputs from runtime dir to output dir'  +'\n'
                       'mkdir $SGE_O_WORKDIR/$JOB_ID'                  +'\n'
                       'rm -r $TMPDIR/parameter_files'                 +'\n'
                       'cp -rf $TMPDIR/* $SGE_O_WORKDIR/$JOB_ID'       +'\n'
                       'cd $SGE_O_WORKDIR'                             +'\n\n'
                       'echo'                                          +'\n'
                       'echo Job $JOB_ID completed on `date`'          +'\n'
                       'echo -----------END OF MY OUTPUT-------------' +'\n'
                       '#Tidy the outputfile into the output directory'+'\n'                                   
                       'mv  $SGE_STDOUT_PATH $SGE_O_WORKDIR/$JOB_ID      \n')                                     
   
        return sun_string


#Format time string for Sun secheduler 
def get_time_string(walltime):
    time_str=''
    for i in range(0,3):
        ws = str(walltime[i])
        if len(ws) <2:
            time_str += '0'+ws+':'
        else:
            time_str +=     ws+':'
    return time_str
