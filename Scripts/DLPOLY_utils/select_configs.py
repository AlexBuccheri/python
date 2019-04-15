#!/usr/bin/env python3

#Import libraries
import numpy as np
import subprocess
import random
import sys
import os 

# Alexander Buccheri April 2019. UoB
# Extract n configurations from a DLPOLY HISTORY file

#Return line numbers of grep matches
def line_numbers_of_string_match(str_to_match,fname):
    grep_str = "grep -n '"+str_to_match+"' "+fname+" | cut -f1 -d:"
    output   = subprocess.check_output(grep_str,shell=True).decode("utf-8").splitlines()
    line_numbers = [int(x) for x in output] 
    return line_numbers

#Read in HISTORY file
def return_trajectory(fname):
    fid = open(file = fname)
    lines = fid.readlines()
    fid.close()
    return np.asarray(lines)

def return_ith_config(config_id):
    i1 = line_numbers[config_id-1]-1
    i2 = line_numbers[config_id]-1
    return history[i1:i2]

def prepend_config_header(prefix_header,first_line):
    first_line = prefix_header + first_line
    first_line = ' '.join(first_line.split())
    return first_line

def all_configids_valid(config_ids, max_config_id):
    if max(config_ids) > max_config_id:
        return False
    else:
        return True 

def mkdir(dir_name, verbose=False):
    if(os.path.isdir(dir_name)==True):
        if(verbose == True): print("Directory already exists: ",dir_name)
        return False
    else:
       os.system("mkdir "+dir_name)
       if(verbose == True): print("Created directory: ",dir_name)
       return True 
    
# Extact configurations from DLPOLY's HISTORY file 
# Call as a function, instead of using as a standalone script
# Arguments
# fname:                         path_2_history
# output_dir:                    File written to output_dir/configi/CONFIG
# Pass number of random configs: nconfigs
# or specifiy configs:           config_ids
#
def extact_configs_from_history(fname,output_dir, nconfigs=None,config_ids=None):

    if nconfigs==None and config_ids==None:
        nconfigs = 5 
    
    fname += '/HISTORY'
    line_numbers = line_numbers_of_string_match('timestep',fname)
    history = return_trajectory(fname)
    max_config_id = len(line_numbers)
    print(max_config_id,'configurations in history file')

    if config_ids:
        if not all_configids_valid(config_ids, max_config_id):
            print('Asked for an ith configuration that exceeds \
                   the nth configuration.')
            sys.exit('Script has stopped')
    elif not config_ids:
        config_ids = []
        for x in range(nconfigs):
            config_ids.append(random.randint(1,len(line_numbers)+1))
            print('Extracting ',nconfigs,'randomly selected configurations:')
            print( config_ids )
            
    #Write new config files
    cnt = 0
    for config_id in config_ids:
        config = return_ith_config(config_id)
        config[0] = prepend_config_header(prefix_header,config[0])

        fid = open(output_dir+'/config'+str(cnt)+'/CONFIG', "w")
        fid.write(config[0]+'\n')
        for i in range(1,len(config)-1):
            fid.write(config[i])

        #Avoid ending file with trailing whiteline     
        fid.write(config[len(config)-1].rstrip())
        fid.close()

        cnt += 1


    
# ---------------------------------------
# Main Routine if calling script standalone
# ---------------------------------------
fname = '/Users/alexanderbuccheri/silicate/1000k/prod'
output_dir = '500k'
prefix_header = 'CONFIG generated from 300 k '
#If blank, code selects configs at random 
config_ids = []
nconfigs = 5



#Set up files 
if not os.path.isdir(output_dir):
    print("Top-level output directory does not exist: ",output_dir)
    sys.exit('Script has stopped')
    
fname += '/HISTORY'
line_numbers = line_numbers_of_string_match('timestep',fname)
history = return_trajectory(fname)
print(len(line_numbers),'configurations in history file')

#Select configurations from this trajectory at random 
if not config_ids:
    config_ids = []
    for x in range(nconfigs):
        config_ids.append(random.randint(1,len(line_numbers)+1))
    print('Extracting ',nconfigs,'randomly selected configurations:')
    print( config_ids )

#Write new config files
cnt = 0 
for config_id in config_ids:
    config = return_ith_config(config_id)
    config[0] = prepend_config_header(prefix_header,config[0])

    subdir = output_dir+'/config'+str(cnt)
    mkdir(subdir, verbose=True)
    fid = open(subdir+'/CONFIG', "w")
    
    fid.write(config[0]+'\n')
    for i in range(1,len(config)-1):
        fid.write(config[i])

    #Avoid ending file with trailing whiteline     
    fid.write(config[len(config)-1].rstrip())
    fid.close()
    cnt += 1 




