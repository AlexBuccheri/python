#!/usr/bin/env python3

#Import libraries
import numpy as np
import subprocess
import random

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


#Main Routine 
fname = '/Users/alexanderbuccheri/silicate/1000k/prod'
prefix_header = 'CONFIG generated from 300 k '
#If blank, code selects configs at random 
config_ids = [1]
nconfigs = 5

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
for config_id in config_ids:
    config = return_ith_config(config_id)
    config[0] = prepend_config_header(prefix_header,config[0])

    fid = open('CONFIG_'+str(config_id), "w")
    fid.write(config[0]+'\n')
    for i in range(1,len(config)-1):
        fid.write(config[i])

    #Avoid ending file with trailing whiteline     
    fid.write(config[len(config)-1].rstrip())
    fid.close()




