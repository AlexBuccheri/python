#!/usr/bin/env python3

#Libraries
import sys
import numpy as np
import os
import subprocess
import shutil
import re
#import pathlib  #Not always available 

#My libraries
sys.path.insert(0, '/panfs/panasas01/chem/ab17369/python_modules')
from Modules.dlpoly_py import dl_classes as dl

# ---------------------------------
# Directory manipulation functions
# ----------------------------------
#Require import pathlib
#def mkdir(dir_name, verbose=False):
#    dir_ = pathlib.Path(dir_name)
#    if(dir_.exists()):
#        if(verbose == True): print("Directory already exists: ",dir_name)
#        return False
#    else:
#       os.system("mkdir "+dir_name)
#       if(verbose == True): print("Created directory: ",dir_name)
#       return True 

#def delete_dir(dir_name):
#    dir_ = pathlib.Path(dir_name)
#    if(dir_.exists()):
#        os.system("rm -r "+dir_name)
#        print("Deleted directory: ",dir_name)
#    return


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

#Split directory path into list containing subdirectories 
def get_subdirectories(wholedir):
    subdir = []
    subdir = wholedir.split('/')
    if( wholedir[0] =='/' ):
        subdir = wholedir[1:].split('/')
    if( wholedir[len(wholedir)-1] =='/' ):
        subdir = wholedir[:len(wholedir)-1].split('/')
    if( wholedir[0] =='/' and wholedir[len(wholedir)-1]  =='/' ):
        subdir = wholedir[1:len(wholedir)-1].split('/')
    return subdir 

#Should modify to take input 
def mpi_run(MPI_np,EXE,output=None):
   #Same line as a new shell is initialised on each os.system call
    if( output == None):
        os.system("export OMP_NUM_THREADS=1 && mpirun -np "+str(MPI_np)+" "+EXE)
    else:
        os.system("export OMP_NUM_THREADS=1 && mpirun -np "+str(MPI_np)+" "+EXE+" > "+output)
    return

#Run DLPOLY calculation - works for any MPI 
def run_dlpoly(run_dir,exe,np):
    root_dir = os.getcwd()
    os.chdir(run_dir)
    mpi_run(np,exe)
    os.chdir(root_dir)
    return 


# ----------------------------------
#  DLPOLY Classes 
# ----------------------------------
class ControlFname:
    def __init__(self,old='CONTROL',new='CONTROL'):
        self.old=old
        self.new=new 

class ConfigFname: 
    def __init__(self,old='CONFIG',new='CONFIG'):
        self.old=old
        self.new=new 
        

#Reads backwards through OUTPUT, until the header of the final
#output summary is found, then takes corresponding volume
#output[7] in DLPOLY V4.09.  output[6] in DLPOLY v4.10
def extract_average_volume(file_name):
    if not os.path.isfile(file_name):
        print('Output file was not found to read volume from')
        print('Exit function and allow script to continue')
        return 0.
    output=[]
    for line in reversed(open(file_name).readlines()):
        output.append(line)   
        if 'step' in line:
            break 
    output=list(reversed(output))
    V_avg = output[6].split()[1] #Formatting always consistent
    return float(V_avg)


def extract_volume_rms(fname):
    if fname[-6:].lower() != 'output':
        print('Filename must end in OUTPUT')
        sys.exit('Script has stopped')
    if not os.path.isfile(file_name):
        print('Output file was not found to read volume RMS from')
        print('Exit function and allow script to continue')
        return 0.  
    grep_str = "grep -a -A 2 'r.m.s.' "+fname 
    output = subprocess.check_output(grep_str,shell=True).decode("utf-8")
    output = output.splitlines()
    vol_rms = float(output[2].split()[1])
    return vol_rms


def output_VvsT_to_file(file_name,V,T,Vrms=None):
    header = '# Temperature (k),  Volume (ang^3)'
    if Vrms != None:
        header += ', Volume R.M.S.'

    if not os.path.isfile(file_name):
        f=open(file_name, mode='w', encoding='utf-8')
        print(header, file=f)  
    else:
        f=open(file_name, mode='a', encoding='utf-8')

    if Vrms == None: print(T,V,      file=f)
    if Vrms != None: print(T,V,Vrms, file=f)
    f.close()


# '+(-)' to count up (down) in T
def get_sign(Tmin,Tmax):
    if Tmin <= Tmax:
        sign = 1
    if Tmin > Tmax:
        sign = -1
    return int(sign) 

# ----------------------------------
# Generate DLPOLY inputs 
# ----------------------------------
#Set up equilibration calculation 
def setup_equilibration_calculation(source_dir,run_dir,control,config,field):
    #Create directory for equilibration run
    if( os.path.isdir(run_dir)==True ):
        print('Run directory for equilibrium job already exists: ',run_dir)
        sys.exit('Script has stopped')
    subdir = get_subdirectories(run_dir)
    directory=subdir[0]
    mkdir(directory, verbose=True)
    for i in range(1,len(subdir)):
        directory+='/'+subdir[i]
        mkdir(directory, verbose=False)
        
    #Directory containing reference FIELD and CONFIG files   
    if( os.path.isdir(source_dir)==False ):
        print('Source directory containing CONFIG and FIELD does not exist: ',source_dir)
        sys.exit('Script has stopped')
    #Copy from prior CONFIG or copy input REVCON to CONFIG
    shutil.copyfile(source_dir+'/'+config.old , run_dir+'/'+config.new)
    shutil.copyfile(source_dir+'/'+field , run_dir+'/'+field)

    #Generate CONTROL file string and write to file 
    control_str = dl.Control_data_to_string(control)
    fid= open(run_dir+'/'+control.name, "w")
    fid.write(control_str)
    fid.close()
    return


#Essentially same function as above with different warning...
def setup_production_calculation(source_dir,run_dir,control,config,field):
    #Create directory for production run
    if( os.path.isdir(run_dir)==True ):
        print('Run directory for production job already exists: ',run_dir)
        sys.exit('Script has stopped')
    subdir = get_subdirectories(run_dir)
    directory=subdir[0]
    mkdir(directory, verbose=False)
    for i in range(1,len(subdir)):
        directory+='/'+subdir[i]
        mkdir(directory, verbose=False)

    #Directory containing reference FIELD and CONFIG files     
    if( os.path.isdir(source_dir)==False ):
        print('Source directory containing CONFIG and FIELD does not exist: ',source_dir)
        sys.exit('Script has stopped')
        
    #Copy from prior CONFIG or copy input REVCON to CONFIG
    shutil.copyfile(source_dir+'/'+config.old , run_dir+'/'+config.new)
    shutil.copyfile(source_dir+'/'+field , run_dir+'/'+field)

    #Generate CONTROL file string and write to file 
    control_str = dl.Control_data_to_string(control)
    fid= open(run_dir+'/'+control.name, "w")
    fid.write(control_str)
    fid.close()
    return

# ------------------
# Run a calculation
# ------------------
def equilibration_calculation(T, source_dir,run_dir,equi_control,equi_config,field,exe,np):
    equi_control.temperature = T
    setup_equilibration_calculation(source_dir,run_dir,equi_control,equi_config,field)
    run_dlpoly(run_dir,exe,np)
    os.system("rm "+run_dir+'/'+"HISTORY")
    source_dir = run_dir 
    return source_dir

def production_calculation(T,source_dir,run_dir,prod_control,prod_config,field,exe,np):
    prod_control.temperature = T
    setup_production_calculation(source_dir,run_dir,prod_control,prod_config,field)
    run_dlpoly(run_dir,exe,np)
    os.system("tar -czf "+run_dir+"/HISTORY.tar.gz "+run_dir+'/'+'HISTORY')
    os.system("rm "+run_dir+'/'+"HISTORY")
    source_dir = run_dir 
    return source_dir

def compute_equi2(T,dT,NT):
    if( int(T/dT)%NT == 0): 
        compute_equi2 = True
    else:
        compute_equi2 = False
    return compute_equi2

def get_slurm_nproc(fname):
    output_bytes = subprocess.check_output("grep 'ppn=' "+fname ,shell=True)
    output = str( output_bytes.decode("utf-8") )
    equals = [i.start() for i in re.finditer('=', output)]
    hashes = [i.start() for i in re.finditer('#', output)]
    #Last equals sign: equals[-1], and last hash sign: hashes[-1]
    i1 = equals[-1] +1
    i2 = len(output)-1
    if( hashes[-1] > equals[-1] ): i2 = hashes[-1] 
    ppn = int( output[i1:i2] )
    return ppn
