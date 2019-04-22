#!/usr/bin/env python3

# --------------------------------------------------------
# DLPOLY V vs T Automation 
# --------------------------------------------------------
#Libraries
import sys
import numpy as np
import os
import subprocess
import shutil
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
def extract_average_volume(file_name):
    output=[]
    for line in reversed(open(file_name).readlines()):
        output.append(line)   
        if 'step' in line:
            break 
    output=list(reversed(output))
    V_avg = output[7].split()[1] #Formatting always consistent 
    return float(V_avg)


def output_VvsT_to_file(file_name,V,T):
    if not os.path.isfile(file_name):
        f=open(file_name, mode='w', encoding='utf-8')
        print('# Volume (Ang^3) vs Temperature (K)', file=f)
        print(V,T, file=f)
        f.close()
    else:
        f=open(file_name, mode='a', encoding='utf-8')
        print(V,T, file=f)
        f.close()

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
    source_dir = run_dir 
    return source_dir

def production_calculation(T,source_dir,run_dir,prod_control,prod_config,field,exe,np):
    prod_control.temperature = T
    setup_production_calculation(source_dir,run_dir,prod_control,prod_config,field)
    run_dlpoly(run_dir,exe,np)
    source_dir = run_dir 
    return source_dir

def compute_equi2(T,dT,NT):
    if( int(T/dT)%NT == 0): 
        compute_equi2 = True
    else:
        compute_equi2 = False
    return compute_equi2


# ------------------------------------------
# Main Routine
# Must come here as Tmin used to set control
# ------------------------------------------
#exe = '/panfs/panasas01/chem/ab17369/codes/clean/dl-poly-master/build-intel16u2-mpi/bin/DLPOLY.Z'
exe = '/Users/alexanderbuccheri/Codes/development/dl-poly-4.10/build-gcc7-mpi/bin/DLPOLY.Z'
np = 1

source_dir='700k/equil'
Tmin = 700
Tmax = 1000
dT = 100
NT = 1


# --------------------------    
# DLPOLY input settings 
# --------------------------

#Initial DLPOLY config and field files
#equi_config  = ConfigFname(old='CONFIG', new='CONFIG') - If starting from input files that have not been ran
equi_config  = ConfigFname(old='REVCON', new='CONFIG')
equi2_config = ConfigFname(old='REVCON', new='CONFIG')
prod_config  = ConfigFname(old='REVCON', new='CONFIG')

field='FIELD'

#Initialise control data for equilibrium calculation 
equi_ensemble = dl.Ensemble(name='npt', etype='hoover', thermostat_relaxation=1.0, barostat_relaxation=0.5)
equi_trajectory = dl.Trajectory(tstart=10000,tinterval=1000,data_level=0)
equi_control = dl.Control(header='Quartz',temperature=Tmin, pressure=0.001, ensemble=equi_ensemble,    \
                     steps=500000, equilibration=500000, scale=7, regauss=13, printout=1000, stack=1000, \
                     stats=1000, vdw_direct=True, rdf=1000, print_rdf=True, trajectory=equi_trajectory, \
                     ewald_precision=1.0e-6, timestep=0.001, rpad=0.4, cutoff=6.0, cap=1.0e3,           \
                     job_time=9.0e5, close_time=2.0e1)

#Initialise control data for 2nd equilibrium calculation, where velocities aren't rescaled 
equi2_control = dl.Control(header='Quartz',temperature=Tmin, pressure=0.001, ensemble=equi_ensemble,    \
                     steps=500000, equilibration=500000, printout=1000, stack=1000, \
                     stats=1000, vdw_direct=True, dump=1000000, \
                     ewald_precision=1.0e-6, timestep=0.001, rpad=0.4, cutoff=6.0, cap=1.0e3,           \
                     job_time=9.0e5, close_time=2.0e1)


#Initialise control data for production calculation 
prod_ensemble = dl.Ensemble(name='npt', etype='hoover', thermostat_relaxation=1.0, barostat_relaxation=0.5)
prod_trajectory = dl.Trajectory(tstart=0,tinterval=1000,data_level=0)
prod_control = dl.Control(header='Quartz',temperature=Tmin, pressure=0.001, ensemble=prod_ensemble,    \
                     steps=1000000, equilibration=0,  printout=1000, stack=1000, \
                     stats=1000, vdw_direct=True, rdf=1000, print_rdf=True, trajectory=prod_trajectory, \
                     dump=1000000, ewald_precision=1.0e-6, timestep=0.001, rpad=0.4, cutoff=6.0, cap=1.0e3,           \
                     job_time=9.0e5, close_time=2.0e1)




#If calculation needs to start at equil2 or production, do so here 
run_dir='700k/equil2'  
source_dir = equilibration_calculation(Tmin,source_dir,run_dir,equi2_control,equi2_config,field,exe,np)
run_dir='700k/prod'  
source_dir = production_calculation(Tmin,source_dir,run_dir,prod_control,prod_config,field,exe,np)
Tmin = Tmin + dT


#DLPOLY calculations for different T
for T in range(Tmin,Tmax+dT,dT):
    print('Temperature:',T)

    #Equilibrium calculation
    run_dir = str(T)+'k/equil'
    source_dir = equilibration_calculation(T, source_dir,run_dir,equi_control,equi_config,field,exe,np)

    #Check equilibrium holds without rescaling velocities, every NT temperature steps
    if( compute_equi2(T,dT,NT) == True):
        print('Check equilibrium holds without scaling at T:',T)
        run_dir = str(T)+'k/equil2'
        source_dir = equilibration_calculation(T, source_dir,run_dir,equi2_control,equi2_config,field,exe,np)
    
    #Production calculation
    run_dir = str(T)+'k/prod'
    source_dir = production_calculation(T,source_dir,run_dir,prod_control,prod_config,field,exe,np)
    V_avg = extract_average_volume(run_dir+'/OUTPUT')
    output_VvsT_to_file('VvsT.dat',V_avg,T)

 
    

