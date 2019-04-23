#!/usr/bin/env python3

# --------------------------------------------------------
# DLPOLY: V vs T Automation 
# --------------------------------------------------------
#Libraries
import sys
import numpy 
import os

#My libraries
sys.path.insert(0, '/panfs/panasas01/chem/ab17369/python_modules')
from Modules.dlpoly_py import dl_classes as dl
from Scripts.DLPOLY_workflows.functions import get_slurm_nproc,ConfigFname,equilibration_calculation,\
     production_calculation,extract_average_volume, output_VvsT_to_file,compute_equi2,\
     get_sign,extract_volume_rms 

#Script-specific functions and objects
from settings import equi_control,equi2_control,prod_control


# ------------------------------------------
# Main Routine
# Range settings 
# ------------------------------------------
exe = '/panfs/panasas01/chem/ab17369/codes/clean/dl-poly-4.10/build-intel16/bin/DLPOLY.Z' 
np = get_slurm_nproc('run.csh')
source_dir='3400k/equil2'
Tmin = 3400
Tmax = 2000
sign = get_sign(Tmin,Tmax) 
dT = sign*100    
NT = 1
alt_start=True

# -------------------------- ----------------------   
# DLPOLY input settings: CONFIG, FIELD and CONTROL 
# -------------------------------------------------
#equi_config  = ConfigFname(old='CONFIG', new='CONFIG')  # - If starting from input files that have not been run
equi_config  = ConfigFname(old='REVCON', new='CONFIG')
equi2_config = ConfigFname(old='REVCON', new='CONFIG')
prod_config  = ConfigFname(old='REVCON', new='CONFIG')
field='FIELD'
equi_control.temperature=Tmin
equi2_control.temperature=Tmin
prod_control.temperature=Tmin


# -------------------------- ----------------------   
# DLPOLY Calculations 
# -------------------------------------------------
#If calculation needs to start at equil2 or production, do so here  
if alt_start:
    if source_dir.split('/')[-1] == 'equil':
        run_dir=str(Tmin)+'k/equil2'
        source_dir = equilibration_calculation(Tmin,source_dir,run_dir,equi2_control,equi2_config,field,exe,np)
    if source_dir.split('/')[-1] == 'equil2':
        run_dir=str(Tmin)+'k/prod'
        source_dir = production_calculation(Tmin,source_dir,run_dir,prod_control,prod_config,field,exe,np)
        V_avg = extract_average_volume(run_dir+'/OUTPUT')
        V_rms =  extract_volume_rms(run_dir+'/OUTPUT')
        output_VvsT_to_file('VvsT.dat',V_avg,T, V_rms)
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
    V_rms =  extract_volume_rms(run_dir+'/OUTPUT')
    output_VvsT_to_file('VvsT.dat',V_avg,T, V_rms)

 
    

