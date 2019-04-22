#!/usr/bin/env python3

import copy
#My libraries
from Modules.dlpoly_py import dl_classes as dl

# ---------------------------------------------------------------
# Control settings for equilibration and production
# classical silicate calculations
# ---------------------------------------------------------------

#Control data for equilibriation calculation 
equi_ensemble   = dl.Ensemble(name='npt', etype='hoover', thermostat_relaxation=1.0, barostat_relaxation=0.5)
equi_trajectory = dl.Trajectory(tstart=10000,tinterval=1000,data_level=0)

equi_control = dl.Control( \
    header='Quartz',\
    temperature=None,\
    pressure=0.001,\
    ensemble=equi_ensemble,\
    steps=500000,\
    equilibration=500000,\
    scale=7,\
    regauss=13,\
    printout=1000,\
    stack=1000,\
    stats=1000,\
    vdw_direct=True,\
    rdf=1000,\
    print_rdf=True,\
    trajectory=equi_trajectory,\
    ewald_precision=1.0e-6,\
    timestep=0.001,\
    rpad=0.4,\
    cutoff=6.0,\
    cap=1.0e3,\
    job_time=9.0e5, close_time=2.0e1)

#Initialise control data for 2nd equilibriation calculation, where velocities aren't rescaled
equi2_control = copy.deepcopy(equi_control)
equi2_control.scale     =None
equi2_control.regauss   =None
equi2_control.rdf       =None
equi2_control.trajectory=None
equi2_control.dump      =1000000


#Initialise control data for production calculation 
prod_ensemble   = dl.Ensemble(name='npt', etype='hoover', thermostat_relaxation=1.0, barostat_relaxation=0.5)
prod_trajectory = dl.Trajectory(tstart=0,tinterval=1000,data_level=0)
prod_control    = dl.Control(
    header='Quartz',\
    temperature=None,\
    pressure=0.001,\
    ensemble=prod_ensemble,\
    steps=1000000,\
    equilibration=0,\
    printout=1000,\
    stack=1000,\
    stats=1000,\
    vdw_direct=True,\
    rdf=1000,\
    print_rdf=True,\
    trajectory=prod_trajectory,\
    dump=1000000,\
    ewald_precision=1.0e-6,\
    timestep=0.001,\
    rpad=0.4,\
    cutoff=6.0,\
    cap=1.0e3,\
    job_time=9.0e5, close_time=2.0e1)
