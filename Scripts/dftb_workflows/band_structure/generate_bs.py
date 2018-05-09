#!/usr/bin/env python3

#########################################################################################
# Alexander Buccheri. University of Bristol. Feb 2018
#
# Sets up directories and dftb+ input scripts for band structure calculations
# First calculation for converged charges (analogous to density in DFT)
# Second calculation for band structure calculation using charges from 1st calculation
#
#########################################################################################

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import re
import sys
import math
import numpy as np
import os
import matplotlib.pyplot as plt
import platform


#Paths to my modules (alternatively put directory in $PATH)
#Blue crystal varies according to specific log-in node
blue3='Linux-2.6.32-642.6.2.el6.x86_64-x86_64-with-redhat-6.4-Carbon'
AlexMac='Darwin-17.4.0-x86_64-i386-64bit'

if platform.platform()[0:12] == blue3[0:12]:
    sys.path.insert(0, '/panfs/panasas01/chem/ab17369/python_modules')
if  platform.platform()[0:5] == AlexMac[0:5]:
    sys.path.insert(0, '/Users/alexanderbuccheri/Python')
    
#My Modules
from Modules.dftb_py.hsd     import hsd_file_string
from Modules.dftb_py.classes import Geometry,Hamiltonian,SlaterKosterFiles,SuperCellFolding,KLines,Ksampling,\
                                    ParserOptions,Analysis,Lattice,ConjugateGradient,Cubic,HS_labels_2_points
from Modules.plotting_py     import gnuplot





#----------------------------------
#Material options
#----------------------------------
structure='cubic'
basis='prim'
elements = np.array(['Si']) #Order should be consistent with anion at (0,0,0) and cation at (0.25,0.25,0.25)
material='Si'
max_ang_momentum= {"Si":'p'}
boundary_conditions='S'
parameter_dir='/panfs/panasas01/chem/ab17369/codes/dftbplus/parameters/siband/siband-1-1/'
kgrid=np.array([[ 24.,0., 0.],[ 0.,24., 0.],[ 0.,0., 24.]])
generate_directories=True



#----------------------------------------------------------------------------------------------
# Initialise geometry (both calculations below with same boundary condition, hence define GEO once)
#----------------------------------------------------------------------------------------------
atomic_structure=Lattice(structure=structure,basis=basis,material=material)

GEO=Geometry(material,elements,boundary_conditions,\
             al=atomic_structure.al,basis_vectors=atomic_structure.basis_vector,\
             lattice_vectors=atomic_structure.lattice_vector,header=atomic_structure.header)


#-------------------------------------
#Do 2-step band structure calculation
#-------------------------------------
if generate_directories == True:
   #Set up directories 
   input_name='dftb_in.hsd'
   root_dir=material+'_bs'
   
   itr=0
   tmp_dir=root_dir
   while os.path.isdir(tmp_dir) == True:
       itr=itr+1
       tmp_dir=str(itr)+'.'+root_dir
        
   root_dir=tmp_dir         
   sub_dir1=root_dir+'/1.converged_charges'
   sub_dir2=root_dir+'/2.bands'
        
   os.mkdir(root_dir)
   os.mkdir(sub_dir1)
   os.mkdir(sub_dir2)



#----------------------------------------------------------------
# 1. Converged charges
#    Initialise option objects for converged charges calculation
#    KDETAILS = superfluous class
#----------------------------------------------------------------
SKF=SlaterKosterFiles(prefix=parameter_dir)
HAM=Hamiltonian(scc='Yes',slaterkosterfiles=SKF,max_ang_momentum=max_ang_momentum,scc_tolerance=1.e-10)
SCELL=SuperCellFolding(structure=structure,kgrid=kgrid)  #If kgrid argument is not specified, object uses default MP grid for cubic
KDETAILS=Ksampling(supercellfolding=SCELL)          
ANALYSE=Analysis(output_DOS=False)
PARS=ParserOptions(version=5)
CG_RELAX=ConjugateGradient()                         #Use all defaults for atomic relaxation

#Generate options string and write to file
header='#Converged SCC DFTB+ calculation on '+GEO.material
file_string = hsd_file_string(header,atomic_structure,GEO,HAM,PARS,KDETAILS,ANALYSE)  # ,CG_RELAX)

#Write string to file or print to screen
if generate_directories == False:
    print(file_string)

if generate_directories == True:
    fid= open(sub_dir1+'/'+input_name, "w")
    fid.write(file_string)
    fid.close()


#----------------------------------
# 2. Band structure calculation
#----------------------------------
HAM=Hamiltonian(scc='Yes',slaterkosterfiles=SKF,max_ang_momentum=max_ang_momentum)

#Would ideally like Cubic to be assigned internally when 'structure' is set, so just choosing Lpoint, Xpoint, etc
HS_points_list=np.array(['L','T','X','U','K','T'])
points_per_line=np.array([1,40,50,20,1,50])

HS_points=HS_labels_2_points(HS_points_list)
KBAND=KLines(atomic_structure, HS_points_list, HS_points, points_per_line, kin_unit=Cubic
.k_units,kout_unit='fractional')
KDETAILS=Ksampling(klines=KBAND)
ANALYSE=Analysis(output_DOS=True,shell_resolved='Yes')

file_string = hsd_file_string(header,atomic_structure,GEO,HAM,PARS,KDETAILS,ANALYSE)

if generate_directories == False:
    print(file_string)

if generate_directories == True:
    fid= open(sub_dir2+'/'+input_name, "w")
    fid.write(file_string)
    fid.close()
    gnuplot.generate_gnuplot_script(KBAND,fname='plot.p')




