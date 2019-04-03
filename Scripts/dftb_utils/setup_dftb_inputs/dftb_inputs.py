#!/usr/bin/env python3

#Modules 
from ase.io import read,write
import sys
import numpy as np
import os

#My Modules
#from Modules.dftb_py.hsd    import hsd_file_string
from Modules.scheduler_py   import slurm
#from Modules.dftb_py.classes import Geometry,Hamiltonian,SlaterKosterFiles,SuperCellFolding,KLines,Ksampling,\
#                                    ParserOptions,Analysis,Lattice,ConjugateGradient,Cubic,HS_labels_2_points

#Number of atoms in supercell from basis
def atoms_in_supercell(nx,ny,nz,basis):
    return int(nx*ny*nz*basis) 


def gen_string(filename):
    string = \
    '''Geometry = GenFormat {
   <<< "'''+filename+'''" 
} \n'''
    return string 


#--------------------------------------------------------------------------
# Generate DFTB+ .hsd and .gen files, and slurm run scripts 
# from .cif input, ASE and my modules
#
# Set up several supercells for non-SCC, total energy and force calculation 
#--------------------------------------------------------------------------

#Supercell integer range 
cells = np.arange(5,10,1)
#Get data
a=read("sio2.cif")
#ASE ordering consistent with cif: '3 Si' then '6 O' = 9 atoms
basis = len( a.arrays['numbers'] )
#Format to output structure in 
output = 'dftb'

hsd_string = \
'''
Hamiltonian = DFTB {
  SCC = Yes
  MaxSCCIterations = 5

  MaxAngularMomentum = {
    Si = "p"
    O  = "p"
  }

  Filling = Fermi {
    Temperature [Kelvin] = 100
  }

  SlaterKosterFiles = Type2FileNames {
    Prefix = "/mnt/storage/home/ab17369/parameters/slater_koster/matsci/matsci-0-3/"
    Separator = "-"
    Suffix = ".skf"
  }

  KPointsAndWeights = SuperCellFolding {
    1 0 0
    0 1 0
    0 0 1
    0.0 0.0 0.0
  }
}

Options = {
   WriteResultsTag=Yes
   WriteDetailedOut=No
   TimingVerbosity=-1
}

Analysis = {
  CalculateForces = Yes
}

ParserOptions {
  ParserVersion = 6
} '''

#DFTB+ hsd settings - easier to just write string and modify structure input 
#max_ang_momentum= {"Si":'p',"O":'p'}
#kgrid=np.array([[1,0, 0],[ 0,1,0],[ 0,0, 1]])
#parameter_dir='/panfs/panasas01/chem/ab17369/codes/dftbplus/parameters/siband/siband-1-1/'
#sk=SlaterKosterFiles(prefix=parameter_dir)
#ham=Hamiltonian(scc='No',slaterkosterfiles=sk,max_ang_momentum=max_ang_momentum)
#ksettings=SuperCellFolding(structure='zb',kgrid=kgrid) 
#kgrid =Ksampling(supercellfolding=ksettings)
#analysis=Analysis()
#parser=ParserOptions(version=6)
#cg=ConjugateGradient()

nodes=1
exe='/mnt/storage/home/ab17369/codes/dftbplus-18.2/_build/prog/dftb+/dftb+'
wt=[1,0,0]
queue='veryshort'

for ncell in cells:
    natoms = atoms_in_supercell(ncell,ncell,ncell,basis)
    print(ncell,natoms)
    
    #DLPOLY supercell
    if output == 'dl':
        write( str(ncell)+"CONFIG",a.repeat(int(ncell)) )
        print('Modify in FIELD file: ATOMS 9, NUMMOLS ',ncell**3)
        
    #DFTB
    if output == 'dftb':
        root = str(natoms)
        os.mkdir(root)

        #geometry file
        gen_file = 'structure_'+str(natoms)+'.gen'
        write(root+'/'+gen_file, a.repeat(int(ncell)) )

        #hsd file
        output_str = gen_string(gen_file)
        output_str += hsd_string
        fid = open(root+'/dftb_in.hsd','w')
        fid.write(output_str)
        fid.close()

        #slurm run script 
        job_name = 'dftb_'+str(natoms)
        slurm_setup = slurm.SlurmJob(exe=exe,nodes=nodes,job_name=job_name,queue=queue,walltime=wt,output='terminal')
        slurm_str = slurm.generate_slurm_script( slurm_setup ) 
        fid = open(root+'/run.csh','w')
        fid.write(slurm_str)
        fid.close()



