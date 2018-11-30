#!/usr/bin/env python3

#Import libraries
import math
import numpy as np
import sys



class Constants:
     boundary_key=np.array(['finite','cubic','orthorhombic','parallel-piped','Not valid integer',\
                            'Not valid integer','xy but not z'])
     atomic_record=np.array(['coorindates', 'coordinates & velocities', 'coordinates, velocities & forces'])

                            
#Data found in DLPOLY CONFIG file 
class Config:
    def __init__(self,header=None,Nlines_per_record=1,boundary_index=1,Natoms=0, \
                      lattice_vector=None,atom_name=None,coord=None,velo=None,force=None,atom_index=None):
                                                  #DL POLY Variable names
       self.header=header        
       self.Nlines_per_record=Nlines_per_record
       self.boundary_index=boundary_index
       self.Natoms=Natoms
       self.lattice_vector=lattice_vector  #Stored row-wise
       self.atom_name=atom_name
       self.coord=coord     #xxx, yyy and zzz
       self.velo=velo       #vxx, vyy and vzz 
       self.force=force     #fxx, fyy and fzz
       self.atom_index=atom_index

       #Does it make sense to initialise shapes?
       if lattice_vector is None:  self.lattice_vector=np.zeros(shape=(3,3))   
       if atom_name is None:       self.atom_name=[]         
       if coord is None:           self.coord=np.zeros(shape=(3,1))
       if atom_index is None:      self.atom_index=np.zeros(shape=(1))  
 
       # print('Boundary conditions:',Constants.boundary_key[boundary_index])
       # print('Each atomic record is comprised of ',Constants.atomic_record[Nlines_per_record-1])

class Ensemble:
    def __init__(self,name='nve',etype=None,thermostat_relaxation=None,barostat_relaxation=None):
        self.name = name
        self.etype = etype
        #Depending on context, times (ps) or speeds (ps^-1)
        self.f1 = thermostat_relaxation
        self.f2 = barostat_relaxation

class Trajectory:
    def __init__(self,tstart=0,tinterval=1,data_level=0):
        self.tstart = tstart
        self.tinterval = tinterval
        self.data_level = data_level

#Data found in DLPOLY CONTROL file - Extend as required 
class Control:
    def __init__(self,header=' ', name='CONTROL', temperature=None, pressure=None, ensemble=None, \
                     steps=None, equilibration=None, scale=None, regauss=None, printout=None,\
                     stack=None, stats=None, vdw_direct=None, rdf=None, print_rdf=None, trajectory=None,\
                     optimise_force=None, ewald_precision=None, timestep=None, rpad=None, cutoff=None,\
                     cap=None, job_time=None, close_time=None, dump=None):
        self.header = header
        self.name = name 
        self.temperature = temperature
        self.pressure = pressure
        self.ensemble = ensemble
        self.steps = steps 
        self.equilibration = equilibration
        self.scale = scale
        self.regauss = regauss
        self.printout = printout
        self.stack = stack
        self.stats = stats
        self.vdw_direct = vdw_direct
        self.rdf = rdf
        self.print_rdf = print_rdf
        self.trajectory = trajectory
        self.optimise_force = optimise_force
        self.ewald_precision = ewald_precision
        self.dump = dump
        self.timestep = timestep
        self.rpad = rpad
        self.cutoff = cutoff
        self.cap = cap
        self.job_time = job_time
        self.close_time = close_time


#One assumes that there is a more sensible way of doing this, like
#being able to iterate through all atributes of the data object
#Something to discuss with Peter 
def Control_data_to_string(data):
        config_string=''
        padding='                       '
        if(data.header      != None): config_string+=data.header+'\n'
        if(data.temperature != None): config_string+='temperature'+padding+str(data.temperature)+'\n'
        if(data.pressure    != None): config_string+='pressure'+padding+str(data.pressure)+'\n'
        if(data.ensemble    != None):
            config_string+='ensemble '+str(data.ensemble.name)
            if(data.ensemble.etype !=None): config_string+=str(data.ensemble.etype)+'  '    
            if(data.ensemble.f1 !=None): config_string+=str(data.ensemble.f1)+' '    #This will break formatting
            if(data.ensemble.f2 !=None): config_string+=str(data.ensemble.f2)+'\n'   #if not f1 and f2
            
        if(data.steps  != None): config_string+='steps'+padding+str(data.steps)+'\n'
        if(data.equilibration != None): config_string+='equilibration'+padding+str(data.equilibration)+'\n'
        if(data.scale   != None): config_string+='scale'+padding+str(data.scale)+'\n'
        if(data.regauss != None): config_string+='regauss'+padding+str(data.regauss)+'\n'
        if(data.printout != None): config_string+='print'+padding+str(data.printout)+'\n'
        if(data.stack != None): config_string+='stack'+padding+str(data.stack)+'\n'
        if(data.stats != None): config_string+='stats'+padding+str(data.stats)+'\n'
        if(data.vdw_direct == True): config_string+='vdw direct    \n'
        if(data.rdf != None): config_string+='rdf'+padding+str(data.rdf)+'\n'
        if(data.print_rdf == True): config_string+='print rdf    \n' 
        if(data.trajectory    != None):
            config_string+='trajectory '+padding+str(data.trajectory.tstart)+' '+str(data.trajectory.tinterval)+' '+str(data.trajectory.data_level)+'\n'
            
        if(data.optimise_force != None): config_string+='optimise_force'+padding+str(data.optimise_force)+'\n'
        if(data.dump !=None):  config_string+='dump'+padding+str(data.dump)+'\n'
        if(data.ewald_precision != None): config_string+='ewald_precision'+padding+str(data.ewald_precision)+'\n'
        if(data.timestep != None): config_string+='timestep'+padding+str(data.timestep)+'\n'
        if(data.rpad != None): config_string+='rpad'+padding+str(data.rpad)+'\n'
        if(data.cutoff != None): config_string+='cutoff'+padding+str(data.cutoff)+'\n'
        if(data.cap != None): config_string+='cap'+padding+str(data.cap)+'\n'
        if(data.job_time  != None): config_string+='job_time'+padding+str(data.job_time)+'\n'
        if(data.close_time  != None): config_string+='close_time'+padding+str(data.close_time)+'\n'
        config_string+='finish'
        return config_string

        

#Take a CONFIG file parsed as a single string and assign data to DL POLY class 
def Config_string_to_data(string, data):

    data.header = string[0]
    Nlines_per_record, boundary_index, Natoms = string[1].split()
    
    data.Nlines_per_record=int(Nlines_per_record)+1
    data.boundary_index=int(boundary_index)
    data.Natoms=int(Natoms)
    #Allocate correct array size for coordinates 
    if data.coord.shape[1] < data.Natoms:
        data.coord.resize((3,data.Natoms))

    #Lattice vectors     
    for i in range(2,5):
        j=i-2
        data.lattice_vector[j,0:3] = string[i].split()  #No arg => whitespace delimiter

    #Atomic records for Natoms, each comprised of 'Nlines_per_record+1' lines
    offset=5
    k=offset
    
    if data.Nlines_per_record == 1:
        for ia in range(0,data.Natoms):
            atom_name,dummy_index = string[k].split()
            data.atom_name.append(atom_name)
            k=k+1
            data.coord[:,ia] = string[k].split()
            k=k+data.Nlines_per_record

    if data.Nlines_per_record == 2:
        print('Needs testing')
        if data.velo.shape(1) < data.Natoms: data.velo.resize((3,data.Natoms))
            
        for ia in range(0,data.Natoms):
            atom_name,dummy_index = string[k].split()
            data.atom_name.append(atom_name)
            k=k+1
            data.coord[:,ia] = string[k].split()
            data.velo[:,ia]  = string[k+1].split()
            k=k+data.Nlines_per_record

    if data.Nlines_per_record == 3:
        print('Needs testing')
        if data.velo.shape(1) < data.Natoms: data.velo.resize((3,data.Natoms))
        if data.force.shape(1) < data.Natoms: data.force.resize((3,data.Natoms))
            
        for ia in range(0,data.Natoms):
            atom_name,dummy_index = string[k].split()
            data.atom_name.append(atom_name)
            k=k+1
            data.coord[:,ia] = string[k].split()
            data.velo[:,ia]  = string[k+1].split()
            data.force[:,ia] = string[k+2].split()
            k=k+data.Nlines_per_record


# Convert DLPOLY CONFIG data into a string            
def config_data_to_string(config):
    config_string=''
    padding='          '

    if isinstance(config.header, str): config_string=config.header+'\n'

    config_string += padding+str(config.Nlines_per_record)+'    '+\
                     str(config.boundary_index)+'      '+\
                     str(config.Natoms)+'\n'

    for ia in range(0,3):
        vec = config.lattice_vector[ia, :]
        lat_str = np.array2string(vec, separator='    ', formatter={'float_kind': lambda vec: "%.8E" % vec})
        config_string += padding + lat_str[1:-1]+'\n'


    if config.Natoms==0: config.Natoms=len(config.coord)

    if config.Nlines_per_record != 1:
        print('Have not written CONFIG to string routine to output when'
              'Nlines_per_record exceeds 1')
        sys.exit('Script has stopped')

    if config.Nlines_per_record==1:
        for ia in range(0,config.Natoms):
            config_string+=config.atom_name[ia]+padding[:-len(config.atom_name[ia])-1]+str(ia+1)+'\n'
            pos = config.coord[ia,:]
            # Using this rather than 'join' to specify the formatting TURN INTO FUNCTION
            pos_str=np.array2string(pos, separator='    ', formatter={'float_kind': lambda pos: "%.8E" % pos})
            config_string += padding+pos_str[1:-1]+'\n'



    return config_string
