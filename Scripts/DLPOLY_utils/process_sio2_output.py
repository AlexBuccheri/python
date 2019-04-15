#!/usr/bin/env python3

# --------------------------------------------------------
# Extract cell vectors from DLPOLY output file, and plot
# cell parameters and volume as a function of temperature
# --------------------------------------------------------
__author__ = "Alexander Buccheri <ab17369@bristol.ac.uk>"
__copyright__ = "Copyright© 2019 A Buccheri"
__license__ = "GPL-3.0-only"
__version__ = "1.0"
# --------------------------------------------------------

#Libraries
import sys
import numpy as np
import os
import subprocess
import shutil
import matplotlib.pyplot as plt


# -------------------------------------
# Functions 
# -------------------------------------

def strlist_to_floatlist(strlist):
    return [float(i) for i in strlist.split()]

def list_to_floatlist(strlist):
     return [float(i) for i in strlist]

# Units of ang. Works for CONFIG and REVCON 
def read_cell_vectors_from_config(fname):
    f = open(fname)
    #Skip junk 
    for i in range(0,2):
        line = f.readline()
    #Get cell vectors in string 
    line = f.readline()
    for i in range(0,2):  
        line += f.readline()
    f.close()
    return line

def read_cell_vectors_from_output(fname):
    if fname[-6:].lower() != 'output':
        print('Filename must end in OUTPUT')
        sys.exit('Script has stopped')
    #-a required to make grep proceed if it detects binary symbols
    output_bytes = subprocess.check_output("grep -a -A 3 'Average cell vectors' "+fname ,\
                                          shell=True)
    #Convert from utf-8 byte literal to string literal
    output = output_bytes.decode("utf-8")
    #Convert to str-list and remove preceding text (fixed formatting)
    output = output.split()[6:]
    return output

#Taken from Alin's statis script
#Example 
#nrdf,nprdf,rdf_data,rdf_labels = readRDF(ROOT+'300k/prod/RDFDAT')
#nrdf,npoints,npoints => x=pair type, y=separation,  z=RDF value 
# rdf_data[0,:,:] = Si-Si
# rdf_data[1,:,:] = Si-O
# rdf_data[2,:,:] = O-O
#
def readRDF(filename="RDFDAT"):
    try:
      title, header, rdfall = open(filename).read().split('\n',2)
    except IOError:
        return 0,0,0,[]
    nrdf,npoints=map(int, header.split())
    b=2*(npoints+1)
    d=np.zeros((nrdf,npoints,2),dtype=float)
    labels=[]
    s=rdfall.split()
    for i in range(nrdf):
      x=s[b*i:b*(i+1)]
      y=np.array(x[2:],dtype=float)
      y.shape= npoints,2
      d[i,:,:]=y
      labels.append(x[0]+" ... "+x[1])
    return nrdf,npoints,d,labels


def get_cell_vectors(fname):
    # Cell vector string 
    line = read_cell_vectors_from_output(fname)
    # 3x(Cell vector + RMS vector) -> split up
    cell_lines=[]
    rms_lines=[]
    index = [0,1,2,6,7,8,12,13,14]
    for i in index:
        cell_lines.append(line[i])
        rms_lines.append(line[i+3])
    cell_vectors = np.asarray( list_to_floatlist(cell_lines) )
    rms_vectors  = np.asarray( list_to_floatlist(rms_lines)  )
    return np.reshape(cell_vectors,(3,3)), np.reshape(rms_vectors,(3,3))

# Get lattice constants 
def get_lattice_constants(cell_vectors):
    a = np.linalg.norm(cell_vectors[0,:])
    b = np.linalg.norm(cell_vectors[1,:])
    c = np.linalg.norm(cell_vectors[2,:])
    return a,b,c

def degrees_to_radians(theta):
    return (theta/360.)*2.*np.pi

# a,b,c are vectors 
def triple_product(a,b,c):
    return np.dot(np.cross(a,b),c)


def get_diffusion_coefficients(fname):
    if fname[-6:].lower() != 'output':
        print('Filename must end in OUTPUT')
        sys.exit('Script has stopped')
    output_bytes = subprocess.check_output("grep -a -A 3 'Approximate 3D Diffusion Coefficients' "+fname ,\
                                          shell=True)
    output = output_bytes.decode("utf-8")
    output = output.splitlines()
    data = []
    label = []
    for i in range(2,4):
        label.append( output[i].split()[0] )
        data.append( output[i].split()[1:3] )
    return label, np.asarray(data)   
                                 



# ------------------------
# Main Script
# ------------------------
# Number of unit cells per dimension of the supercell 
ncells = 5.

# ------------------------
# Cell constants 
# ------------------------
Ti = 300
Tf = 1300
dT = 100
Nt = int((Tf-Ti)/dT)+1

sabc = np.zeros( shape=(3,Nt) )
sabc_rms = np.zeros( shape=(3,Nt) )
svol = np.zeros( shape=(Nt) )
ROOT = '../silicate/'


#One needs to extract cell parameters from OUTPUT, not REVCON
for i in range(0,Nt):
    T = Ti + i*dT
    #Production output
    fname = ROOT+str(T)+'k/prod/OUTPUT'
    scell,rms = get_cell_vectors(fname)
    sabc[:,i] = get_lattice_constants(scell)
    sabc_rms[:,i] = get_lattice_constants(rms)
    svol[i] = triple_product(scell[0,:],scell[1,:],scell[2,:])
    #print(T, sabc[:,i]/ncells, svol[i]/ncells**3.)

    
#Production (unit) cell constants vs temperature
plot_me = False
if plot_me:
    T = np.arange(Ti, Tf+dT, dT)
    plt.xlabel("Temperature (k)")
    plt.ylabel("Cell constants (angstrom)")
    plt.errorbar(T, sabc[0,:]/ncells, xerr=None, yerr=sabc_rms[0,:]/ncells, fmt='ro--', label='a')
    plt.errorbar(T, sabc[1,:]/ncells, xerr=None, yerr=sabc_rms[1,:]/ncells, fmt='bs--', label='b')
    plt.errorbar(T, sabc[2,:]/ncells, xerr=None, yerr=sabc_rms[2,:]/ncells, fmt='g^--', label='c')
    plt.legend(loc='upper left')
    plt.show()

#Production volume (unit cell)
plot_me = False
percentage = True
if plot_me:
    T = np.arange(Ti, Tf+dT, dT)
    if percentage:
        svol = (svol[:]-svol[0])*100./svol[0]#for % increase
        plt.ylabel("Unit Cell Volume Increase (%)")
        plt.xlabel("Temperature (k)")
        plt.plot(T,svol , 'ro--')
        plt.show()
    else:
        plt.ylabel("Unit Cell Volume (ang^3)")
        plt.xlabel("Temperature (k)")
        plt.plot(T,svol/ncells**3. , 'ro--')
        plt.show()


# ----------------------------------------------
# Experimental Cell constants
# See exp_abc.dat file. Refs [2] and [6] removed 
# ----------------------------------------------
exp_temp = np.array( [291,293,293,293,293,293,293] )
exp_abc = np.array( [ [4.9130,  4.9130,   5.4047 ],\
                      [4.914,   4.914,    5.406  ],\
                      [4.921,   4.921,    5.4163 ],\
                      [4.9158,  4.9158,   5.4091 ],\
                      [4.91444, 4.91444,  5.40646],\
                      [4.921,   4.921,    5.416  ],\
                      [4.9148,  4.9148,   5.4062] ] )
                        
#Compare room temp/ambient pressure exp lattice constants to MD structure
plot_me = False
if plot_me:
    # a,b
    plt.plot(exp_temp, exp_abc[:,0], 'ro', label='exp')
    plt.plot(300, sabc[0,0]/ncells, 'bs', label='production')
    plt.legend(loc='upper right')
    plt.xlabel("Temperature (k)")
    plt.ylabel("Cell constants a & b (ang)")
    plt.show()
    # c
    plt.plot(exp_temp, exp_abc[:,2], 'ro', label='exp')
    plt.plot(300, sabc[2,0]/ncells, 'bs', label='production')
    plt.xlabel("Temperature (k)")
    plt.ylabel("Cell constant c (ang)")
    plt.legend(loc='upper right')
    plt.show()


# ---------------------------------------
# RDFs - Production output and experimental  
# ---------------------------------------
# Average RDFs from the database (which ref did I take this from?)
exp_rdf={ 'si-si':3.0523,'si-o':1.6012, 'o-o':2.5977}
prod_rdf={ 'si-si':3.077,'si-o':1.625, 'o-o':2.65} 


#Max RDF value w.r.t. temperature
bond_si_si = np.zeros(shape=(Nt))
bond_si_o  = np.zeros(shape=(Nt))
bond_o_o   = np.zeros(shape=(Nt))

for i in range(0,Nt):
    # Read RDF data 
    T = Ti + i*dT
    fname = ROOT+str(T)+'k/prod/RDFDAT'
    nrdf,nprdf,rdf_data,rdf_labels = readRDF(fname)
    #Could turn into function and make more readable 
    max_peak_index = np.argmax(rdf_data[0,:,1]) 
    bond_si_si[i] = rdf_data[0,max_peak_index,0] 
    max_peak_index = np.argmax(rdf_data[1,:,1]) 
    bond_si_o[i] = rdf_data[1,max_peak_index,0] 
    max_peak_index = np.argmax(rdf_data[2,:,1]) 
    bond_o_o[i] = rdf_data[2,max_peak_index,0] 
    
#Plot the peak RDF value against temperature
plot_me = False 
if plot_me:
    T = np.arange(Ti, Tf+dT, dT)
    plt.xlabel("Temperature (k)")
    plt.ylabel("Max RDF value (arb)")
    plt.plot(T,bond_si_si, 'ro--', label='Si-Si')
    plt.plot(T,bond_si_o, 'bs--', label='Si-O')
    plt.plot(T,bond_o_o, 'g^--', label='O-O')
    plt.legend(loc='upper left')
    plt.show()

#Plot DLPOLY RDF for 300K, and add experimental data peaks 


# ---------------------------------------
# Diffusion Coefficients 
# ---------------------------------------

nspecies = 2
d_label = np.zeros(shape=(Nt,nspecies), dtype=str)
# Temperature, species (Si, O), (D, error)
dcoeff  = np.zeros(shape=(Nt,nspecies,2))  
    
for i in range(0,Nt):
    T = Ti + i*dT
    fname = ROOT+str(T)+'k/prod/OUTPUT'
    d_label[i,:], dcoeff[i,:,:] = get_diffusion_coefficients(fname)
    #print( d_label[i,:], dcoeff[i,:,0:] )


plot_me = True  
if plot_me:
    T = np.arange(Ti, 3500+dT, dT)
    N = int((3500-Ti)/dT)+1
    plt.xlabel("Temperature (k)")
    (r'$\alpha > \beta$')
    plt.ylabel(r'Diffusion coefficient ($10^{-9} m^2 s^{-1}$)')
    plt.plot(T, dcoeff[0:N,0,0], 'ro--', label='Si')
    plt.plot(T, dcoeff[0:N,1,0], 'bo--', label='O')
    plt.legend(loc='upper left')
    plt.show()

    T = np.arange(3500, Tf+dT, dT)
    plt.xlabel("Temperature (k)")
    plt.ylabel(r'Diffusion coefficient ($10^{-9} m^2 s^{-1}$)')
    plt.plot(T, dcoeff[N-1:,0,0], 'ro--', label='Si')
    plt.plot(T, dcoeff[N-1:,1,0], 'bo--', label='O')
    plt.legend(loc='upper left')
    plt.show()
