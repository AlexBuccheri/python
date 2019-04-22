#!/usr/bin/env python3

# --------------------------------------------------------
# Read in data from several silicate runs and average
# --------------------------------------------------------
__author__ = "Alexander Buccheri <ab17369@bristol.ac.uk>"
__copyright__ = "Copyright 2019 A Buccheri"
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


configs = ['1064','1127','1502','1644','1703','1875',\
           '333','527','728','73']

configs = ['config'+x for x in configs]
Nconfigs = float(len(configs))

#Initial data point
fname = configs[0]+'/VvsT.dat'
data = np.loadtxt(fname, skiprows=1)
npt = np.size(data,0)
tdata = np.zeros(shape=(npt))
tdata = data[:,0] 
avg_vdata = np.zeros(shape=(npt))
avg_vdata = data[:,1]

#All other data points 
configs.pop(0)
print(configs)
for config in configs:
    fname = config+'/VvsT.dat'
    print(config,fname)
    data = np.loadtxt(fname, skiprows=1)
    avg_vdata += data[:,1]
avg_vdata = avg_vdata/Nconfigs

print('#Temperature (k), avg volume (ang^3):')
for i in range(0,npt):
    print(tdata[i], avg_vdata[i])
