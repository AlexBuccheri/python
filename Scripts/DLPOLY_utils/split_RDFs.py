#!/usr/local/bin/python3

################################################################
# Process DLPOLY Output 'RDFDAT' into separate files
# Alexander Buccheri. Sept 2018 
################################################################

#Import libraries
import math
import numpy as np
import sys
import os

def get_Nlines(fname):
    line=os.popen('wc -l '+fname).read()
    return int(line.lstrip().split(' ')[0])


#Main Routine
f='RDFDAT'
Nl=get_Nlines(f)
data = np.genfromtxt(f, dtype=str,  delimiter='\n')

#Number of diatomic pairs in RDF file 
Npairs=int( data[1].split()[0] )
#Number of bins per atomic pair
Nbin=int( data[1].split()[1] )
header='# '+data[0]

#Remove non-repeating header lines
data = np.delete(data, [0,1])

for ipair in range(0,Npairs):
    i1=(Nbin+1)*ipair
    
    element = data[i1].split()
    pair_label = str(element[0]+'_'+element[1])
    fname=f+'_'+pair_label+'.dat'
    fout=open(fname, mode='wt', encoding='utf-8')
    print(header,file=fout)
    print('# '+data[i1],file=fout)
    
    for ibin in range(i1+1,i1+Nbin+1):
        print(data[ibin],file=fout)
        
    fout.close()



    

# Check this command out
# import string
# ALPHA = string.ascii_letters
# if line.startswith(tuple(ALPHA)) == False:
#    print (line,file=f2)
    
        


    
