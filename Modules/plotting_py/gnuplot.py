#!/usr/bin/env python3

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import numpy as np
import os

#sys.path.insert(0, '/panfs/panasas01/chem/ab17369/2ndbs_tests')
#from Modules.classes import 


def generate_gnuplot_script(kband,fname='plot.p'):
    plot_string='set key off  \n'

    #K-point labels on x-axis
    #syntax=set xtics ("L" 1, "T" 40)
    plot_string += 'set xtics ('
    i=0
    ikpt=0
    for point in kband.HS_points_list:
        ikpt=ikpt+kband.Nkpts_per_line[i]
        #if kband.Nkpts_per_line[i]-kband.Nkpts_per_line[i-1]==1:
        plot_string += '"'+point+'" '+str(ikpt)+', '    
        i=i+1

    #Remove final trailing comma
    plot_string=plot_string[:-2]
    plot_string+=')  \n'

    plot_string+="set ylab'Energy (eV?)'   \n"

    #Plot bands
    plot_string += "plot 'band_tot.dat' u 1:2 w l  \n"

    for i in range(1,9):
        plot_string += "replot 'band_tot.dat' u 1:"+str(i+2)+" w l  \n"

    #print(plot_string)
    fid= open(fname, "w")
    fid.write(plot_string)
    fid.close()
