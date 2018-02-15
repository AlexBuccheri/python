#!/usr/bin/env python3

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import math
import numpy as np


#Concatenates several string outputs to form the "dftb_in.hsd" file
def hsd_file_string(header,atomic_structure,GEO,HAM,PARS,KDETAILS,ANALYSE,CG_RELAX='None'):

    #Initalise strings 
    geometry_string, cg_string, hamiltonian_string, analysis_string, parser_string= '','','','',''
    
    #Geometry 
    if GEO.boundary_conditions.lower() == 'c':
        geometry_string = hsd_geometry_string_cluster(GEO) 
    if GEO.boundary_conditions.lower() == 's' or GEO.boundary_conditions.lower() == 'f':
        geometry_string = hsd_geometry_string_periodic(GEO)

    #Geometry relaxation
    if CG_RELAX !='None':
        cg_string = hsd_cg_string(CG_RELAX)
        
    #Hamiltonian (includes k-vectors)
    hamiltonian_string = hsd_hamiltonian_string(atomic_structure,GEO,HAM,KDETAILS)

    #Post-processing
    analysis_string=hsd_analysis_string(ANALYSE,GEO)
    
    #Parser version
    parser_string = hsd_parser_string(PARS)

    file_string=  header + "\n" + geometry_string + "\n" + cg_string + "\n" + hamiltonian_string + "\n" + analysis_string \
                         + "\n" + parser_string
    return file_string




#For periodic boundary conditions  
def hsd_geometry_string_periodic(GEO):

    geometry_string = GEO.header

    function_tag='hsd_geometry_string: '
    print(function_tag,'Multiplying basis & lattice vectors by lattice constant')
    
    #Multiply basis vectors by lattice constant and convert to strings
    al=GEO.al
    NBasis=len(GEO.basis_vectors)
    b_str=[]
    
    for i in range(0,NBasis):
        b=(al*GEO.basis_vectors[i,:])
        #Using this rather than 'join' to specify the formatting  
        b_str.append( np.array2string( b, separator=' ',formatter={'float_kind':lambda b: "%.8E" % b} ) )        

    l_str=[]
    for i in range(0,3):
        lvec=(al*GEO.lattice_vectors[i,:])
        l_str.append( np.array2string(lvec, separator=' ',formatter={'float_kind':lambda lvec: "%.8E" % lvec}) )
     
    #Generate gen format block     
    geometry_string=geometry_string+ """\n Geometry = GenFormat { 
         """+str(NBasis) + "  "+GEO.boundary_conditions+"""
         """ +' '.join(GEO.elements)

    for i in range(0,NBasis):
        if len(GEO.elements) ==1:       
            atom_index='1'
        else:
            atom_index=str(i+1)
        geometry_string=geometry_string + '\n         ' + str(i+1)+'  '+atom_index+'  '+b_str[i][1:-1]


    #Explicitly given origin as (0,0,0)
    geometry_string=geometry_string+ '\n         0.00000000E+00 0.00000000E+00 0.00000000E+00'

    for i in range(0,3):
        geometry_string=geometry_string + '\n         ' + l_str[i][1:-1]

    geometry_string=geometry_string + '\n   }'

        
    return geometry_string



def hsd_geometry_string_cluster(GEO):
    print('Write geometry string generation for finite boundary conditions')



def hsd_cg_string(CG_RELAX):
    cg_string =  \
    """ Driver = ConjugateGradient { 
         MovedAtoms = """+CG_RELAX.moved_atoms+""" 
         MaxForceComponent = """+str(CG_RELAX.max_force_component)+"""  
         MaxSteps = """+str(CG_RELAX.max_steps)+"""  
         OutputPrefix = """+CG_RELAX.output_prefix+"""
       }   """
    
    return cg_string


def hsd_hamiltonian_string(atomic_structure,GEO,HAM,KDETAILS):

    #SCC details 
    if KDETAILS.klines == None:
        SCC_string= """SCC = """ + HAM.scc +  """
    SCCTolerance = """ + str(HAM.scc_tolerance)+' \n'
    
    #Use charges from prior run if doing band structure (and therefore technically not SCC)
    if KDETAILS.klines != None:
        SCC_string= """SCC = Yes
    ReadInitialCharges = Yes
    MaxSCCIterations = 1 \n"""

    
    hamiltonian_string = """Hamiltonian = DFTB { 
    """+SCC_string+ \
    """    SlaterKosterFiles = Type2FileNames { 
        Prefix = """ + '"'+HAM.slaterkosterfiles.prefix+'"' + """
        Separator = """+'"'+HAM.slaterkosterfiles.separator+'"'+ """
        Suffix = """ +'"'+ HAM.slaterkosterfiles.suffix+'"' + """
    } 
    MaxAngularMomentum { \n"""
    
    for element in GEO.elements:
        hamiltonian_string = hamiltonian_string + '       '+element+' = "'+HAM.max_ang_momentum[element]+'" \n'   

    hamiltonian_string =  hamiltonian_string + """    }
    Filling = Fermi {
       Temperature ["""+HAM.tmp_unit+"""] = """ + str(HAM.tmp_f) + """
    } \n"""

    #Monkhorst-Pack k-grid 
    if KDETAILS.supercellfolding != None:
        hamiltonian_string =  hamiltonian_string + """    KPointsAndWeights = SupercellFolding {
        """+ str(KDETAILS.supercellfolding.kgrid[0,:])[1:-1]+ """
        """+ str(KDETAILS.supercellfolding.kgrid[1,:])[1:-1]+ """
        """+ str(KDETAILS.supercellfolding.kgrid[2,:])[1:-1]+ """
        """+ str(KDETAILS.supercellfolding.sweight)[1:-1]+"""
    } \n"""
        
    #K-vectors for band structure 
    if KDETAILS.klines != None:        
        #Unit conversion: Convert k-points from kin_unit to k_out_unit
        KDETAILS.klines.convert_kvector_units(atomic_structure,'reduced',KDETAILS.klines.kout_unit)
        hamiltonian_string += '#k-points output in units: '+KDETAILS.klines.kout_unit+' \n'
        kid=kunit_identifer(KDETAILS.klines.kout_unit)

        hamiltonian_string =  hamiltonian_string + """    KPointsAndWeights """+kid+""" = KLines {  \n """
        NKvectors=len(KDETAILS.klines.HS_points)
        for i in range(0,NKvectors):
           hamiltonian_string =  hamiltonian_string + "      "+ str(KDETAILS.klines.Nkpts_per_line[i]) \
           + '   '+str(KDETAILS.klines.HS_points[i,:])[1:-1]+'  \n'
        hamiltonian_string =  hamiltonian_string +'    } \n'
        
    #Explicit k-point grid specified
    if KDETAILS.explicitkpoints != None:
        print('Need to write parser for explicit k-grid settings. Skipping this')

    #Close Hamiltonian options
    hamiltonian_string =  hamiltonian_string +'} \n'
    
    #print (hamiltonian_string)
    return hamiltonian_string 



def kunit_identifer(kout_unit):
    if kout_unit=='absolute':
        kunit_identifer_string='[absolute]'
    if kout_unit=='reduced':
        print('Warning: Reduced units for outputted k-points is not valid for DFTB+')
        kunit_identifer_string=''
    if kout_unit=='fractional':
        kunit_identifer_string='[relative]'
    return kunit_identifer_string


def hsd_analysis_string(ANALYSE,GEO):

    dos_string=''
  
    if ANALYSE.output_DOS == True:
      dos_string= """ ProjectStates { \n"""
      
      for element in GEO.elements:
             dos_string = dos_string + """  Region {
     Atoms = """+element+"""
     ShellResolved = """+ANALYSE.shell_resolved+"""
     Label = "dos_"""+element+'" \n    } \n'
      dos_string=dos_string+'   }'
             
    analysis_string = """ Analysis {  \n""" + \
    dos_string+ """
  }  \n"""

    return analysis_string




def hsd_parser_string(PARS):
    parser_string=""" ParserOptions {
    ParserVersion = """ +str(PARS.version)+"""
  }"""
    return parser_string

