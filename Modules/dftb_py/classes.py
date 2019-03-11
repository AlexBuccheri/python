#!/usr/bin/env python3

#----------------------------------
#Libraries and modules
#----------------------------------
#Libraries 
import sys
import math
import numpy as np

#----------------------------------
#Classes
#----------------------------------

      
#---------------------------------------------------------------
# Container class for cubic 
# fcc real space/bcc reciprocal space lattice
# All vectors in crystalographic (reduced) units of the lattice
#---------------------------------------------------------------
class Cubic:
    header='#Basis and lattice vectors for primitive ZB cell, in units of Ang (DFTB+ default)'

    #Primitive basis and lattice vectors
    prim_basis_vectors = np.array(  [ [0.  ,0.  ,0.  ],
                                      [0.25,0.25,0.25] ])       

    prim_lattice_vectors = 0.5*( np.array( [ [1.,  1.,  0.],
                                             [0.,  1.,  1.],
                                             [1.,  0.,  1.] ]) )  
                            
    #High symmetry points in units of 2 pi/al (reduced crystalographic)
    k_units='reduced'
    X=np.array([0.,1.,0.])
    T=np.array([0.,0.,0.])
    L=np.array([0.5,0.5,0.5])
    U=np.array([0.25,1.,0.25])
    K=np.array([0.75,0.75,0.])
    W=np.array([0.5,1.,0.])
    #Corresponding dictionary
    HSdic={'X':X,'T':T,'L':L,'U':U,'K':K,'W':W}

    #Good reference for lattice constants
    #http://7id.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html

    #Lattice constant (Angstrom - default for geo format)
    al= {"C":3.56683,"ZnS":5.420,"Si":5.43095,"GaAs":5.6533}
    bond_length={"C":1.54}
    
    def bond_length2lattice_constant(bond_length):
        return 4.*bond_length/(3.**0.5)


    
#Assigns the lattice parameters and vectors based on the material, structure and basis choice   
class Lattice:
    def __init__(self,structure,basis,material):
        self.select_vectors(structure,basis,material)
        
    def select_vectors(self,structure,basis,material):
        #Cubic materials
        if (structure.lower()=='cubic' or structure.lower()=='zb') and basis[0:4].lower()=='prim':
            self.lattice_vector=Cubic.prim_lattice_vectors
            self.basis_vector=Cubic.prim_basis_vectors
            self.reciprocal_lattice_vector= self.reciprocal_lattice()
            if material not in Cubic.al:
                print(material,' lattice constant is not tabulated in dftb classes.')
                sys.exit('Script has stopped')
            else:
                self.al=Cubic.al[material]
            self.header=Cubic.header
            self.basis=basis

        #Other materials
        if structure.lower()!='cubic' or basis[0:4].lower()!='prim':
            print('This combination of basis and lattice vectors has not been coded up')
            print('Address the Vectors class')
            sys.exit('Code has stopped')
        
    #Compute reciprocal lattice vectors
    def reciprocal_lattice(self):
        a= self.lattice_vector
        b= np.zeros(shape=(3,3))
        mag= np.dot(a[0,:],  np.cross(a[1,:],a[2,:]))
        Lattice.check_magnitude(mag)
        b[0,:] = np.cross(a[1,:],a[2,:]) /mag
        b[1,:] = np.cross(a[2,:],a[0,:]) /mag
        b[2,:] = np.cross(a[0,:],a[1,:]) /mag
        return b

    def check_magnitude(mag):
        if mag == 0:
            print('Magnitude of a1.(a2^a3) = 0')
            print('Issue with lattice vectors')
            sys.exit('Code has stopped')


#------------------------------------------------------------------------------------
# All information required by the Geometry section of the dftb_in.hsd
#---------------------------------------------------------------------------------------         
# Finite system.  
# Inputs:
#         elements[1:Natoms]        List of elements associated with each atomic position
#         material                  Material name 
#         boundary condition       'c'luster
#         position[1:3,1:Natoms]    Atomic positions (Ang should be used)
# Outputs: 
#         Natoms                    Number of atoms in system         
#         uni_elements              List of unique elements/atomic species in system
#         uni_elements_index        Dictionary pairing unique element with unique index            
#
# Periodic system. 
# Inputs:
#         elements[1:basis]         List of elements consistent with basis atoms          
#         material                  Material name 
#         boundary condition       's'upercell or 'f'ractional
#                                   => calculation is either for bulk or supercell
#         position[1:3,1:Natoms]    Atomic positions in supercell (Ang should be used)
#         basis_vectors[?,?]        Positions of basis atoms for unit cell
#         Lattice vectors[3,3]      Defines unit cell or supercell 
# Outputs: 
#         Natoms                    Number of atoms in supercell
#         uni_elements              List of unique elements in supercell
#                                   when system==bulk, elements==uni_elements
#         uni_elements_index        Dictionary pairing unique element with unique index 
#---------------------------------------------------------------------------------------
class Geometry(object):
    def __init__(self,material,elements,boundary_conditions,\
                      position=None, al=None, basis_vectors=None, lattice_vectors=None, header=None ):

        self.material = material 
        self.elements = elements
        self.boundary_conditions = boundary_conditions.lower()
        self.Natoms = len(elements)
        self.position=position
        self.al= al
        self.basis_vectors = basis_vectors
        self.lattice_vectors = lattice_vectors
        self.header=header

        #Should in principle add these to the argument list too
        uni_elements=[]
        uni_elements_index={}
        cnt=0
        for ele in elements:
            if ele not in uni_elements:
                uni_elements.append(ele)
                cnt+=1
                uni_elements_index.update({ele:cnt})
                
        self.uni_elements = uni_elements
        self.uni_elements_index=uni_elements_index

        #Flags 
        if self.boundary_conditions=='c' or self.boundary_conditions=='s':
            if len(elements) != len(position):
                print("Elements list must correspond to position list")
                #Assume dimensions for position[1:Natoms,1:3]
                if len(position[0])> len(position):
                    print('DFTB geometry class expects dimensions for position[1:Natoms,1:3]')
                sys.exit('Script has stopped')

        if self.boundary_conditions=='s' or  self.boundary_conditions=='f':
            if al==None:
                print('Note, using using periodic boundary conditions but not defined lattice constant in DFTB+ object')
            #Probably a unit cell
            if self.Natoms<=10 and basis_vectors == None:
                print("Require basis positions for bulk calculations")
                sys.exit('DFTB class has stopped')       
            if lattice_vectors.all()==None:
                print("Require lattice vectors for periodic bulk and supercell calculations")
                sys.exit('DFTB class has stopped')
               
        elif self.boundary_conditions not in ['c','s','f']:
            print("Choice of boundary condition,",boundary_conditions,", for DFTB+ is not valid")
            sys.exit('DFTB class has stopped')
            


class Hamiltonian:
    def __init__(self,scc,slaterkosterfiles,max_ang_momentum,scc_tolerance=1e-5,tmp_f=0.,tmp_unit='Kelvin'):
        self.scc =scc
        self.scc_tolerance = scc_tolerance
        self.slaterkosterfiles=slaterkosterfiles
        self.max_ang_momentum=max_ang_momentum
        self.tmp_f=tmp_f
        self.tmp_unit=tmp_unit
        
class SlaterKosterFiles:
    def __init__(self,prefix,separator="-",suffix=".skf"):
        if prefix[-1] != '/':
           prefix = prefix+'/'
        self.prefix=prefix
        self.separator = separator
        self.suffix=suffix

        
#Have the 3 optional arguments and select appropriately 
class Ksampling:
    def __init__(self,supercellfolding=None, explicitkpoints=None, klines=None):
        self.supercellfolding=supercellfolding
        self.explicitkpoints=explicitkpoints
        self.klines=klines 
        Ksampling.check_conflicting_options(supercellfolding,explicitkpoints,klines)
        
    def check_conflicting_options(supercellfolding,explicitkpoints,klines):
        conflict=False
        
        if supercellfolding != None:
            if explicitkpoints != None or klines != None:
                conflict=True
        if explicitkpoints != None:
            if supercellfolding != None or klines != None:
                conflict=True
        if klines != None:
            if supercellfolding != None or explicitkpoints != None:
                conflict=True

        if conflict == True:
            print('Conflict: More than one k-grid option object is set:')
            print('  supercellfolding = ',supercellfolding)
            print('  explicitkpoints = ',explicitkpoints)
            print('  klines = ',klines)
            sys.exit('Code has stopped')
                  

    
class ExplicitKPoints:
    def __init__(self):
        print('Need to write class defining explicit k-points')

        
#-------------------------------------------------------
# Class containing functions for use with k-vectors
#-------------------------------------------------------
class KFunctions(object):
    def __init__(self):
        pass

    #------------------------------------------------------------------
    #Convert reduced kpoint coordinates (units of 2pi/a) to fractional
    #coordinates, where b are the primitive reciprocal space lattice
    #vectors and G are k-vectors  
    #------------------------------------------------------------------
    def reduced_to_fractional(self,b,G,basis,printt='no'):
        if basis[0:4].lower() !='prim':
            print('Warning, lattice vectors provided must be defined for the primitive cell \
                   to transform from absolute to fractional coordinates.')
            sys.exit('Code has stopped')
        
        #High symmetry points MUST be stored columnwise, in a matrix of (3,N)
        G=G.transpose() 
        f= np.empty(shape=(3,len(G)))
        #(b)^-1 i.e. the inverse of the reciprocal lattice vectors gives real space vectors
        t=np.linalg.inv(b.transpose())    
        if printt =='yes':
            print('transformation matrix:')
            print(t)
            print('K-point:')
            print(G)
        f= np.dot(t,G)
        #Fractional coordinates returned columnwise, hence transpose to rowwise
        f=f.T
        return f


    def convert_kvector_units(self,lattice,kin_unit,kout_unit):
        if kin_unit=='reduced':
                if kout_unit=='reduced':
                     return self.HS_points
                if kout_unit=='absolute':
                     self.HS_points=(2.*math.pi/lattice.al)*(self.HS_points)
                     return
                if kout_unit=='fractional':
                     self.HS_points=self.reduced_to_fractional(lattice.reciprocal_lattice_vector, \
                                                                     self.HS_points,lattice.basis,printt='no')
                     return
        if kin_unit=='absolute':
            if kout_unit=='reduced' or  kout_unit=='fractional':
                print('Conversion not written for: ',kin_unit,':',kout_unit)
                sys.exit('Code has stopped')

        if kin_unit=='fractional':
            if kout_unit=='absolute' or  kout_unit=='fractional':
                print('Conversion not written for: ',kin_unit,':',kout_unit)
                sys.exit('Code has stopped')


                
#High symmetry points that form k-vectors, number of sampling points per k-vector
#and units used for input k-points
class KLines(KFunctions):
    def __init__(self,lattice,HS_points_list,HS_points,points_per_line,kin_unit,kout_unit='fractional'):
        super(KLines,self).__init__()

        self.HS_points_list=HS_points_list
        self.HS_points = self.check_internal_kvector_units(lattice,HS_points,kin_unit)
        self.Nkpts_per_line = points_per_line
        self.kin_unit = kin_unit
        self.kout_unit = kout_unit

    #Internal units of kvectors should always be 'reduced'
    def check_internal_kvector_units(self,lattice,HS_points,kin_unit):
            if kin_unit=='reduced':
                    return HS_points
            elif kin_unit=='absolute':
                    HS_points = KLines.convert_kvector_units(self,lattice,kin_unit,kout_unit='reduced')
                    return HS_points
            elif kin_unit=='fractional':
                    HS_points = KLines.convert_kvector_units(self,lattice,kin_unit,kout_unit='reduced')
                    return HS_points
            else:
                print('Input units for kvectors not recognised: ',kin_unit)
                sys.exit('Code has stopped')
  

def HS_labels_2_points(HS_points_list):
   HS_points=np.empty(shape=(len(HS_points_list),3))
   i=0
   for point in HS_points_list:
      HS_points[i,:]=Cubic.HSdic[point]
      i=i+1
   return HS_points

        
class SuperCellFolding:
   #Reasonable sampling for Monkhorst-Pack grid 
   default_grid = np.array([[12., 0., 0.],[0., 12., 0.], [0., 0., 12.]])
   default_sweight = np.array([.5, .5, .5])

   def check_grid_validity(structure):
         if structure.lower() != 'zb'and structure.lower() !='cubic' :
               print('MP k-grid weights only valid for ZB structure, not ',structure)
               sys.exit('Script has stopped')

   #Initialise object of this class
   def __init__(self,structure, sweight=default_sweight,kgrid=default_grid):
         SuperCellFolding.check_grid_validity(structure)
         
         #Reasonable sampling for Monkhorst-Pack grid 
         default_grid = np.array([[12, 0, 0],[0, 12, 0], [0, 0, 12]])
         
         self.kgrid=kgrid.astype(int)
         if (kgrid == default_grid).all() == True:
             self.sweight=sweight
         if (kgrid ==default_grid).all() == False:
             self.sweight=SuperCellFolding.init_sweights(kgrid)

             
   #Valid for diagonal grid for cubic materials 
   def init_sweights(kgrid):
        sweight=np.empty(shape=(3))
        
        if kgrid[0,0] % 2 == 0:
            sweight[0]=0.5
        if kgrid[0,0] % 2 != 0:
            sweight[0]=0.
            
        if kgrid[1,1] % 2 == 0:
            sweight[1]=0.5
        if kgrid[1,1] % 2 != 0:
            sweight[1]=0.
            
        if kgrid[2,2] % 2 == 0:
            sweight[2]=0.5
        if kgrid[2,2] % 2 != 0:
            sweight[2]=0.

        return sweight


class ConjugateGradient:
      def __init__(self,max_force_component=1.e-4, max_steps=100, output_prefix='geom.out', moved_atoms='1:-1'):
        self.max_force_component = max_force_component
        self.max_steps = max_steps
        self.output_prefix = output_prefix
        self.moved_atoms = moved_atoms

    
class Analysis:
    def __init__(self,output_DOS=False,shell_resolved=True):
        self.output_DOS=output_DOS
        self.shell_resolved=shell_resolved

        
class ParserOptions:
    def __init__(self,version=5):
        self.version=version


        









    
