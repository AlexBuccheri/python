import pytest
import numpy as np

from process.optimised_basis import parse_species_string, create_lo_label, latex_lo_basis_labels, \
    table_of_lo_energies

basis_string = """
 
<?xml version="1.0" encoding="UTF-8"?>
<spdb xsi:noNamespaceSchemaLocation="../../xml/species.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <sp chemicalSymbol="Zr" name="zirconium" z="-40.0000" mass="166291.1791">
    <muffinTin rmin="0.100000E-05" radius="2.0000" rinf="26.3269" radialmeshPoints="600"/>
    <atomicState n="1" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="2" l="1" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="0" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="3" l="1" kappa="1" occ="2.00000" core="true"/>
    <atomicState n="3" l="1" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="2" kappa="2" occ="4.00000" core="true"/>
    <atomicState n="3" l="2" kappa="3" occ="6.00000" core="true"/>
    <atomicState n="4" l="0" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="4" l="1" kappa="1" occ="2.00000" core="false"/>
    <atomicState n="4" l="1" kappa="2" occ="4.00000" core="false"/>
    <atomicState n="4" l="2" kappa="2" occ="2.00000" core="false"/>
    <atomicState n="5" l="0" kappa="1" occ="2.00000" core="false"/>
    <basis>

      <default type="lapw" trialEnergy="0.1500" searchE="false"/>
      <custom l="0" type="lapw" trialEnergy="-2.0" searchE="true"/>  
      <lo l="0"> 
	<wf matchingOrder="0" trialEnergy="-2.0" searchE="true"/> 
	<wf matchingOrder="1" trialEnergy="-2.0" searchE="true"/> 
      </lo> 
      <lo l="0"> 
	<wf matchingOrder="1" trialEnergy="-2.0" searchE="true"/> 
	<wf matchingOrder="2" trialEnergy="-2.0" searchE="true"/> 
      </lo> 

       <lo l="0">
        <wf matchingOrder="0" trialEnergy="3.42" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="3.42" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="11.86" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="11.86" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="23.56" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="23.56" searchE="false"/>
       </lo>
        
       <lo l="0">
        <wf matchingOrder="0" trialEnergy="38.22" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="38.22" searchE="false"/>
       </lo>
       
         <custom l="2" type="lapw" trialEnergy="0.33" searchE="true"/>  
      <lo l="2"> <wf matchingOrder="0" trialEnergy="0.33" searchE="true"/> <wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> </lo>
      <lo l="2"> 
    	<wf matchingOrder="1" trialEnergy="0.33" searchE="true"/> 
	    <wf matchingOrder="2" trialEnergy="0.33" searchE="true"/> 
      </lo> 

      
       <lo l="2">
        <wf matchingOrder="0" trialEnergy="5.63" searchE="false"/> 
        <wf matchingOrder="1" trialEnergy="5.63" searchE="false"/>
       </lo>


"""

def test_label_by_parsing_basis_string():
    basis_los = parse_species_string([0, 2], basis_string)
    basis_labels = create_lo_label(basis_los)
    basis_str = "".join(s for s in basis_labels)
    assert basis_str == '6s 3d ', "basis string"





zr_l0 = np.array([
 -956.668857441521,
 -90.4809096033025,
 -14.2958579803110,
 -1.39085408704669,
  3.41674973498977,
  11.8576335922119,
  23.5553488286487,
  38.2225139492148,
  55.7183934928493,
  75.9631434089540,
  98.9062300331486,
  124.506902544629,
  152.736321370449,
  183.572269369561,
  216.996197372533,
  252.993631242249,
  291.552350804872,
  332.661952981571,
  376.313557393332,
  422.499336330021,
  471.212424048644])

zr_l1 = np.array([
  -80.6169740080916,
  -11.0985082112073,
 -0.510197715196521,
   4.25309289714343,
   12.6683041144307,
   24.1783693122313,
   38.5543138792572,
   55.6730176346185,
   75.4659088967503,
   97.8905521589912,
   122.912605026692,
   150.509391934944,
   180.663611187383,
   213.361293673377,
   248.592066934076,
   286.347157691748,
   326.619318498780,
   369.402476756102,
   414.691420217843,
   462.481804168748,
   512.769999174939])

zr_l2 = np.array([
  -5.85379645274240,
  0.905888268806381,
   5.63201190642294,
   13.6742717607037,
   24.5438162077841,
   38.1197697668803,
   54.3281060520894,
   73.1368245214196,
   94.5250199257870,
   118.473406611852,
   144.970357053250,
   174.004502271232,
   205.566711556200,
   239.649662781681,
   276.246762591029,
   315.352828089062,
   356.963536377846,
   401.075359311671,
   447.685378181735,
   496.791009824918,
   548.389885069050])



def test_basis_table_generation():
    basis_los = parse_species_string([0, 2], basis_string)
    zr_lorecommendations = [zr_l0, zr_l1, zr_l2]
    latex_lo_basis_labels(basis_los, zr_lorecommendations)
    # TODO Add an assertion


def test_table_of_lo_energies():
    basis_los = parse_species_string([0, 2], basis_string)
    zr_lorecommendations = [zr_l0, zr_l1, zr_l2]
    table_of_lo_energies(basis_los, zr_lorecommendations)