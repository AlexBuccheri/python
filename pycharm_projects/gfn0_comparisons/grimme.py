elements = ['H',
'He',
'Li',
'Be',
'B',
'C',
'N',
'O',
'F',
'Ne',
'Na',
'Mg',
'Al',
'Si',
'P',
'S',
'Cl',
'Ar',
'K',
'Ca',
'Sc',
'Ti',
'V',
'Cr',
'Mn',
'Fe',
'Co',
'Ni',
'Cu',
'Zn',
'Ga',
'Ge',
'As',
'Se',
'Br',
'Kr']


#Make a bunch of xyz files, run grimme code for each

# Make xyz
for i,atom in enumerate(elements):
    fname = str(i)+'_'+atom+"_"+atom+".xyz"
    fid = open(fname, "w+")
    xyz_string = "2\n\n"+atom+" 0 0 0 \n"+atom+" 1 0 0"
    fid.write(xyz_string)
    fid.close()


# entos input file
# fname = "gfn0_homo_dimmers.xyz"
# fid = open(fname, "w+")
# string =''
# for atom in elements:
#    string += "xtb(structure(xyz=[[" + atom + ", 0,0,0], [" + atom + ",1,0,0]]) version=gfn0) \n"
# fid.write(string)
# fid.close()

#GFN run

