#!/usr/bin/env python3

#------------------------------------------------------------
# Alexander Buccheri. UoB March 2019
# Pattern matching and string replacement in files: 
# Find and captialise keywords in fortran code in the style
# specific to DLPOLY 
#------------------------------------------------------------
import re
import sys
import os


def first_character(line):
    return line.replace(" ", "")[0]

def contains_string(line):
    if ('\'' in line) or ('\"' in line):
        return True
    else:
        return False 

def capitalise_write(line):
    if line.strip()[0:5].lower() == 'write':
        #Only want to match first occurance i.e. write(
        string_column = r"{}".format('write\(')  
        line = re.sub(string_column, 'Write(', line)
    else:
        pass  
    return line
    
def process_file(fname,fkeywords,bkeywords):
    output_file_string = ''
    with open(fname, 'r') as input_file:
        for line in input_file:

            #Avoid comments, processing and strings 
            if first_character(line) != f90_comment and \
               first_character(line) != fpp_hash and \
               contains_string(line) == False:
                line = captialise_fkeywords(line,fkeywords)
                line = split_statementblock_keywords(line,bkeywords)
            #If write statement, capitalise first occurance of write()
            if contains_string(line):
                line = capitalise_write(line)
                
            output_file_string += line
    input_file.close()
    return output_file_string

def captialise_fkeywords(line,fkeywords):
    for keyword in fkeywords:
        string_column = r"\b{}\b".format(keyword)
        line = re.sub(string_column, keyword.capitalize(), line)
    return line
    
#elseif -> Else If, etc 
def split_statementblock_keywords(line,bkeywords):
    for key,value in bkeywords.items():
        string_column = r"\b{}\b".format(key)
        line = re.sub(string_column, value, line,flags=re.IGNORECASE)
    return line

bkeywords = {'elseif':'Else If','endif':'End If','enddo':'End Do'}

fkeywords = ['if','allocate','allocated','enddo','endif','present','then',
             'interface','module','subroutine','public','private','end',
             'intent','type','elseif','else','call','integer','real','logical',
             'character','deallocate','write','allocatable','do','use','size',
             'only','implicit','none','optional','class','read','inquire',
             'open','close']

f90_comment = '!'
fpp_hash='#'

# --------------------------
# Main Routine
# --------------------------

root='.'
fnames=['dftb_mpi.F90','dftb_api.F90','dftb_dummytypes.F90',
        'dftb_utils.F90']

for fname in fnames:
    fpath=root+'/'+fname
    file_string = process_file(fpath,fkeywords,bkeywords)
    nfname = root+'/'+str(fname.split('.')[0])+'_fmt.F90'
    fid= open(nfname, "w+")
    fid.write(file_string)
    fid.close()
