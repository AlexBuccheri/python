#!/usr/bin/env python3

#------------------------------------------------------------
# Alexander Buccheri. UoB
# Classes and methods 
# Pattern matching and string replacement in files
#------------------------------------------------------------
import re
import sys
import os


#---------------------------------------------------------------------------------------
#Replace 'replace' between 'before' and 'after' in 'string' using regex
#---------------------------------------------------------------------------------------
def replace_string_section_regex(before,after,replace,string):
    if type(replace) is not str:
        replace=str(replace)
    reg = "(?<=%s).*?(?=%s)" % (before,after)
    r = re.compile(reg,flags=re.DOTALL | re.IGNORECASE)
    newstring = r.sub('  '+replace+'  ', string)
    return newstring 


#---------------------------------------------------------------------------------------
# Replace 'replace' between 'before' and 'after' in 'string'
# Match from the start of the line, rather than only one string
# backwards w.r.t. the string to replace.
# Only evaluating first character of 'after' string
# For example, if it's '#', but there's no space in '#comment', then the evaluation breaks.
#---------------------------------------------------------------------------------------
def replace_string_section(before,after,replace,padding,string):
    if type(replace) is not str:
        replace=str(replace)
    if  string.split()[0].lower() == before.lower()  \
    and string.split()[2].lower()[0] ==after.lower()[0]:
        newstring=before.title()+ padding +replace+ padding+after+" ".join(string.split()[3:])+'\n'
        return newstring
    else:
        return string



class String_replacement:
    def __init__(self,before,after,replace,padding):
        self.before=before
        self.after=after
        self.replace=replace
        self.padding=padding


#-------------------------------------------------
# Modify 1 or more specific strings in input file 
# fin='Input file path'
# fout='Output file path'
# string = Object of class String_replacement 
#-------------------------------------------------
def modify_variables(fin,fout,string):
    finput=fin
    foutput=fout
    
    #Iterate over N different strings to modify
    for i in range(0,len(string.before)):
        before=str(string.before[i])
        replace=str(string.replace[i])
        after=str(string.after[i])
        padding=str(string.padding[i])
        output_file_string=''
        
        with open(finput, 'r') as input_file:
            for line in input_file:
                if line.isspace() is False:
                    newline=replace_string_section(before,after,replace,padding, line)
                elif line.isspace() is True:
                    newline=line
                output_file_string += newline
            
        input_file.close()
        fid= open(foutput, "w")
        fid.write(output_file_string)
        fid.close()
        finput=foutput  


#Input is dir+filename 
def add_directory(fullpath): 
    directory=fullpath.split('/')[:-1]
    directory="/".join(directory)
    if not os.path.isdir(directory):
        os.makedirs(directory)
        print('Written directory:',directory)
