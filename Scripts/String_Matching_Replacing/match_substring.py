#!/usr/bin/env python3


import re
import sys
import os

#Regex function to extract a substring 
#https://stackoverflow.com/questions/4666973/how-to-extract-a-substring-from-inside-a-string-in-python
def find_string_section_regex(before,after,string):
    m = re.search(before+'(.+?)'+after, string)
    if m: substring = m.group(1)
    else: substring=''
    return substring
    
some_string="Volume:  12000 m^3"
substring=find_string_section_regex(':','m',some_string)
print(int(substring.lstrip()))
