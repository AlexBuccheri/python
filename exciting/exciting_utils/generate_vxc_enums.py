"""
Generate Vxc enums (parameters) module by scrubbing the ground state schema.


"""
import re

# XC functionals, only
start = 653
end = 2183

fid = open("/Users/alexanderbuccheri/Codes/exciting/xml/schema/groundstate.xsd")
file = fid.readlines()[start:end]
fid.close()

open_tag = False
vxc_enums = {}
for line in file:
    if ('<xs:enumeration' in line) or ('<!--xs:enumeration' in line):
        xc = re.findall('"([^"]*)"', line)[0]
        open_tag = True
    #if open_tag and ('oldnr' in line):
    if 'oldnr' in line:
        oldnr = int(re.findall(r'\d+', line)[0])
        # All 0 keys overwrite oneanother
        vxc_enums[oldnr] = xc
        open_tag = False

# Sort the keys
sorted_keys = [key for key in vxc_enums.keys()]
sorted_keys.sort()

for oldnr in sorted_keys:
    doc_string = '  !> ' + vxc_enums[oldnr] + ' index \n'
    enum_string = '  integer, public, parameter :: ' + vxc_enums[oldnr] + " = " + str(oldnr)
    print(doc_string + enum_string)



# print(oldnr)
# print(len(oldnr))
# print(len(set(oldnr)))
