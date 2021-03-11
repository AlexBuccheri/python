# More effort to do it properly
# import xml.etree.ElementTree as ET
# tree = ET.parse('parse/example_data/Zr.xml')
# root = tree.getroot()
#
# for child in root:
#     print(child.tag, child.attrib)

import re

def parse_basis_as_string(file_name:str):
    """
    Another disgusting parser, but does the job
    :param file_name: species.xml
    :return: basis string with tags that can be substituted with .format
    """

    fid = open(file_name, "r")
    file = fid.readlines()
    fid.close()

    basis_string = ''
    for line in file[:-3]:
        # If explicitly defined lo
        if 'custom l' in line:
            digits = re.findall(r'\d+', line)
            l_value = int(digits[0])
            # Add custom tag such that .format can be used to
            # insert optimised lo's
            if l_value > 0:
                leading_ws = re.match(r"\s*", line).group()
                basis_string += leading_ws + "{custom_l" + str(l_value-1) + "} \n\n"

        basis_string += line

    basis_string += leading_ws + "{custom_l" + str(l_value ) + "} \n\n"

    # Close tags
    for line in file[-3:]:
        basis_string += line

    assert ('</spdb>' in basis_string.splitlines()[-3:]), "basis xml tag not closed correctly"

    return basis_string
