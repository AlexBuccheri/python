import json

# Given a JSON input file with data scraped from the Bilbao Crystalographic server,
# generatea python file containing all of the crystal Wyckoff positions, for all point groups.
#
# http://www.cryst.ehu.es/cryst/point_wp.html
#
# File naming convention according to shoenflies symbols.
#


# regex would be cleaner
def split_ops_string(ops_string):

    ops = []
    indices = [i for i, x in enumerate(ops_string) if x == "]"]
    i = 0
    for j in indices:
        ops.append(ops_string[i:j+2])
        i = j+2
    return ops


def function_definition_string(point_group):

    shoenflies = point_group["Shoenflies"].lower()
    ops_string = point_group["ops"]
    # Should just modify brackets n "extract_symmetry.py"
    ops_string = ops_string.translate(str.maketrans("{", "[") ) 
    ops_string = ops_string.translate(str.maketrans("}", "]") ) 
    ops = split_ops_string(ops_string)
    white_space = " " * 25

    function_str = "def wyckoff_positions_" + shoenflies + "(x, y, z): \n"
    function_str += "    positions = np.array(" + ops[0] + "\n"
    for i, op in enumerate(ops[1:]):
        function_str += white_space + op + "\n"

    function_str = function_str.rstrip() + ") \n"
    function_str += "    return positions \n\n"

    return  function_str 



# python file defining a function with switch statement 
def generate_crystal_point_groups(json_filename, output_dir=""):
    with open(json_filename) as fid:
        json_data = json.load(fid)

    import_string = "import numpy as np \nimport enum \n\n"

    # enum class
    enum_class_string = "class PointGroup(enum.Enum):\n"
    for index in json_data:
        point_group = json_data[index]
        shoenflies = point_group["Shoenflies"].lower()
        enum_class_string += "    " + shoenflies.capitalize() + " = enum.auto\n"
    enum_class_string += "\n"

    # Each function
    func_def_string = ''
    for index in json_data:
        point_group = json_data[index]
        func_def_string += function_definition_string(point_group)
    
    # Switch statement via dictionary 
    function_string = "def wyckoff_positions(point_group: PointGroup, point):\n"
    function_string += "    assert len(point) == 3\n"
    function_string += "    x = point[0]\n    y = point[1]\n    z = point[2]\n\n"
    function_string += "    switch = {\n"

    for index in json_data:
        point_group = json_data[index]
        shoenflies = point_group["Shoenflies"]
        enum_string = shoenflies.capitalize()
        function_name = "wyckoff_positions_" + shoenflies.lower()
        function_string += "        PointGroup." + enum_string + ": " + function_name +  ",\n"
    function_string +="        }\n"
    function_string +="    return switch.get(point_group)(point) \n\n"

    file_string = import_string + func_def_string + enum_class_string + function_string 
    fid = open(output_dir + "crystal_point_groups.py", 'w')
    fid.write(file_string)
    fid.close()

    return


# Run script
json_file = 'point_groups.json'
generate_crystal_point_groups(json_file)
