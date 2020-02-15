
from modules.electronic_structure.structure import atoms

def structure(molecule, unit):
    assert(type(unit) == atoms.CoordinateType)

    string = 'structure(  \n'
    
    if unit == atoms.CoordinateType.XYZ:
        prefix = "   xyz=[["
    elif unit == atoms.CoordinateType.FRACTIONAL:
        prefix = "   fractional=[["
        
    string += prefix + molecule[0].species + ', '
    pos = molecule[0].position
    string += str(pos[0]) + ', ' + str(pos[1]) + ', ' + str(pos[2]) + '],\n'

    prefix = '        ['
    for atom in molecule[1:]:
        pos = atom.position
        string += prefix + atom.species + ', ' + str(pos[0]) + ', ' + str(pos[1]) + ', ' + str(pos[2]) + '],\n'

    string = string[:-2] +'] \n )'
    return string