import numpy as np

# Build a box of length a
# Put CO2 in centre of the box

cell_length = 50
cell_centre = [0.5*length for length in [cell_length,cell_length,cell_length]]
offset = 23
lattice_opts = LatticeOpt(bravais_type = 'cubic', a = cell_length)

supercell = []
for atom in carbon_dioxide:
    supercell.append(Atom(atom.type, atom.position + cell_centre + offset))

supercell = []
for atom in acetylene:
    supercell.append(Atom(atom.type, atom.position + cell_centre + offset))

structure_string = structure_string(supercell, lattice_opts)
print(structure_string)



#Cubic box with water in centre

water = atoms.water()
al = 50
cubic_lattice = bravais.simple_cubic(al)
cell_centre = [0.5 * length for length in [al, al, al]]

#Could make a builder which fills supercell with molecules with given offsets
#i.e. just adding to the list
supercell = []
for atom in water:
    supercell.append(atoms.Atom(atom.type, atom.position + cell_centre))


#xyz_file = xyz_string(supercell, "Water molecule in a box")

print(entos_structure_string(supercell))