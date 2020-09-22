# Silicon cells for gamma point calculation with Grimme gfn1

import numpy as np

from modules.electronic_structure.structure import atoms, bravais, crystals, supercell
from modules.electronic_structure.basis import gfn1

from modules.parameters import crystals as params
from modules.fileio import write
from modules import entos

# Silicon cubic unit cell in angstrom
al = params.silicon['lattice_constant']['angstrom']
fractional_basis_positions = crystals.silicon('conventional')
lattice = bravais.simple_cubic(al)

unit_cell = []
for atom in fractional_basis_positions:
    pos_angstrom = np.matmul(lattice, atom.position)
    unit_cell.append(atoms.Atom(atom.species, pos_angstrom))


def silicon_supercell(n, centred_on_zero=None):
    if centred_on_zero == None:
        centred_on_zero = (n[0]*n[1]*n[2]) % 2 != 0

    translations = supercell.translation_vectors(lattice, n, centred_on_zero=centred_on_zero)
    super_cell = supercell.build_supercell(unit_cell, translations)

    # Information
    print('Number of atoms in supercell:', len(super_cell))

    # Find central cell
    icentre = supercell.flatten_supercell_limits([0, 0, 0], n, centred_on_zero=centred_on_zero)

    cells = supercell.list_global_atom_indices_per_cells(unit_cell, translations)
    central_cell_atom_indices = cells[icentre]

    print("Central cell contains atoms with this indices:", central_cell_atom_indices)
    central_cell = []
    for iatom in central_cell_atom_indices:
        central_cell.append(super_cell[iatom])

    return super_cell




#ns = [[1,1,1], [2,2,2], [3, 3, 3], [4,4,4], [5,5,5]]
ns = [[1,1,1]]
for n in ns:
    super_cell = silicon_supercell(n)
    lattice_opts = bravais.LatticeOpt(a=n[0]*al, b=n[1]*al, c=n[2]*al, alpha=90, beta=90, gamma=90, bravais_type='cubic')
    mole_string = write.turbomole_riper_periodic(super_cell, lattice_opts)
    fid = open("silicon"+str(len(super_cell))+".in", "w")
    print(mole_string, file=fid)
    fid.close()






