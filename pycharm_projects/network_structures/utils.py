import numpy as np


# Convert ASE atom data into spglib data.
# All tuples apart from atomic_numbers
#
def ase_atom_to_spg_atom(ase_data):

    lattice = []
    for vector in ase_data.cell:
        lattice.append(tuple(vector))

    # ASE converts basis to cartesian, so convert back
    inv_lattice = np.linalg.inv(np.transpose(np.asarray(ase_data.cell)))

    basis = []
    for pos in ase_data.positions:
        frac_pos = np.matmul(inv_lattice, pos)
        basis.append(tuple(frac_pos))

    atomic_numbers = []
    for an in ase_data.numbers:
        atomic_numbers.append(an)

    molecule = (lattice, basis, atomic_numbers)

    return molecule