import numpy as np

# Convert ASE containers into SPGLIB containers

def ase_to_spglib(ase_data):
    # Get ASE into spglib format. All tuples apart from atomic_numbers
    lattice = []
    for vector in ase_data.cell:
        lattice.append(tuple(vector))

    # ASE converts basis to cartesian, so convert back to fractional
    inv_lattice = np.linalg.inv(np.transpose(np.asarray(ase_data.cell)))

    basis = []
    for pos in ase_data.positions:
        frac_pos = np.matmul(inv_lattice, pos)
        basis.append(tuple(frac_pos))

    atomic_numbers = []
    for an in ase_data.numbers:
        atomic_numbers.append(an)

    # SGLIB data storage for a crystal is a tuple, ordered as:
    molecule = (lattice, basis, atomic_numbers)
    return molecule
