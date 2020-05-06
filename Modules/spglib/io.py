
# Helper functions for printing, reading and writing of SPGLIB data
# https://github.com/atztogo/spglib/blob/master/python/examples/example.py

# Libraries
from ase.io import write as ase_write

#My modules
from modules.spglib.cell import spglib_to_ase


def show_symmetry(symmetry, n_symmetries=None):
    if n_symmetries == None:
        for i in range(symmetry['rotations'].shape[0]):
            print("  --------------- %4d ---------------" % (i + 1))
            rot = symmetry['rotations'][i]
            trans = symmetry['translations'][i]
            print("  rotation:")
            for x in rot:
                print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
            print("  translation:")
            print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))
    else:
        for i in range(0, n_symmetries):
            print("  --------------- %4d ---------------" % (i + 1))
            rot = symmetry['rotations'][i]
            trans = symmetry['translations'][i]
            print("  rotation:")
            for x in rot:
                print("     [%2d %2d %2d]" % (x[0], x[1], x[2]))
            print("  translation:")
            print("     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2]))


def show_lattice(lattice):
    print("Basis vectors:")
    for vec, axis in zip(lattice, ("a", "b", "c")):
        print("%s %10.5f %10.5f %10.5f" % (tuple(axis, ) + tuple(vec)))


def show_cell(lattice, positions, numbers):
    print("Number of atoms in primitive: ", len(numbers))
    show_lattice(lattice)
    print("Atomic points:")
    for p, s in zip(positions, numbers):
        print("%2d %10.5f %10.5f %10.5f" % ((s,) + tuple(p)))


# Not in my prefered notation:
# https://en.wikipedia.org/wiki/Crystallographic_point_group#Hermann–Mauguin_notation
def show_spg_symmetry_info(dataset, wyckoff=False, equivalent_atoms=False):
    print("  Spacegroup is %s (%d)." % (dataset['international'], dataset['number']))
    print("  Pointgroup is %s." % (dataset['pointgroup']))
    print("  Hall symbol is %s (%d)." % (dataset['hall'], dataset['hall_number']))
    if wyckoff:
        print("  Wyckoff letters are: ", dataset['wyckoffs'])
    if equivalent_atoms:
        print("  Mapping to equivalent atoms are: ")
        for i, x in enumerate(dataset['equivalent_atoms']):
            print("  %d -> %d" % (i + 1, x + 1))
    return


# A wrapper for ASE's write
def write(f_name, spg_molecule, pbc = None):
    ase_primitive_cell = spglib_to_ase(spg_molecule)

    if pbc:
        assert(isinstance(pbc, tuple))
        ase_primitive_cell.set_pbc(pbc)

    ase_write(f_name, ase_primitive_cell)
    return



