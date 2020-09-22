import numpy as np
import os

from modules.fileio import grep

# Test details
n_atoms = 568
n_md_steps = 1
dftb_prefix = 'standalone_dftb_'


def read_dl_forces(n_atoms, n_md_steps):
    dl_forces = np.zeros(shape=(3, n_atoms, n_md_steps))
    for imd in range(0, n_md_steps):
        fname = 'dl_forces_' + str(imd) + '.dat'
        dl_forces[:, :, imd] = np.loadtxt(fname, usecols=(1, 2, 3), skiprows=1, unpack=True)
    return dl_forces


def read_dftb_forces(n_atoms, n_md_steps):
    dftb_forces = np.zeros(shape=(3, n_atoms, n_md_steps))

    for imd in range(0, n_md_steps):
        forces_string = grep.extract('Total Forces', dftb_prefix + str(imd) + '/detailed.out',
                                     n_lines_after=n_atoms).splitlines()[1:]
        forces = []
        for line in forces_string:
            forces_line = line.split()[1:]
            forces.append([float(x) for x in forces_line])

        dftb_forces[:, :, imd] = np.transpose(np.asarray(forces))
    return dftb_forces


def compare_forces(dl_forces, dftb_forces):
    norm = np.zeros(shape=(n_atoms))
    for imd in range(0, n_md_steps):
        for iatom in range(0, n_atoms):
            diff = np.abs(dl_forces[:, iatom, 0] - dftb_forces[:, iatom])
            norm[iatom] = np.linalg.norm(diff)
        avg_err = np.sum(norm) / float(n_atoms)
        print("MD Step ", imd, ". Min, Max, Avg error:",
              np.amin(norm), np.amax(norm), avg_err)
    return

def compare_forces_for_step(dl_forces, dftb_forces):
    assert(dl_forces.shape[0] == dftb_forces.shape[0])
    assert(dl_forces.shape[1] == dftb_forces.shape[1])
    n_atoms = dl_forces.shape[1]
    norm = np.zeros(shape=(n_atoms))

    print("Min, Max, Avg error:")
    for iatom in range(0, n_atoms):
        diff = np.abs(dl_forces[:, iatom,] - dftb_forces[:, iatom])
        norm[iatom] = np.linalg.norm(diff)
    avg_err = np.sum(norm) / float(n_atoms)
    print(np.amin(norm), np.amax(norm), avg_err)
    return


def add_geometry_to_hsd(imd):
    geometry_string = 'Geometry = GenFormat {' + '\n'
    geometry_string += '  <<< "structure_' + str(imd) + '.gen"' + '\n'
    geometry_string += '}' + '\n\n'
    fid = open('dftb_in.hsd')
    hsd_string = ''
    for line in fid:
        hsd_string += line
    fid.close()
    return geometry_string + hsd_string


# DL_POLY should do this
def add_lattice_to_structure(imd):
    fid = open('structure_' + str(imd) + '.gen')
    structure_string = ''
    for line in fid:
        structure_string += line
    fid.close()
    structure_string += '\n\n'
    # Origin and lattice vectors for periodic systems
    # Extracted from CONFIG. Row-wise (ang) in both dl and DFTB inputs
    lattice_vectors = \
        ' 0.0000000000  0.0000000000   0.0000000000' + '\n' + \
        '50.0036047000  0.0000000000   0.0000000000' + '\n' + \
        '0.0000000000  19.500003000   0.0000000000' + '\n' + \
        '0.0000000000   0.000000000  67.5017716800'
    return structure_string + lattice_vectors


def set_up_dftb(n_md_steps):
    for imd in range(0, n_md_steps):
        dftb_dir = dftb_prefix + str(imd)
        os.mkdir(dftb_dir)

        hsd_string = add_geometry_to_hsd(imd)
        fid = open(dftb_dir + '/dftb_in.hsd', 'w')
        fid.write(hsd_string)
        fid.close()

        structure_string = add_lattice_to_structure(imd)
        fid = open(dftb_dir + '/structure_' + str(imd) + '.gen', 'w')
        fid.write(structure_string)
        fid.close()
    return


def run_dftb(n_md_steps):
    # Run DFTB+ instances
    for imd in range(0, n_md_steps):
        dftb_dir = dftb_prefix + str(imd)
        os.chdir(dftb_dir)
        os.system("export OMP_NUM_THREADS=1 && mpirun -np 2 $DFTB_ROOT/bin/dftb+")
        os.chdir('../')
    return


# Write me
def run_dl_poly(n_md_steps):
    return []
