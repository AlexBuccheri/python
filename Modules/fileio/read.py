import typing
import numpy as np
import os

# Read  xyz format
def xyz(fname: str ) -> tuple:

    if not os.path.exists(fname):
        exit("File: "+fname+" does not exist")

    fid = open(fname, "r")
    data = fid.read()
    fid.close()

    natoms = int(data.splitlines()[0])
    header = data.splitlines()[1]
    names = []
    pos = []
    for line in data.splitlines()[2:]:
        names.append(line.split()[0])
        pos.append([float(x) for x in line.split()[1:]])

    assert(len(pos) == len(names))

    return (names, pos)



# Experiment with reading
# Read  xyz format
# def xyz(fname: str ) -> tuple:
#
#     if not os.path.exists(fname):
#         exit("File: "+fname+" does not exist")
#
#     fid = open(fname, "r")
#     data = fid.read()
#
#     natoms = int(data.splitlines()[0])
#     header = data.splitlines()[1]
#     # Remove natoms and header from data
#     # is this general?
#     n = len(data.splitlines()[0]) + len(data.splitlines()[1]) + 2 # \n
#     data = data[n:]
#     # Use fixed format of .xyz to extract names and positions separately
#     # Only keep first element per 4
#     names = data.split()[::4]
#     #Issue line - Delete first element .. deletes all?
#     pos = np.delete(x, np.arange(0, x.size, 1))
#
#     fid.close()
#     return