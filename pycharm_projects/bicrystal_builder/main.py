
import os

# List of xyz data with no header or atom type (as all Cu)
def read_plane(fname: str):
    if not os.path.exists(fname):
        exit("File: "+fname+" does not exist")

    fid = open(fname, "r")
    data = fid.readlines()
    fid.close()

    positions = []
    for line in data:
        xyz = [float(i) for i in line.split()]
        positions.append(xyz)
    return positions


root = "/Users/alexanderbuccheri/Desktop/Work 2020/bicrystal/cu_planes/"
z = 6
positions = read_plane(root + "zplane" + str(z) + ".dat")
print(positions)
