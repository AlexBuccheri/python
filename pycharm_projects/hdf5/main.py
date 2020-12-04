"""
Look at parsing with HDF5

pip3 install h5py pyhdfview h5json h5sparse xmltodict
This would probably make more sense in a python notebook


https://www.pythonforthelab.com/blog/how-to-use-hdf5-files-in-python/
"""
import numpy as np
import h5py
import xmltodict  # Doesn't do a great job
from collections import OrderedDict
from pathlib import Path
import os

def read_exciting_xml_output(data_root=None, file_name=None) -> OrderedDict:
    """
    Parse Exciting XML output into a dictionary
    """
    if data_root is None:
        data_root = "/users/sol/abuccheri/tutorials/exciting/1.diamond/"
    if file_name is None:
        file_name = "info.xml"

    with open(data_root + file_name) as fd:
        info = xmltodict.parse(fd.read())
    fd.close()

    return info


def write_hdf5_datasets(hdf5_fname='info.hdf5'):
    """
    Write various data to hdf5 via dataset
    :return: None
    """

    # HDF5 complains if file exists. 'w'' write option didn't seem to deal with it

    if Path(hdf5_fname).is_file():
        os.system("rm " + hdf5_fname)

    hdf5_file = h5py.File(hdf5_fname, 'w')

    # Cannot pass python dictionary directly to hdf5
    # hdf5_file.create_dataset('volumes', data=volumes)

    hdf5_file.close()
    return


# Main Routine
info = read_exciting_xml_output()
ground_state = info["info"]['groundstate']['scl']
structure = ground_state['structure']

# Inspect the data associated with the crystal tag
for key, value in structure['crystal'].items():
    if key[0] == '@':
        key = key[1:]
    print(key, value)

# Probably want to do one loop over the XML keys to create the groups from keys
# Then do a second pass and fill the datasets into the groups as scalar types or np types
volumes = {'unit_cell': structure['crystal']['@unitCellVolume'],
           'brillouin_zone': structure['crystal']['@BrillouinZoneVolume']
           }








# GROUPS
# Groups are the container mechanism by which HDF5 files are organized.
# Key = group member name
# Value = group member
# crystal_group = hdf5_file.create_group('crystal')
# volume_group = crystal_group.create_group('volume')
