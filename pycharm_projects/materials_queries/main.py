import pandas as pd
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure

dfs: dict = pd.read_excel('ct9b00322_si_002.xlsx')

n_materials = 471
# Skip header
# for i in range(1, n_materials+1):
#     print(dfs['Material'][i], dfs['MP_ID'][i])

with MPRester("Q53K5A4N4Xkjt22X6o") as m:
    structure: Structure = m.get_structure_by_material_id(dfs['MP_ID'][1])
    print(structure.lattice)
    print(structure.species)
    print(structure.cart_coords)
