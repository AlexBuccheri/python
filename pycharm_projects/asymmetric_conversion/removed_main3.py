quit()
# -------------------------------------------------------
# Index atoms in central cell
# -------------------------------------------------------

# Central cell index valid when translation vectors centred on zero
central_cell_index = int(0.5 * len(translations))
central_cell_atom_indices = np.arange(n_atoms_prim * central_cell_index,
                                      n_atoms_prim * (central_cell_index + 1))

check_central_cell = False
if check_central_cell:
    central_cell = []
    for ia in central_cell_atom_indices:
        central_cell.append(s_cell[ia])
    alex_xyz(output_dir + '/' + "aei_central_cell", central_cell)


# -----------------------------------------------------------------------
# Find neighbours to atoms of central cell that exist in adjacent cells
# -----------------------------------------------------------------------

def atoms_neighbouring_central_cell(central_cell_atom_indices, super_cell, radius, flatten=True):

    positions = [atom.position for atom in super_cell]
    d = spatial.distance_matrix(positions, positions)

    neighbours_of_central_cell = []
    for ia in central_cell_atom_indices:
        indices = np.where((d[ia, :] > 0.) & (d[ia, :] <= radius))[0]
        neighbours_of_central_cell.append(indices.tolist())

    if flatten:
        neighbours_of_central_cell = flatten_list(neighbours_of_central_cell)
        # Remove duplicates
        neighbours_of_central_cell = set(neighbours_of_central_cell)
        neighbours_of_central_cell = list(neighbours_of_central_cell)

    return neighbours_of_central_cell

def remove_central_cell_atoms(central_cell_atom_indices, neighbours_of_central_cell):

    neigbours_outside_cell = []
    for iN in neighbours_of_central_cell:
        if iN not in central_cell_atom_indices:
            neigbours_outside_cell.append(iN)

    return neigbours_outside_cell


neighbours_of_central_cell = atoms_neighbouring_central_cell(central_cell_atom_indices, s_cell, radius=1.8, flatten=True)
neighbours_outside_cell = remove_central_cell_atoms(central_cell_atom_indices, neighbours_of_central_cell)

neighbours = []
for iN in neighbours_outside_cell:
    position = s_cell[iN].position
    species = s_cell[iN].species
    neighbours.append(atoms.Atom(position=position, species=species))

unit_cell_and_neighbours = unit_cell + neighbours


# -----------------------------------------------------------------------
# Fold these atoms into central cell and see if they coincide with free-floating oxygens
# -----------------------------------------------------------------------

folded_neighbours = []

# Convert into fractional coordinates

# If so, then unit_cell_and_neighbours defines the set of atomic positions with which to work with

# Print out
#alex_xyz(output_dir + '/' + "aei_central_NN_cell", unit_cell_and_neighbours)



