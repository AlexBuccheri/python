
# TODO(Alex)
# This is more work than I would like:
# Need to ensure the cell origin is at the corner
# Would prefer fractional coordinates
# Need to test that the expression works on a simple 2D example I can work out the result to
# dr - np.rint(dr / cell_lengths) * cell_lengths

# def minimum_image_distance_matrix(a: np.ndarray, lattice: np.ndarray):
#
#     assert a.shape[1] == 3, "Expect vectors to be Euclidean"
#     n_vectors = a.shape[0]
#     cell_lengths = np.asarray([np.linalg.norm(lattice[i, :]) for i in range(0, 3)])
#
#     d = np.empty(shape=(n_vectors, n_vectors))
#     np.fill_diagonal(d, 0.)
#
#     # Upper triangle
#     for i in range(0, n_vectors):
#         for j in range(i + 1, n_vectors):
#             dr = a[j, :] - a[i, :]
#             print(dr - np.rint(dr / cell_lengths) * cell_lengths)
#             d[i, j] = np.linalg.norm(dr - np.rint(dr / cell_lengths) * cell_lengths)
#             d[j, i] = d[i, j]
#
#     return d
#
#
# def distance_matrix(a: np.ndarray):
#
#     assert a.shape[1] == 3, "Expect vectors to be Euclidean"
#     n_vectors = a.shape[0]
#     d = np.empty(shape=(n_vectors, n_vectors))
#     np.fill_diagonal(d, 0.)
#
#     # Upper triangle
#     for i in range(0, n_vectors):
#         for j in range(i + 1, n_vectors):
#             d[i, j] = np.linalg.norm(a[j, :] - a[i, :])
#             d[j, i] = d[i, j]
#
#     return d
#
# def optimal_muffin_tin_radii(positions: np.ndarray, lattice: np.ndarray, atomic_numbers: List[int]):
#
#     # Use rgkmax as proxy for MT radius (seem to be correlated)
#     species_min_mt = np.amin([fixed_precision_rgkmax(x) for x in atomic_numbers])
#
#     dm = scipy_distance_matrix(positions, positions)
#     dm2 = distance_matrix(positions)
#     dm3 = minimum_image_distance_matrix(positions - positions[0, :], lattice)
#     print(dm - dm2)
#     print(dm2)
#     print(dm3)
#