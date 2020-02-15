# Reference:
# https://stackoverflow.com/questions/35746806/how-to-get-indices-of-non-diagonal-elements-of-a-numpy-array
def off_diagonal_indices(n):
    return np.where(~np.eye(n, dtype=bool))
