import numpy as np 

def off_diagonal_indices(n: int):
    """
    Get indices of all off-diagonal elements of a square matrix, with dimensions n. 

    Parameters
    ----------
    n : int
        Dimensions of square matrix

    Returns
    -------
    c : int array 
        Rototation matrix .shape(n,n)

    Notes
    -----
    Based on stack exchange example:
    https://stackoverflow.com/questions/35746806/
    how-to-get-indices-of-non-diagonal-elements-of-a-numpy-array
  
    """   
    return np.where(~np.eye(n, dtype=bool))
