"""
Test functions in post_process.py
"""

from gw_benchmark_outputs.post_process_utils import get_basis_labels

def test_get_basis_labels:


    # Get basis labels from dictionary to type (per element at least)
    # Differentiate between default basis and optimised basis

    settings= {'rgkmax': 8,
                         'l_max_values': [D([('zr', 3), ('o', 2)])],
                         'n_img_freq': 32,
                         'q_grid': [2, 2, 2],
                         'n_empty_ext': [2000],
                         'max_energy_cutoffs': max_energy_exts_set4_spd
                        }

    basis_labels = get_basis_labels(root: str, settings: dict, verbose=false)

