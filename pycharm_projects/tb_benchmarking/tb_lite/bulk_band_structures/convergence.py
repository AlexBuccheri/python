"""
Tabulated data for k-convergence.
All energies tabulated in eV
"""
# TODO(Alex) Turn this into JSON and generate automatically

silicon = [
    {'k_grid': [4, 4, 4], 'total_energy': -100.1538},
    {'k_grid': [6, 6, 6], 'total_energy': -100.1547},     # Converged
    {'k_grid': [8, 8, 8],    'total_energy': -100.1547},  # Use
    {'k_grid': [10, 10, 10], 'total_energy': -100.1547}
]

ge = [
    {'k_grid': [4, 4, 4], 'total_energy': -92.6272},
    {'k_grid': [6, 6, 6], 'total_energy': -92.6334},
    {'k_grid': [8, 8, 8],    'total_energy': -92.6338},
    {'k_grid': [10, 10, 10], 'total_energy': -92.6339},  # Converged
    {'k_grid': [12, 12, 12], 'total_energy': -92.6339}   # Use
]

diamond = [
    {'k_grid': [4, 4, 4], 'total_energy': -116.2287},
    {'k_grid': [6, 6, 6], 'total_energy': -116.2284},  # Converged
    {'k_grid': [8, 8, 8], 'total_energy': -116.2284}   # Use
]

# anatase = [
#     {'k_grid': [6, 6, 6], 'total_energy': },
#     {'k_grid': [8, 8, 8], 'total_energy': },
#     {'k_grid': [10, 10, 10], 'total_energy': },
#     {'k_grid': [12, 12, 12], 'total_energy': }
# ]
