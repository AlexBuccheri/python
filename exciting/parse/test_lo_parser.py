import pytest

from lorecommendations_parser import parse_lorecommendations

def test_parse_lorecommendations():
    """
    Test parsing lorecommendations file.

    RTF means if its structure changes, this will break.
    """

    species = ['zr', 'o']
    # Path assumes run from project root
    basis_per_species = parse_lorecommendations('parse/example_data/lorecommendations.dat', species)
    assert len(basis_per_species) == len(species), "One entry per species"

    l_max = 6
    n_l_channels = l_max + 1
    assert len(basis_per_species['zr']) == n_l_channels, "~r: All l-channels are stored per species"
    assert len(basis_per_species['o']) == n_l_channels, "O: All l-channels are stored per species"

    print(basis_per_species['zr'])

    n_nodes = 21
    for radial_solutions in basis_per_species['zr']:
        assert len(radial_solutions) == n_nodes, "Zr: Each l-channel contains n_nodes radial functions"

    for radial_solutions in basis_per_species['o']:
        assert len(radial_solutions) == n_nodes, "O: Each l-channel contains n_nodes radial functions"

