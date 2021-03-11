import pytest

from lorecommendations_parser import parse_lorecommendations
from parse_linengy import parse_lo_linear_energies

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

    n_nodes = 21
    for radial_solutions in basis_per_species['zr']:
        assert len(radial_solutions) == n_nodes, "Zr: Each l-channel contains n_nodes radial functions"

    for radial_solutions in basis_per_species['o']:
        assert len(radial_solutions) == n_nodes, "O: Each l-channel contains n_nodes radial functions"



def test_parse_lo_linear_energies():

    linear_energies_species = parse_lo_linear_energies(file_path="parse/example_data")

    assert len(linear_energies_species) == 2, "Two species in LINENGY.OUT"
    species = [key for key in linear_energies_species.keys()]
    assert 'Zr' in species, "Expect Zr in species list"
    assert 'O' in species, "Expect Zr O species list"

    linear_energies_Zr = linear_energies_species['Zr']
    linear_energies_O = linear_energies_species['O']

    # Reference values found by inspecting the parse/example_data/LINENGY.OUT
    l_channels = [l_channel for l_channel in linear_energies_Zr.keys()]
    assert len(linear_energies_Zr) == 5, "5 l-channels for Zr"
    assert max(l_channels) == 4, "l_max = 4 for Zr"

    assert linear_energies_Zr == {0: [-1.39, -1.39, -1.39, -1.39],
                                  1: [-0.51, -0.51, -0.51, -0.51],
                                  2: [0.33, 0.33, 0.33, 0.33],
                                  3: [1.0, 1.0, 1.0, 1.0],
                                  4: [1.0, 1.0, 1.0, 1.0]
                                  }

    l_channels = [l_channel for l_channel in linear_energies_O.keys()]
    assert len(linear_energies_O) == 4, "4 l-channels for O"
    assert max(l_channels) == 3, "l_max = 3 for O"

    assert linear_energies_O == {0: [-0.04125, -0.04125, -0.04125, -0.04125],
                                 1: [0.1, 0.1, 0.1, 0.1],
                                 2: [1.0, 1.0, 1.0, 1.0],
                                 3: [1.0, 1.0, 1.0, 1.0]
                                 }

    return
