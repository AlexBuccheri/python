import pytest

from basis import parse_lo_linear_energies

def test_parse_lo_linear_energies():

    linear_energies_species =  parse_lo_linear_energies("process/example_data/LINENGY.OUT",
                                                        filter_duplicate_species=True)

    assert len(linear_energies_species) == 2, "Two species in LINENGY.OUT"

    species = [key for key in linear_energies_species.keys()]
    assert 'Zr' in species, "Expect Zr in species list"
    assert 'O' in species, "Expect Zr O species list"
    assert len(species) == 2

    linear_energies_Zr = linear_energies_species['Zr']
    linear_energies_O = linear_energies_species['O']

    # Order of l-channels doesn't matter

    assert linear_energies_Zr ==  {0: [-1.39, -1.39, -1.39, -1.39],
                                   1: [-0.51, -0.51, -0.51, -0.51],
                                   2: [0.33, 0.33, 0.33, 0.33],
                                   3: [1.0, 1.0, 1.0, 1.0]
                                   }
    # MISSING L=4 for linear_energies_Zr

    assert linear_energies_O ==  {0: [-0.04125, -0.04125, -0.04125, -0.04125],
                                  1: [0.1, 0.1, 0.1, 0.1],
                                  2: [1.0, 1.0, 1.0, 1.0],
                                  3: [1.0, 1.0, 1.0, 1.0]
                                  }



    return
