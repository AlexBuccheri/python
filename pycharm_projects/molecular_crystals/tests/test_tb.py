from tb_lite.src.tb import parse_qcore_structure, parse_qcore_settings


def test_parse_qcore_input():
    # Assumes pytest run from the root
    structure: dict = parse_qcore_structure('data/entos/acetic.cif.in')
    # TODO(Alex) Add an assert
    #print(structure)


def test_parse_qcore_settings():
    data = parse_qcore_settings('data/entos/acetic.cif.in')
    assert set(data.keys()) == {'k_points', 'electronic_temperature'}
    assert data['k_points'] == [4, 4, 4]
    assert data['electronic_temperature'] == 0.0


