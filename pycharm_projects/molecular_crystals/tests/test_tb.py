from tb_lite.src.tb import parse_qcore_structure


def test_parse_qcore_input():

    # Assumes pytest run from the root
    structure: dict = parse_qcore_structure('data/entos/acetic.cif.in')
    print(structure)




