import json

from excitingparser import ExcitingParser
from nomad.datamodel import EntryArchive


# from nomad.parsing.file_parser.xml_parser import XMLParser



def parse_exciting_results(directory: str):
    """
    Parse all exciting outputs in 'directory' into a dictionary, via
    NOMAD's exciting parser.

    :param directory:
    :return: dict
    """
    archive = EntryArchive()
    parser = ExcitingParser()
    parser.parse(directory + '/INFO.OUT', archive, None)
    data_as_json_str = archive.m_to_json()
    data_as_dict = json.loads(data_as_json_str)['section_run'][0]
    # dict_keys(['program_basis_set_type', 'program_name', 'program_version', 'section_frame_sequence', 'section_method', 'section_sampling_method', 'section_single_configuration_calculation', 'section_system'])
    print(data_as_dict)


def parse_dos_xml(directory: str):
    from excitingparser.exciting_parser import DOSXMLParser
    p = DOSXMLParser()
    p.mainfile = directory + '/dos.xml'
    print(p.keys())
    p.get('total')


def parse_info_out(directory: str):
    from excitingparser.exciting_parser import ExcitingInfoParser
    p = ExcitingInfoParser()
    p.mainfile = directory + '/INFO.OUT'
    print(p['initialization'].keys())
    print(p['groundstate'].keys())
    # Also works as
    p['initialization']
    print(p.get('initialization').lattice_vectors.magnitude)
    print(p.get('initialization').lattice_vectors.units)
    print(p.get('initialization').lattice_vectors.to('m'))





# def parse_gw_info(directory: str):
#     from excitingparser.exciting_parser import GWInfoParser
#     p = GWInfoParser()
#     p.mainfile = directory + '/GW_INFO.OUT'
#     print(p, vars(p))




parse_info_out('groundstate_example')
#parse_exciting_results('groundstate_example')