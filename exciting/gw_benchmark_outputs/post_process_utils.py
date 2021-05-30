"""
Post-processing utilities that are common between various plot scripts
used for processing gw benchmark outputs
"""
from collections import OrderedDict
import numpy as np
import re

from process.optimised_basis import parse_species_string, create_lo_label
from parse.parse_gw import parse_gw_info, parse_gw_evalqp
from process.process_gw import process_gw_gamma_point

from gw_benchmark_inputs.input_utils import restructure_energy_cutoffs


def directory_to_string(l_max: OrderedDict) -> str:
    """
    Top level basis convergence directory string of the
    form:
        species1_lmax1_species2_lmax2_...speciesn_lmaxn
    where n = number of species.

    :return str, for example zr_4_o_3
    """
    path = ''
    for species, l in l_max.items():
        path += species + '_lmax' + str(l) + '_'

    return path


def l_max_dict_to_string(l_max: OrderedDict) -> str:
    """
    :param l_max:
    :return: str of the form (4,3)
    """
    l_max_str = '('
    for l in l_max.values():
        l_max_str +=  str(l) + ','
    return  l_max_str[:-1] + ')'


def max_energy_ext_per_directory(energy_cutoffs):
    """
    Get the energy extension string for each run directory,
    for a given set of calculations

    :param energy_cutoffs: of the form

    {'zr': {0: np.linspace(60, 120, num=4),
            1: np.linspace(60, 120, num=4),
            2: np.linspace(90, 300, num=4),
            3: np.linspace(60, 120, num=4)},
     'o':  {0: np.linspace(60, 120, num=4),
            1: np.linspace(60, 120, num=4),
            2: np.linspace(60, 120, num=4)}
                       }
    :return: max_energy_exts
    """
    max_energy_exts = []
    for energy in restructure_energy_cutoffs(energy_cutoffs):
        max_energy_per_species = [max(energy_per_l_channel.values()) for energy_per_l_channel in
                                  energy.values()]
        max_energy_exts.append( str(int(max(max_energy_per_species))) )
    return max_energy_exts


def parse_gw_results(root: str, settings: dict) -> dict:
    """

    QP direct-gap (relative to the KS gap)
    and self-energies of the band edges at Gamma.

    Almost easier to just path full directive strings.

    :return: dictionary containing the above.
    """

    # Basis settings
    rgkmax = settings['rgkmax']
    l_max_values = settings['l_max_values']
    max_energy_exts = settings['max_energy_cutoffs']

    # GW settings
    n_empty_ext = settings['n_empty_ext']
    q_grid = settings['q_grid']
    n_img_freq = settings['n_img_freq']

    # Data to parse and return
    delta_E_qp = np.empty(shape=(len(max_energy_exts), len(l_max_values)))
    re_self_energy_VBM = np.empty(shape=delta_E_qp.shape)
    re_self_energy_CBm = np.empty(shape=delta_E_qp.shape)
    q_str = "".join(str(q) for q in q_grid)

    # Lmax in LO basis
    for i, l_max in enumerate(l_max_values):
        basis_root = root + '/' + directory_to_string(l_max) + 'rgkmax' + str(rgkmax)
        gw_root = basis_root + "/gw_q" + q_str + "_omeg" + str(n_img_freq) + "_nempty" + str(n_empty_ext[i])

        # Max energy cut-off of LOs in each l-channel
        for ienergy, energy in enumerate(max_energy_exts):
            file_path = gw_root + '/max_energy_' + str(energy)

            gw_data = parse_gw_info(file_path)
            qp_data = parse_gw_evalqp(file_path)

            print('Reading data from ', file_path)
            results = process_gw_gamma_point(gw_data, qp_data)
            delta_E_qp[ienergy, i] = results['E_qp'] - results['E_ks']
            re_self_energy_VBM[ienergy, i] = results['re_sigma_VBM']
            re_self_energy_CBm[ienergy, i] = results['re_sigma_CBm']

    return {'delta_E_qp': delta_E_qp,
            're_self_energy_VBM': re_self_energy_VBM,
            're_self_energy_CBm': re_self_energy_CBm
            }


def get_basis_labels(root: str, settings: dict, verbose=False) -> dict:
    """
    For each calculation, with the directory hierarchy:
        zr_lmax4_o_lmax3_rgkmax7
            gw_q222_omeg32_nempty1000
                max_energy_ext

    parse the species basis files (xml) and generate a label for each.
    Return a dictionary with keys like:
      basis_labels['(4, 3)']['zr'] = lo_basis_label_strings

    :param root: str, root directory for calculations.
    :param settings: dict of calculation settings.
    :param verbose: bool, print names of parsed files.

    :return: A dictionary of the general form:
     basis_labels['(lmax_1, l_max2)'][species] = lo_basis_label_strings

    where len(lo_basis_label_strings) = len(max_energy_exts)
    """

    # Basis settings
    rgkmax = settings['rgkmax']
    l_max_values = settings['l_max_values']
    max_energy_exts = settings['max_energy_cutoffs']

    # GW settings that define directory names
    n_empty_ext = settings['n_empty_ext']
    n_img_freq = settings['n_img_freq']
    q_grid = settings['q_grid']
    q_str = "".join(str(q) for q in q_grid)

    basis_labels = OrderedDict()

    # Lmax in LO basis, for example {Zr:4, O:3}
    for i, l_maxs in enumerate(l_max_values):
        basis_root = root + '/' + directory_to_string(l_maxs) + 'rgkmax' + str(rgkmax)
        gw_root = basis_root + "/gw_q" + q_str + "_omeg" + str(n_img_freq) + "_nempty" + str(n_empty_ext[i])
        lmaxs_str = l_max_dict_to_string(l_maxs)
        basis_labels[lmaxs_str] = OrderedDict()

        # Max energy cut-off of LOs in each [species, l-channel]
        for species, l_max in l_maxs.items():
            basis_labels[lmaxs_str][species] = []
            l_values = [l for l in range(0, l_max + 1)]

            for ienergy, energy in enumerate(max_energy_exts):
                file_path = gw_root + '/max_energy_' + str(energy)
                if verbose:
                    print('Reading basis from ', file_path)

                fid = open(file_path + '/' + species.capitalize() + '.xml', 'r')
                basis_string = fid.read()
                fid.close()

                basis_los = parse_species_string(l_values, basis_string)
                basis_str = "".join(string for string in create_lo_label(basis_los))
                basis_labels[lmaxs_str][species].append(basis_str)

    return basis_labels


def combine_species_basis_labels(basis_labels: dict, species_per_line=False) -> dict:
    """
    Combine basis labels of all species, per energy cut-off entry/

    Convert from this form:
         basis_labels['(lmax_1, l_max2)'][species] = lo_basis_label_strings
    to this form:
         basis_labels['(lmax_1, l_max2)'] = lo_basis_label_strings

    where lo_basis_label_strings[energy_index] = 'Zr:(7s 7p 8d 7f). O:(6s 6p 6d)'

    :param basis_labels: dictionary of form: basis_labels['(lmax_1, l_max2)'][species] = lo_basis_label_strings
    :param species_per_line: bool, Put each species basis LOs on a new line. 
    :return: dictionary of form: basis_labels['(lmax_1, l_max2)'] = lo_basis_label_strings
    """

    new_line = '\n' if species_per_line else ''

    basis_per_lmaxpair = OrderedDict()
    for lmaxs, labels_by_species in basis_labels.items():
        basis_per_lmaxpair[lmaxs] = []

        # Syntax to select a sub-selection of an ordered dict is verbose
        first_entry = OrderedDict(list(labels_by_species.items())[:1])
        remaining_entries = OrderedDict(list(labels_by_species.items())[1:])

        # First species key
        for species, lo_labels in first_entry.items():
            for lo_label in lo_labels:
                combined_label = species.capitalize() + ':(' + lo_label.rstrip() + '). ' + new_line
                # Because we don't know what the size of basis_per_lmaxpair should be
                basis_per_lmaxpair[lmaxs].append(combined_label)

        # All other species keys
        for species, lo_labels in remaining_entries.items():
            for i, lo_label in enumerate(lo_labels):
                basis_per_lmaxpair[lmaxs][i] += species.capitalize() + ':(' +  lo_label.rstrip() +'). ' + new_line

    return basis_per_lmaxpair


def n_local_orbitals(basis_labels: dict) -> dict:
    """
    Compute the number of LOs per species
    TODO Alex. All this label stuff should be a class

    For example, for basis_labels =
    OrderedDict([('(4,3)', OrderedDict([ ('zr', ['7s 7p 8d 7f 6g ',
                                                 '8s 8p 11d 8f 7g ',
                                                 '9s 9p 13d 9f 8g ',
                                                 '9s 9p 15d 10f 9g ']),
                                         ('o', ['6s 6p 6d 5f ',
                                                '7s 7p 6d 6f ',
                                                '7s 8p 7d 7f ',
                                                '8s 8p 8d 7f '])
                                      ])
                )])

    the routine would return
    n_basis_orbitals = OrderedDict([('(4,3)', OrderedDict([
                                      ('zr', [35, 42, 48, 52]),  # where the digits of 7s 7p 8d 7f 6g sum to 35, etc.
                                      ('o', [23, 26, 29, 31])])
                                  )])

    :param basis_labels:
    :return: n_basis_orbitals
    """
    n_basis_orbitals = OrderedDict()

    for l_maxs_str, labels_by_species in basis_labels.items():
        n_basis_orbitals[l_maxs_str] = OrderedDict()

        # Label of LOs in each [species, l-channel]
        for species, lo_labels in labels_by_species.items():
            n_basis_orbitals[l_maxs_str][species] = []

            # Per LO energy cutoff
            for basis_label in lo_labels:
                #print(species, basis_label)
                n_los = sum([int(s) for s in re.findall(r'\d+', basis_label)])
                n_basis_orbitals[l_maxs_str][species].append(n_los)

    return n_basis_orbitals
