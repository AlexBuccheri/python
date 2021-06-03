"""
Plots
"""
# External libraries
import matplotlib.pyplot as plt
from typing import List
import re

# My modules
from units_and_constants.unit_conversions import ha_to_mev
from gw_benchmark_outputs.post_process_utils import parse_gw_results_two, get_basis_labels, combine_species_basis_labels, \
    n_local_orbitals, process_basis_numbers, sum_los_per_species, get_basis_labels_AHH


def change_in_qp_gap(delta_E_qp, max_energy_exts: list):
    """
    Print change in QP gap
    """
    assert delta_E_qp.size == len(max_energy_exts), 'Should be same number of QP energies as there are basis energy cutoffs'

    for ie, energy in enumerate(max_energy_exts[1:], start=1):
        change_in_delta_E_qp = (delta_E_qp[ie] - delta_E_qp[ie-1])
        print("Change in Delta E_QP from max energy param: ")
        print(max_energy_exts[ie-1].replace("\n", " "), "to ")
        print(energy.replace("\n", " "), ":")
        print(change_in_delta_E_qp, "(meV)")
    return


def combine_basis_label_wrt_species(basis_labels: List[dict]) -> List[str]:

    concat_basis_labels = []
    for entry in basis_labels:
        assert isinstance(entry, dict), "basis_labels must contain dictionaries"
        label_str = "".join(key.capitalize() + ':(' + value.rstrip() + ')\n'
                            for key, value in entry.items())
        concat_basis_labels.append(label_str)

    return concat_basis_labels


def n_local_orbitals(basis_labels: list) -> List[dict]:
    """
    Compute the number of LOs per species

    :param basis_labels:
    :return: n_basis_orbitals
    """
    n_los = []
    for basis_label in basis_labels:
        # Name better. LOs for this calculation
        stuff = {}
        for species, label in basis_label.items():
            stuff[species] = sum([int(s) for s in re.findall(r'\d+', label)])
        n_los.append(stuff.copy())

    return n_los


def plot_convergence(x, y, labels):

    fig, ax = plt.subplots()
    ax.set_xlabel('N LOs')
    ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)')

    ax.plot(x, y, color='blue', marker='o', markersize=8)
    for i, txt in enumerate(labels):
        ax.annotate(txt, (x[i], y[i]))

    plt.show()

    return


# def convergence_per_channel()


def gw_basis_convergence(root: str):

    info = """Converge QP gap w.r.t each l-channel LOs of Zr"""
    print(info)

    # Read in each one and plot separately
    channel_dirs = ['s_channel', 'p_channel', 'd_channel', 'f_channel']
    s_channel = ['i0', 'i1', 'i2', 'i3', 'i4']
    p_channel = ['i0', 'i1', 'i2', 'i3', 'i4']
    d_channel = ['i0', 'i1', 'i2', 'i3', 'i4']
    f_channel = ['i0', 'i1', 'i2', 'i3', 'i4']

    basis_labels_s = get_basis_labels_AHH(root + '/s_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in s_channel])
    n_los = n_local_orbitals(basis_labels_s)
    basis_labels_s = combine_basis_label_wrt_species(basis_labels_s)
    data_set6_s = parse_gw_results_two(root + '/s_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in s_channel])

    delta_E_qp_set6_s = data_set6_s['delta_E_qp'] * ha_to_mev

    print("s channel")
    print(delta_E_qp_set6_s)

    # x = NLOs for Zr, y = QP - KS and labels are the full LO basis labels
    plot_convergence([n['zr'] for n in n_los], delta_E_qp_set6_s, basis_labels_s)


    quit()

    data_set6_p = parse_gw_results_two(root + '/p_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in p_channel])
    data_set6_d = parse_gw_results_two(root + '/d_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in d_channel])
    data_set6_f = parse_gw_results_two(root + '/f_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in f_channel])

    delta_E_qp_set6_p = data_set6_p['delta_E_qp'] * ha_to_mev
    delta_E_qp_set6_d = data_set6_d['delta_E_qp'] * ha_to_mev
    delta_E_qp_set6_f = data_set6_f['delta_E_qp'] * ha_to_mev

    # Give some changes in QP gap w.r.t. calculations

    change_in_qp_gap(delta_E_qp_set6_s, s_channel)

    print("p channel")
    print(delta_E_qp_set6_p)
    change_in_qp_gap(delta_E_qp_set6_p, p_channel)

    print("d channel")
    print(delta_E_qp_set6_d)
    change_in_qp_gap(delta_E_qp_set6_d, d_channel)

    print("f channel")
    print(delta_E_qp_set6_f)
    change_in_qp_gap(delta_E_qp_set6_f, f_channel)

    # For each channel, plot
    # TODO Get this done this evening and show andris
    # QP - KS vs NLOs, and add basis labels
    # Note, also need to refactor get_basis_labels so I can just pass it the directory




    basis_labels = get_basis_labels(root, settings_set4_spd, verbose=True)
    basis_labels_set4 = combine_species_basis_labels(basis_labels, species_per_line=True)['(3,2)']
    n_los_set4 = n_local_orbitals(basis_labels)['(3,2)']
    # Total number of LOs per calculation (i.e. sum Zr and O basis sizes together)
    n_los_set4 = sum_los_per_species(n_los_set4)

    plot_convergence(x=nlos, y=delta_E_qp_set6_s, labels=basis_labels)



gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8")