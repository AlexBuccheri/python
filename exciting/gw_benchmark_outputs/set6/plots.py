"""
Plots
"""
# External libraries
import matplotlib.pyplot as plt
from typing import List
import re

# My modules
from units_and_constants.unit_conversions import ha_to_mev
from gw_benchmark_outputs.post_process_utils import parse_gw_results_two, n_local_orbitals, get_basis_labels_AHH


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


def plot_convergence(x, y, labels, name: str, save_plot=False):
    """
    x = NLOs for Zr, y = QP - KS and labels are the full LO basis labels
    """
    fig, ax = plt.subplots()
    ax.set_xlabel('N LOs')
    ax.set_ylabel('Quasiparticle Gap - KS Gap at Gamma (meV)')

    ax.plot(x, y, color='blue', marker='o', markersize=8)
    for i, txt in enumerate(labels):
        ax.annotate(txt, (x[i], y[i]), fontsize=8)

    if save_plot:
        ax.set_aspect(0.7)
        #fig.subplots_adjust(bottom=0.15)
        #plt.subplots_adjust(right=0.5)
        #fig.tight_layout(rect=[1, 1, 3.0, 1])
        #fig.set_size_inches(10, 20)
        fig.tight_layout(pad=2)
        plt.savefig(name + '.jpg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)

    plt.show()

    return


def convergence_per_channel(path: str, directories: list):

    # Return basis in the form: basis_labels[directory_index] = {'species_A': basis_str, 'species_B': basis_str}
    basis_labels = get_basis_labels_AHH(path, directories)

    # Total number of LOs per species: n_los = [{'zr': 50, 'o': 45}]
    n_los = n_local_orbitals(basis_labels)

    # Combine species into one label: basis_labels[directory_index] = 'Zr:(5, 5, 10, 10)\n O:(6, 6, 6, 6)\n'
    combined_basis_labels = combine_basis_label_wrt_species(basis_labels)

    # Parse various GW data
    data_set6 = parse_gw_results_two(path, directories)

    # QP - KS in meV
    delta_E_qp = data_set6['delta_E_qp'] * ha_to_mev

    return {'x': [n['zr'] for n in n_los], 'y': delta_E_qp, 'labels': combined_basis_labels}


def gw_basis_convergence(root: str):

    info = """Converge QP gap w.r.t each l-channel LOs of Zr"""
    print(info)
    save_plot = False

    # Read in each one and plot separately
    s_channel = ['i0', 'i1', 'i2', 'i3', 'i4']
    p_channel = ['i0', 'i1', 'i2', 'i3', 'i4', 'i5', 'i6']
    d_channel = ['i0', 'i1', 'i2', 'i3', 'i4', 'i5', 'i6', 'i8']  # i7 failed
    #f_channel = ['i0', 'i1', 'i2', 'i3', 'i4', 'i5', 'i6']   # These runs were erroneous due to where the basis tag was
    f_channel = ['i0', 'i1', 'i3', 'i4']


    # print("s channel")
    # data_s = convergence_per_channel(root + '/s_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in s_channel])
    # print(data_s['x'], data_s['y'], data_s['labels'])
    # change_in_qp_gap(data_s['y'], s_channel)
    # plot_convergence(data_s['x'], data_s['y'], data_s['labels'], name='schannel_conv_lmax32', save_plot=save_plot)

    print("p channel")
    data_p = convergence_per_channel(root + '/p_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in p_channel])
    print(data_p['x'], data_p['y'], data_p['labels'])
    change_in_qp_gap(data_p['y'], p_channel)
    plot_convergence(data_p['x'], data_p['y'], data_p['labels'], name='pchannel_conv_lmax32', save_plot=save_plot)

    # print("d channel")
    # data_d = convergence_per_channel(root + '/d_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in d_channel])
    # print(data_d['x'], data_d['y'], data_d['labels'])
    # change_in_qp_gap(data_d['y'], d_channel)
    # plot_convergence(data_d['x'], data_d['y'], data_d['labels'], name='dchannel_conv_lmax32', save_plot=save_plot)

    # print("f channel")
    # data_f = convergence_per_channel(root + '/f_channel/gw_q222_omeg32_nempty2000', ['max_energy_' + ext for ext in f_channel])
    # print(data_f['x'], data_f['y'], data_f['labels'])
    # change_in_qp_gap(data_f['y'], f_channel)
    # plot_convergence(data_f['x'], data_f['y'], data_f['labels'], name='fchannel_conv_lmax32', save_plot=save_plot)


gw_basis_convergence("/users/sol/abuccheri/gw_benchmarks/A1_more_APW/set6/zr_lmax3_o_lmax2_rgkmax8")