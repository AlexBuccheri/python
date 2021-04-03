"""
Process results for pure MPI scaling tests
"""
import numpy as np
import matplotlib.pyplot as plt

from parse.parse_gw import parse_gw_timings


def timing_units(units:str) -> tuple:
    if units == 'seconds':
        unit_conversion = 1.
    elif units == 'minutes':
        unit_conversion = 1./60.
    elif units == 'hours':
        unit_conversion = 1./3600.
    else:
        print('Choice of units not valids')
        quit()
    return unit_conversion


def process_mpi_results(scaling_root: str):

    time_unit = 'hours'
    save_plot = True
    units = timing_units(time_unit)

    # Nodes to use in scaling tests
    # waiting for 1 node to finish
    # /users/sol/abuccheri/gw_benchmarks/scaling/pure_mpi/set2/n_nodes_1
    mpi_processes_per_node = 36
    nodes = np.arange(2, 10+1)
    processes = nodes * mpi_processes_per_node

    # Report most expensive parts of the calculation
    # Sum of the time in various routines exceeds total calculation time,
    # suggesting that there’s some overlap of the timings i.e. dielectric function and calcminm
    t_calcminm = []
    t_calcselfc = []
    t_total = []

    for node_count in nodes:
        job_dir = scaling_root + '/n_nodes_' + str(node_count)

        print("Reading file:", job_dir)
        timings = parse_gw_timings(job_dir)

        t_calcminm.append(timings['calcminm'] * units)
        t_calcselfc.append(timings['calcselfc'] * units)
        t_total.append(timings['Total'] * units)


    # Plot scaling
    print('Single thread in all cases')
    print('Each marker corresponds to a node - ran in multiples of 36')
    print("MPI-parallelised over q-points only")
    print("q = [8, 8, 8] implies that the calculations did not reach one MPI process per q-point")
    print('Waiting for 1-node job to finish')

    fig, ax = plt.subplots()
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_xlabel('MPI Processes', fontsize=13)
    ax.set_ylabel('Time (' + str(time_unit) + ')', fontsize=13)

    ax.plot(processes, t_calcminm, color='blue', marker='o',
            linestyle='solid', linewidth=3, markersize=6, label='Expansion Coeffs M')
    ax.plot(processes, t_calcselfc, color='green', marker='o', linestyle='solid',
            linewidth=3, markersize=6, label='Re{Sigma_c}')
    ax.plot(processes, t_total, color='red', marker='o', linestyle='solid',
            linewidth=3, markersize=6, label='Total')

    ax.legend(loc='upper right')

    if save_plot:
        plt.savefig('pure_mpi_scaling.jpeg', dpi=300, facecolor='w', edgecolor='w',
                    orientation='portrait', transparent=True, bbox_inches=None, pad_inches=0.1)
    plt.show()

    return

process_mpi_results("/users/sol/abuccheri/gw_benchmarks/scaling/pure_mpi/set2")