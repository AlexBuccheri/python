from initial_gw_data_processing.zr_lmax3_o_lmax2.zr import zr_basis
from initial_gw_data_processing.zr_lmax3_o_lmax2.o import  o_basis

energy_cutoffs = [100.]   #[20., 40., 60., 80., 100]
max_matching_order = 1

for energy_cutoff in energy_cutoffs:
    print("ZR. Trial energy cutoff: " + str(energy_cutoff) + " max_matching_order:" + str(max_matching_order))
    print(zr_basis(energy_cutoff, max_matching_order))

    print("OXYGEN. Trial energy cutoff: " + str(energy_cutoff) + " max_matching_order:" + str(max_matching_order))
    print(o_basis(energy_cutoff, max_matching_order))
