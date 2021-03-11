"""




"""

import matplotlib.pyplot as plt

import gw_data



# Convergence in groundstate properties
# ADD ME
# Comment on rgkmax
# Briefly comment on basis definition
# Adding high-energy lo's - what does this mean? See some of my jupyter notes

# Using all unoccupied states. Varies per k-point (as atoms can have different basis sets?)
# COMMENT ON CHOICE OF UNOCCUPIED VALUES
# l_max => the max orbital angular momentum of the local orbitals in the muffin tin part of the basis
# 32 imaginary frequency points
# KS gap w.r.t.  max orbital angular momentum. doesn't change w.r.t. the GW options


# Quasiparticle Gap w.r.t l_max, for all energy parameters. Direct gap
# Note, 20 Ha looks erroneous for (3,2) and 100 Ha has not finished running
# Comment on a) what this energy parameter is and b) how I control it
fig, ax = plt.subplots()
gw_data.quasiparticle_gap_wrt_lmax_and_trial_energies(gw_data.gw_results, fig, ax, [40, 60, 80, 100])
legend = ax.legend(loc='lower right', title='Energy Param Cutoff (Ha)')
plt.title("Quasiparticle Gap vs MT Basis Size")
#plt.savefig('E_GW_vs_lmax.pdf', dpi=300)
plt.show()


# Quasiparticle Gap - KS Gap (recorded from the GW output) for 80 Ha. Direct gap
# fig, ax = plt.subplots()
# gw_data.quasiparticle_gap_minus_ks_gap(gw_data.gw_results, trial_energy=80, fig=fig, ax=ax)
# legend = ax.legend(loc='lower right', title='Energy Param Cutoff (Ha)')
# #plt.xlim((1,3))
# #plt.ylim((1.95,2.07))
# plt.show()

# Change in Quasiparticle Gap w.r.t l_max, for 80 Ha. Direct gap
fig, ax = plt.subplots()
gw_data.delta_quasiparticle_gap_wrt_lmax(gw_data.gw_results, trial_energy=100, fig=fig, ax=ax)
legend = ax.legend(loc='lower right', title='Energy Param Cutoff (Ha)')
#plt.savefig('delta_E_GW_vs_lmax.pdf', dpi=300)
plt.show()

data = gw_data.gw_results
energy100 = gw_data.trial_energy_to_index[100]
Eg_GW_43 = data[(4,3)]['E_GW_CBB'][energy100] - data[(4,3)]['E_GW_VBT'][energy100]
Eg_GW_54 = data[(5,4)]['E_GW_CBB'][energy100] - data[(5,4)]['E_GW_VBT'][energy100]
ha_to_mev = gw_data.ha_to_ev * 1000
print("For l_max values of (4,3), the quasiparticle gap is converged to:")
print((Eg_GW_54 - Eg_GW_43) * ha_to_mev)

# fig, ax = plt.subplots()
# gw_data.delta_quasiparticle_gap_wrt_lmax(gw_data.gw_results, trial_energy=80, fig=fig, ax=ax)
# legend = ax.legend(loc='lower right', title='Energy Param Cutoff (Ha)')
# plt.show()

# Real Self-energy VBM
fig, ax = plt.subplots()
gw_data.real_sigma_vbt_lmax_and_trial_energies(gw_data.gw_results, fig, ax, [40, 60, 80, 100])
legend = ax.legend(loc='upper right', title='Energy Param Cutoff (Ha)')
#plt.savefig('Re_sigma_VBM.pdf', dpi=300)
plt.show()

# Change is relatively large, however the Re(Self-energy) at VBM is not particularly sensitive to the energy cut-off
fig, ax = plt.subplots()
gw_data.real_sigma_vbt_lmax_and_trial_energies(gw_data.gw_results, fig, ax, [40, 60, 80, 100])
legend = ax.legend(loc='upper right', title='Energy Param Cutoff (Ha)')
plt.xlim((3, 3))
plt.ylim((3.8, 3.83))
#plt.savefig('Re_sigma_VBM_zoomed.pdf', dpi=300)
plt.show()

# Real Self-energy CBm
fig, ax = plt.subplots()
gw_data.real_sigma_cbb_lmax_and_trial_energies(gw_data.gw_results, fig, ax, [40, 60, 80, 100])
legend = ax.legend(loc='upper right', title='Energy Param Cutoff (Ha)')
#plt.savefig('Re_sigma_CBm.pdf', dpi=300)
plt.show()

# q-point convergence is running
# Add some details to the presentation