"""
Plot the QP energy w.r.t. systematic reduction the of LOs per channel,
such that a minimal basis for converged calculations can be found.

1 meV change in energy w.r.t converged basis from set9/A1_zr_o.py should
be the maximum allowed, as there are 13 channels => Could introduce 13 meV error
ontop of whatever the error is associated with converged basis (i.e. one might get
a few more meV from increasing rgkmax, moving to (7,6), increasing the energy cut-off etc.
"""