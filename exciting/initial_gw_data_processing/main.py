from parse.parse_gw import parse_gw_evalqp, parse_gw_info, parse_gw_timings

gw_data = parse_gw_info(file_path="./parse")
qp_data = parse_gw_evalqp("./parse", nempty=600, nkpts=3)
gw_timings = parse_gw_timings(file_path="./parse")




# qp - ks at gamma
process_gw_gamma_point(gw_data, qp_data)



# TODOs
# 	⁃	parse lorecommendations.
# 	⁃	    The main thing to do with avoiding low n states, is to know the nodal structure associated with the functions of a given l, already in the basis
# 	⁃	    See question to Andris ==> CORRECTLY SET UP OPTIMISED BASIS
# 	⁃	Label basis not with trial energies but with number of los per l-channel
# 	⁃	    See the downloaded paper and then ask
# 	⁃	Generate slurm or pbs input
# 	⁃	    Write my own or do with Aiida?
# 	⁃	GW input from ground state
# 	⁃	Parse ground state, change ground state fromscratch to fromfile, add GW options
# 	⁃	Some post-processing to get the data formats for plotting