from modules.electronic_structure.grids import monkhorst_pack

#TODO(Alex) Move this file into its own project

# Original Monkhorst pack grid sampling
def print_mk_grid():
    print("Even integer sampling does not capture Gamma, odd integers do:")

    indices = monkhorst_pack.mk_function(1)
    print(indices)

    indices = monkhorst_pack.mk_function(2)
    print(indices)

    indices = monkhorst_pack.mk_function(3)
    print(indices)

    indices = monkhorst_pack.mk_function(4)
    print(indices)

# Monkhorst-pack shifted grid
# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.155109
# https://www.c2x.org.uk/mp_kpoints.html

i = 5

indices = monkhorst_pack.mk_function(i)
print(indices)

indices = monkhorst_pack.mk_function2(i)
print(indices)


# Other grid schemes