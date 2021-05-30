"""
Run CP2K's xTB periodic implementation,
write out the radius cut-off for the potential of the form sqrt(1/(r^2+c)) - 1/r
(requires write(*,*) ... added to ...)
and tabulate
"""

import subprocess

from modules.parameters.elements import an_to_symbol
from modules.fileio.substitution import safe_substitute

# We don't care about most of these settings:
# As long as the calculation calls the short-range potential we can extract
# the required information
# Homo-nuclear molecule in a box
xtb_base_str = """
&FORCE_EVAL
  &DFT
   LSD
    &QS
      METHOD xTB
      &xTB
         DO_EWALD T
         CHECK_ATOMIC_CHARGES F
     &END XTB
    &END QS
    &KPOINTS
       SCHEME  GAMMA
    &END
    &SCF
      EPS_SCF 1.e-7
      MAX_SCF 2
      &MIXING
          METHOD DIRECT_P_MIXING
          ALPHA   0.40
      &END
    &END SCF 
  &END DFT
  &SUBSYS
    &TOPOLOGY
      &GENERATE
        BONDLENGTH_MAX = 7
      &END GENERATE	
    &END TOPOLOGY
    &CELL
      ABC 10.0 10.0 10.0
    &END CELL
    &COORD
    $X1   5.000000    5.000000   5.000000  sys
    $X2   5.000000    5.000000   6.000000  sys
   &END COORD
  &END SUBSYS
&END FORCE_EVAL
&GLOBAL
  PROJECT element-box
  PRINT_LEVEL LOW
  RUN_TYPE ENERGY
&END GLOBAL
"""

cp2k_exe = "/home/alex/cp2k.exe/local/cp2k"
input = "homo_rcut.inp"
n_elements = 86


def generate_input(ispecies: int, input: str):
    species = an_to_symbol[ispecies]
    input_str = safe_substitute(xtb_base_str, {'X1': species, 'X2': species})
    fid = open(input, "w+")
    fid.write(input_str)
    fid.close()
    return


def run_cp2k(cp2k_exe: str, input: str, mode="sdbg") -> str:
    assert mode in ["sdbg"]
    cp2k_command = [cp2k_exe + "." + mode, input]
    try:
        cp2k_result = subprocess.check_output(cp2k_command, stderr=subprocess.STDOUT).decode("utf-8")
        return cp2k_result
    except subprocess.CalledProcessError as error:
        #print("subprocess error:", error.returncode, "found:", error.output)
        return "NULL"


def get_rcut(cp2k_output: str) -> float:
    """
    Based on the format I've chosen for printing
    Not fast but does the job
    """
    for line in cp2k_output.split('\n'):
        if "species" in line:
            rcut_a = float(line.split()[2])
            return rcut_a
    return 0


def tabulate_cutoffs():
    for ispecies in range(1, 3):
        generate_input(ispecies, input)
        cp2k_output = run_cp2k(cp2k_exe, input)
        rcut_a = get_rcut(cp2k_output)
        print(an_to_symbol[ispecies], rcut_a)


tabulate_cutoffs()


