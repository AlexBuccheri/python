"""
Generate DFTB+/TB lite input file
"""


class Driver:
    def __init__(self,
                 type='ConjugateGradient',
                 lattice_option='No'):
        self.type = type
        self.lattice_option = lattice_option


class Hamiltonian:
    def __init__(self,
                 method='GFN1-xTB',
                 temperature=0.0,
                 scc_tolerance=1.e-6,
                 k_points=None,
                 max_scf=50):
        self.method = method
        self.temperature = temperature
        self.scc_tolerance = scc_tolerance
        self.k_points = [4, 4, 4] if k_points is None else k_points
        self.k_weights = self.set_k_weights()
        self.max_scf = max_scf

    def set_k_weights(self) -> list:
        weights = [0.0, 0.0, 0.0]
        # 0.0 if odd, 0.5 if even
        weights = [0.5 for k in self.k_points if k % 2 == 0]
        return weights


class Options:
    def __init__(self, timing_level=1):
        self.timing_level = timing_level


def generate_dftb_hsd(driver: Driver, ham: Hamiltonian, options: Options) -> str:
    """ Generate dftb_in.hsd input string

    GEN format done separately, with ASE
    """
    k1, k2, k3 = ham.k_points
    w1, w2, w3 = ham.k_weights

    dftb_in_template = f"""
Geometry = GenFormat {{
    <<< "geometry.gen"
}}

Driver = {driver.type} {{
  LatticeOpt = {driver.lattice_option}
}}

Hamiltonian = xTB {{
  Method = "{ham.method}"
  SCC = Yes
  SCCTolerance = {ham.scc_tolerance}
  Filling = Fermi {{
    Temperature [Kelvin] = {ham.temperature}
  }}
  KPointsAndWeights = SuperCellFolding {{
    {k1} 0 0
    0 {k2} 0
    0 0 {k3}
    {w1} {w2} {w3}
  }}
}}

Options = {{
   TimingVerbosity = {options.timing_level}
}}

Parallel {{
  UseOmpThreads = Yes
}}
    """
    return dftb_in_template
