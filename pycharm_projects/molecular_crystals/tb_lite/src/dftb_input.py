"""
Generate DFTB+/TB lite input file
"""
from typing import Optional


class DftbInput:
    """ Class to collate all DFTB+ input classes
    """

    def __init__(self,
                 driver=None,
                 hamiltonian=None,
                 options=None):
        self.driver = driver
        self.hamiltonian = hamiltonian
        self.options = options

    def generate_dftb_hsd(self) -> str:
        """ Generate dftb_in.hsd input string

        GEN format done separately, with ASE
        """
        driver_string = self.driver.to_string() if isinstance(self.driver, Driver) else ''
        ham_string = self.hamiltonian.to_string() if isinstance(self.hamiltonian, Hamiltonian) else ''
        options_string = self.options.to_string() if isinstance(self.options, Options) else ''

        dftb_in_template = f"""
    Geometry = GenFormat {{
        <<< "geometry.gen"
    }}

    {driver_string}

    {ham_string}

    {options_string}

    Parallel {{
      UseOmpThreads = Yes
    }}
        """
        return dftb_in_template


class Driver:
    def __init__(self,
                 type='ConjugateGradient',
                 lattice_option='No'):
        self.type = type
        self.lattice_option = lattice_option

    def to_string(self):
        string = f"""Driver = {self.type} {{
  LatticeOpt = {self.lattice_option}
}}"""
        return string


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

    def to_string(self):
        k1, k2, k3 = self.k_points
        w1, w2, w3 = self.k_weights
        string = f"""Hamiltonian = xTB {{
  Method = "{self.method}"
  SCC = Yes
  SCCTolerance = {self.scc_tolerance}
  Filling = Fermi {{
    Temperature [Kelvin] = {self.temperature}
  }}
  KPointsAndWeights = SuperCellFolding {{
    {k1} 0 0
    0 {k2} 0
    0 0 {k3}
    {w1} {w2} {w3}
  }}
}}
        """
        return string


class Options:
    def __init__(self, timing_level=1):
        self.timing_level = timing_level

    def to_string(self):
        string = f"""Options = {{
   TimingVerbosity = {self.timing_level}
}}"""
        return string
