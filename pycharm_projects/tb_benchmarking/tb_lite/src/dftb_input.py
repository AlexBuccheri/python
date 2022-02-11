"""
Generate DFTB+/TB lite input files
"""
import numpy as np
import ase


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
        ham_string = self.hamiltonian.to_string() if isinstance(self.hamiltonian, (Hamiltonian, BandStructureHamiltonian)) else ''
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


class KGrid:
    def __init__(self, k_grid: list, k_weights=None):
        """ k-grid sampling class
        :param k_grid: 3 integers
        """
        self.k_grid = k_grid
        self.k_weights = self.set_k_weights(k_weights)

    def set_k_weights(self, k_weights) -> list:
        if k_weights is not None:
            return k_weights
        weights = [0.0, 0.0, 0.0]
        # 0.0 if odd, 0.5 if even
        weights = [0.5 for k in self.k_grid if k % 2 == 0]
        return weights

    def to_string(self) -> str:
        k1, k2, k3 = self.k_grid
        w1, w2, w3 = self.k_weights
        string = f"""  KPointsAndWeights = SuperCellFolding {{
    {k1} 0 0
    0 {k2} 0
    0 0 {k3}
    {w1} {w2} {w3}
  }}
        """
        return string


class ExplicitKPoints:
    def __init__(self, k_points: np.ndarray, k_weights=None):
        """

        :param k_points:  Expect .shape = (n_k_points, 3)
        :param k_weights:
        """
        self.k_points = k_points
        self.n_k_points = k_points.shape[0]
        self.k_weights = self.set_k_weights(k_weights)

    def set_k_weights(self, k_weights) -> list:
        if k_weights is None:
            return [1.0] * self.n_k_points
        else:
            return k_weights

    def to_string(self) -> str:
        string = "  KPointsAndWeights = {\n"
        weight = "1.0"

        for ik in range(0, self.n_k_points):
            k_str = np.array2string(self.k_points[ik, :], precision=8, separator=' ', suppress_small=False)[1:-1]
            string += "    " + k_str + " " + weight + "\n"

        string += "    }\n"

        return string


class Hamiltonian:
    def __init__(self,
                 method='GFN1-xTB',
                 temperature=0.0,
                 scc_tolerance=1.e-6,
                 k_grid=None,
                 k_weights=None,
                 max_scf=50):
        if k_grid is None:
            k_grid = [4, 4, 4]
        self.method = method
        self.temperature = temperature
        self.scc_tolerance = scc_tolerance
        # Should initialise outside class and pass in object if being proper
        self.k_grid = KGrid(k_grid, k_weights)
        self.max_scf = max_scf

    def to_string(self):
        string = f"""Hamiltonian = xTB {{
  Method = "{self.method}"
  SCC = Yes
  SCCTolerance = {self.scc_tolerance}
  Filling = Fermi {{
    Temperature [Kelvin] = {self.temperature}
  }}
  {self.k_grid.to_string()}
}}
        """
        return string


class BandStructureHamiltonian:
    def __init__(self,
                 k_points,
                 method='GFN1-xTB',
                 k_weights=None,
                 ):
        self.k_points = ExplicitKPoints(k_points, k_weights)
        self.method = method

    def to_string(self):
        string = f"""Hamiltonian = xTB {{
  Method = "{self.method}"
  SCC = Yes
  ReadInitialCharges = Yes
  MaxSCCIterations = 1
  {self.k_points.to_string()}
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


def generate_band_structure_input(lattice_vectors, method: str) -> str:
    """ Generate DFTB+ input file string for a band structure calculation,
    using a band path as standardised by ASE.

    :param lattice_vectors: Crystal lattice vectors
    :param method: Calculation method
    :return: Input file string
    """
    assert method in ['GFN1-xTB', 'GFN2-xTB'], "Method is not valid"
    cell = ase.atoms.Cell(lattice_vectors)
    band_path = cell.bandpath()
    h_bands = BandStructureHamiltonian(band_path.kpts, method=method)
    return DftbInput(hamiltonian=h_bands).generate_dftb_hsd()
