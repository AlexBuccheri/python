"""
Set GW input
"""
from typing import Optional

class GWInput:

    # All GW input options
    gw_string_template = \
        """  
           <gw
            taskname="{taskname}"
            nempty="{nempty}"
            ngridq="{ngridq}"
            skipgnd="{skipgnd}"
            >

            <mixbasis
              lmaxmb="{lmaxmb}"
              epsmb="{epsmb}"
              gmb="{gmb}"
            ></mixbasis>

            <freqgrid
              nomeg="{nomeg}"
              freqmax="{freqmax}"
            ></freqgrid>

            <barecoul
              pwm="{pwm}"
              stctol="{stctol}"
              barcevtol="{barcevtol}"
            ></barecoul>

            <selfenergy
              actype="{actype}"
              singularity="{singularity}"
            ></selfenergy>

           </gw>
        """

    def __init__(self, taskname: str, nempty: int, ngridq: list, skipgnd: bool,
                 n_omega: int, freqmax: Optional[float] = 1.0,
                 lmaxmb: Optional[int] = 4, epsmb: Optional[float] = 1.e-3, gmb: Optional[float] = 1.0,
                 pwm: Optional[float] = 2.0, stctol: Optional[float] = 1.e-16, barcevtol: Optional[float] = 0.1,
                 actype: Optional[str] = 'pade', singularity: Optional[str] = 'mpb'
                 ):
        """
        Set GW input class
        :param taskname: GW method
        :param nempty: Number of unoccupied states
        :param ngridq: k-grid
        :param skipgnd: Skip ground state calculation
        :param n_omega: Number of imaginary frequency grid points
        """
        self.taskname = taskname
        self.nempty = nempty
        self.ngridq = ngridq
        self.skipgnd = skipgnd
        self.nomeg = n_omega
        self.freqmax = freqmax
        self.lmaxmb = lmaxmb
        self.epsmb = epsmb
        self.gmb = gmb
        self.pwm = pwm
        self.stctol = stctol
        self.barcevtol = barcevtol
        self.actype = actype
        self.singularity = singularity
        self.string = ''
        self.set_input()

    def dict_for_format(self):
        """
        Dictionary for use with .format
        Requires special handling of bools and lists
        :return:
        """
        d = {}
        for key, value in self.__dict__.items():
            if type(value) == bool:
                d[key] = str(value).lower()
            elif type(value) == list:
                d[key] = str(value)[1:-1].replace(",", "")
            else:
                d[key] = value
        return d

    def set_input(self):
        """
        Set input string for GW, only
        :return: self.input_string
        """
        self.string = self.gw_string_template.format(**self.dict_for_format())


def set_gw_input_string(gs_input: str, gw_input: GWInput, gw_template):
    """

    Given a converged ground state input, set it
    to repeat the ground state calculation from file (due
    to the additions to the basis) and add the GW inputs

    Note, both replace and format are not inplace

    :return: GW calculation input string
    """

    gs_input = gs_input.replace('do="skip"', 'do="fromfile"')
    gw_input = gw_template.format(**gw_input.dict_for_format())

    return gs_input.format(GW_INPUT=gw_input)
