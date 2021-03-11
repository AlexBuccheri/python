"""
Set GW input for A2 system
"""
from benchmark_inputs.A1_groundstate import converged_ground_state_input
from benchmark_inputs.A1_gw import gw_string_template

class GWInput:
    def __init__(self, taskname:str, nempty:int, ngridq:list, skipgnd:bool):
        """
        Set GW input class
        :param taskname: GW method
        :param nempty: Number of unoccupied states
        :param ngridq: k-grid
        :param skipgnd: Skip ground state calculation
        """
        self.taskname = taskname
        self.nempty = nempty
        self.ngridq = ngridq
        self.skipgnd = skipgnd

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


def set_gw_input_string(gs_input:str, gw_input: GWInput):
    """

    Given a converged ground state input, set it
    to repeat the ground state caculation from file (due
    to the additions to the basis) and add the GW inputs

    Note, both replace and format are not inplace

    :return: GW calculation input string
    """

    gs_input = gs_input.replace('do="skip"', 'do="fromfile"')
    gw_input = gw_string_template.format(**gw_input.dict_for_format())

    return gs_input.format(GW_INPUT=gw_input)

