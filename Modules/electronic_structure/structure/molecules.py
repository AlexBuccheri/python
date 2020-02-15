import numpy as np

from atoms import Atom
from parameters import molecules as param

# Should have something selecting bohr or ang + degrees or radian
def water(r1=None, r2=None, theta_radian=None):

    #O-H bond length
    _r_bohr = param.water['o-h']['bohr']
    _r_angstrom = param.water['o-h']['angstrom']
    _theta_degrees = param.water['theta']['degrees']
    _theta_radian = _theta_degrees * (np.pi/180.)

    def _water(r1 , r2 , theta):
        molecule = [Atom('O', [0, 0, 0]),
                    Atom('H', [0,  r1 * np.sin(theta / 2), r1 * np.cos(theta / 2)]),
                    Atom('H', [0, -r2 * np.sin(theta / 2), r2 * np.cos(theta / 2)])]
        return molecule

    molecule = []
    if (r1 and r2 and theta_radian) == None:
        molecule = _water(_r_angstrom, _r_angstrom, _theta_radian)
    else:
        molecule = _water(r1, r2, theta_radian)

    return molecule


def carbon_dioxide(bond_length=None):
    if bond_length == None:
        bl = param.carbon_dioxide['c-h']['angstrom']
    else:
        bl = bond_length
    return [Atom('C', [0,0,0]), Atom('O', [0,0,-bl]), Atom('O',[0,0,bl])]


def acetylene(bond_length_ch=None, bond_length_cc=None):

    if bond_length_ch == None:
        bl_ch = param.acetylene['c-h']['angstrom']
    else:
        bl_ch = bond_length_ch

    if bond_length_cc == None:
        bl_cc = param.acetylene['c-c']['angstrom']
    else:
        bl_cc = bond_length_ch

    # Centred on [0,0,0]
    return [Atom('H',[0,0,-bl_ch-0.5*bl_cc]), Atom('C', [0,0,-0.5*bl_cc]),
                Atom('C',[0,0,0.0.5*bl_cc]), Atom('H',[0,0,bl_ch+0.5*bl_cc])]
