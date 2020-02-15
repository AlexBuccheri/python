# Common crystals, with atomic basis positions for primitive and conventional
# cell types, in fractional coordinates

from .atoms import Atom

# Silicon unit cell, with basis atoms defined in fractional coordinates
def silicon(basis_type='primitive'):
    assert(basis_type == 'primitive' or basis_type =='conventional')
    # Fractional
    primitive_basis = [Atom('Si',[0,0,0]), Atom('Si',[0.25,0.25,0.25])]
    conventional_basis = [Atom('Si',[0,  0,  0  ]),
                          Atom('Si',[0,  0.5,0.5]),
                          Atom('Si',[0.5,0., 0.5]),
                          Atom('Si',[0.5,0.5,0  ]),
                          Atom('Si',[0.25, 0.25, 0.25]),
                          Atom('Si',[0.25, 0.75, 0.75]),
                          Atom('Si',[0.75, 0.25, 0.75]),
                          Atom('Si',[0.75, 0.75, 0.25])]
    if basis_type == 'primitive':
        return primitive_basis
    elif basis_type == 'conventional':
        return conventional_basis


def magnesium_oxide(basis_type='primitive'):
    assert(basis_type == 'primitive' or basis_type =='conventional')
    quit('Magnesium oxide basis not finished.')
    # Fractional
    #TODO(Alex) Add basis atoms for both cell types
    primitive_basis = [Atom('Mg',[0,0,0]), Atom('O',[0.,0.,0.])]
    conventional_basis = []
    if basis_type == 'primitive':
        return primitive_basis
    elif basis_type == 'conventional':
        return conventional_basis
