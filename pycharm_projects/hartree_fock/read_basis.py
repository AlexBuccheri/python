import json
import typing

an_to_label = {1:'H', 2:'He'}

def sto_3g(fname: str) -> tuple:

    fid = open(fname, "r")
    data = json.load(fid)
    fid.close()

    print(data['elements']['1']['electron_shells'][0]['exponents'])
    for element in data['elements']:
        an = int(element)
        print(element['electron_shells'])
        # ang_mom = element['electron_shells']['angular_momentum']
        # exponents = element['electron_shells']['exponents']
        # coefficients = element['electron_shells']['coefficients']

    #print(data)
    return