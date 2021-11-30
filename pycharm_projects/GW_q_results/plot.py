"""
Plot GW results
"""
import matplotlib.pyplot as plt

ha_to_mev = 27211.396132

def plot_gamma_gamma(band_gap_type='gw'):
    assert band_gap_type in ['gw', 'ks']

    # T -> T
    gamma_to_gamma = {'ks': {'222': 0.14204,
                             '444': 0.14204
                             },
                      'gw': {'222': 0.21533,
                             '444': 0.21301
                             }
                      }


    y = [y * ha_to_mev for y in gamma_to_gamma[band_gap_type].values()]
    print('QP gap for Gamma - Gamma (meV):', y)
    plt.xlabel('q-points')
    plt.ylabel('Quasiparticle Gap (meV)')
    plt.plot([8, 64], y, 'bo', markersize=15)
    plt.show()


def plot_x_x(band_gap_type='gw'):
    assert band_gap_type in ['gw', 'ks']

    # X -> X
    x_to_x = {'ks': {'222': 0.13774,
                     '444': 0.13774
                     },
              'gw': {'222': 0.20071,
                     '444': 0.2015
                     }
              }

    y = [y * ha_to_mev for y in x_to_x[band_gap_type].values()]
    print('QP gap for X - X (meV):', y)

    fig, ax = plt.subplots()
    ax.set_xlabel('q-points')
    ax.set_ylabel('Quasiparticle Gap (meV)')
    ax.plot([8, 64], y, color='blue', marker='o', markersize=6)
    plt.show()


def plot_indirect(band_gap_type='gw'):
    """
    where X is VBT and Gamma is CBB
    :param band_gap_type:
    :return:
    """
    assert band_gap_type in ['gw', 'ks']

    # T -> X
    gamma_to_x = {'ks': {'222': 0.12202 ,
                         '444': 0.12202},
                  'gw': {'222': 0.19395,
                         '444': 0.19175}
                  }

    y = [y * ha_to_mev for y in gamma_to_x[band_gap_type].values()]
    print('QP indirect gap for X (VBT) to Gamma (CBB) (in meV):', y)

    fig, ax = plt.subplots()
    ax.set_xlabel('q-points')
    ax.set_ylabel('Quasiparticle Gap (meV)')
    ax.plot([8, 64], y, color='blue', marker='o', markersize=6)
    plt.show()


plot_gamma_gamma()
plot_x_x()
plot_indirect()
