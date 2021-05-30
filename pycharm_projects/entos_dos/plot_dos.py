import matplotlib.pyplot as plt
import numpy as np

from read_dos import resolve_energy_grid

class PlotSettings:
    """
    Labels and settings for density of state plot
    """
    def __init__(self, xlabel='Energy (Ha)', ylabel='Density of States (1/Ha)',
                 xmin=None, xmax=None, ymin=0, ymax=None, line_style='b'):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.line_style = line_style


def plot_total_density_of_states(dos: dict, settings: PlotSettings):
    """
    Plot the density of states
    :param dos: density of states  JSON output from Qcore
    :param settings : Plot settings object
    :return: fig and ax plot objects
    """
    energies = resolve_energy_grid(dos)
    dos_values = dos['total_density_of_states']
    assert len(energies) == len(dos_values)

    xmin = min(energies) if settings.xmin is None else settings.xmin
    xmax = max(energies) if settings.xmax is None else settings.xmax
    ymin = 0 if settings.ymin is None else settings.ymin
    ymax = 1.1 * np.max(dos_values) if settings.ymax is None else settings.ymax

    fig, ax = plt.subplots()
    ax.set_xlabel(settings.xlabel)
    ax.set_ylabel(settings.ylabel)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.plot(energies, dos_values)
    plt.axvline(x=dos['fermi_level'], color='black')

    return fig, ax
