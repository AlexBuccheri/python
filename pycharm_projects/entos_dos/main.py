"""
Example script demonstrating plotting the density of states output
from ENTOS Qcore
"""
import matplotlib.pyplot as plt

from read_dos import get_dos
from plot_dos import plot_total_density_of_states, PlotSettings

dos = get_dos('example/bulk_carbon.json')
fig, ax = plot_total_density_of_states(dos, PlotSettings())
plt.show()
