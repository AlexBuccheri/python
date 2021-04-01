"""
Plot classes that 'should' reduce boilerplate code.
If it doesn't reduce code, it's a waste of time.
"""
import matplotlib.pyplot as plt


class Plot:
    """
    Plotting wrapper to allow the overlay of several sets of data on the same plot.
    TODO Be able to pass (line, marker) - [type, size, colour]
    """

    def __init__(self, x, y, x_label='', y_label='', xticklabels=None, yticklabels=None, legend_loc='upper right', legend_title=None):
        self.x = x
        self.y = y
        self.x_label = x_label
        self.y_label = y_label
        self.xticklabels = xticklabels
        self.yticklabels = yticklabels
        self.legend_loc = legend_loc
        self.legend_title = legend_title
        self.initialise_plot()

    def initialise_plot(self):
        """
        Initialise plot axes and labels
        :return:
        """
        self.fig , self.ax = plt.subplots()
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel(self.y_label)

        if self.xticklabels is not None:
            self.ax.set_xticks(self.x)
            self.ax.set_xticklabels(self.xticklabels)

        if self.yticklabels is not None:
            self.ax.set_yticks(self.y)
            self.ax.set_yticklabels(self.yticklabels)

        return

    def plot_data(self, x, y, label=None):
        """
        Basic plot wrapper
        :param x: x data
        :param y: y data
        :param label: optional label for legend
        """
        self.legend_label = label

        if label is None:
            self.ax.plot(x, y)
        else:
            self.ax.plot(x, y, label=str(label))
        return

    def show(self, file_name=None):
        if self.legend_label is not None:
            self.ax.legend(loc=self.legend_loc, title=self.legend_title)

        if file_name is not None:
            plt.savefig(file_name, dpi=300)
        else:
            plt.show()

        return
