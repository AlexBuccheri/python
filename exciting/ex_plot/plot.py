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

    def __init__(self, x, y, x_label='', y_label='', xticklabels=None, yticklabels=None, legend_loc='upper right',
                 legend_title=None, label_size=14, font_size=14, linewidth=1):
        self.x = x
        self.y = y
        self.x_label = x_label
        self.y_label = y_label
        self.xticklabels = xticklabels
        self.yticklabels = yticklabels
        self.legend_loc = legend_loc
        self.legend_title = legend_title
        self.label_size = label_size
        self.font_size = font_size
        self.line_width = linewidth
        self.initialise_plot()

    def initialise_plot(self):
        """
        Initialise plot axes and labels
        :return:
        """
        self.fig , self.ax = plt.subplots()
        self.ax.set_xlabel(self.x_label, fontsize=self.font_size)
        self.ax.set_ylabel(self.y_label, fontsize=self.font_size)

        if self.xticklabels is not None:
            self.ax.set_xticks(self.x)
            self.ax.set_xticklabels(self.xticklabels)

        if self.yticklabels is not None:
            self.ax.set_yticks(self.y)
            self.ax.set_yticklabels(self.yticklabels)

        self.ax.tick_params(axis='both', which='major', labelsize=self.label_size)

        return

    def plot_data(self, x, y, label=None):
        """
        Basic plot wrapper
        :param x: x data
        :param y: y data
        :param label: optional label for legend
        # TODO Probably a better way to pass options as dictionary,
        # as with save
        """
        self.legend_label = label

        if label is None:
            self.ax.plot(x, y, linewidth = self.line_width)
        else:
            self.ax.plot(x, y, linewidth = self.line_width, label=str(label))
        return

    def save(self, save_options:dict):
        if self.legend_label is not None:
            self.ax.legend(loc=self.legend_loc, title=self.legend_title)

        file_name = save_options.pop('file_name')
        plt.savefig(file_name, **save_options)

    def show(self):
        if self.legend_label is not None:
            self.ax.legend(loc=self.legend_loc, title=self.legend_title)
        plt.show()

