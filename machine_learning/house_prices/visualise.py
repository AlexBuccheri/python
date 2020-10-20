import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import scatter_matrix


def view_colour_map(housing: pd.DataFrame):

    # s = size of scatter points. Determined from population
    # c = colour map. Determined by the label 'median_house_price' in 'housing'

    housing.plot(kind="scatter", x="longitude", y="latitude", alpha=0.4, s=housing["population"] / 100,
                 label="population", figsize=(10, 7), c="median_house_value", cmap=plt.get_cmap("jet"), colorbar=True)
    plt.legend()
    plt.show()
    return


def view_scatter_plots(housing: pd.DataFrame):
    """
    Look at correlation between a subset of attributes
    (Don't want to do so for all attributes as 11 x 11 = 121 plots)
    """

    # scatter_matrix() plots every numerical attribute against every other numerical attribute
    # Use a subset of attributes

    # Target attribute: median_house_value
    # From these plots, clearly the only correlated attribute with target attribute is median_income

    attributes = ["median_house_value", "median_income", "total_rooms", "housing_median_age"]
    scatter_matrix(housing[attributes], figsize=(12, 8))
    plt.show()
    housing.plot(kind="scatter", x="median_income", y="median_house_value", alpha=0.1)
    plt.show()
    return
