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


def inspecting_data(strat_train_set):
    """
    Look at correlations in the data
    Used to inform how one could combine attributes

    :param strat_train_set:
    :return:
    """
    # ----------------------------------------------
    # Inspecting at the data for linear correlation
    # ----------------------------------------------
    housing = strat_train_set.copy()

    # 9x9 matrix
    corr_matrix = housing.corr()

    # The correlation coefficient ranges from –1 to 1
    # The correlation coefficient ONLY measures linear correlations
    print(corr_matrix["median_house_value"].sort_values(ascending=False))

    view_plots = False
    if view_plots:
        view_colour_map(housing)
        view_scatter_plots(housing)

    # Conclusions
    # From view_scatter_plots, the only correlated attribute with target attribute is median_income
    # This could be due to how the data is prepared - see next section.
    #
    # Need to remove patterns from the data that could affect the learning
    # i.e. the cap at $500,000 and the other straight lines visible in fig 2.16 on page 61

    # -------------------------------------------
    # Attribute Combinations
    # Prepare the data in a more sensible way
    # pg 61 - 62
    # -------------------------------------------
    housing["rooms_per_household"] = housing["total_rooms"] / housing["households"]
    housing["bedrooms_per_room"] = housing["total_bedrooms"] / housing["total_rooms"]
    housing["population_per_household"] = housing["population"] / housing["households"]
    corr_matrix = housing.corr()
    print(corr_matrix["median_house_value"].sort_values(ascending=False))

    # Conclusions
    # 'rooms_per_household' and 'bedrooms_per_room' now much more strongly
    # correlate and anti-correlation, respectively, with median_house_value
