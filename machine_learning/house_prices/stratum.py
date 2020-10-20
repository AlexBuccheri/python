"""
# stratified sampling

* Strata = homogeneous subgroup

* The right number of instances are sampled from each stratum to guarantee
  that the test set is representative of the overall population.

* In our example, median income is an important attribute for predicting median house price

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedShuffleSplit


def stratify_median_income_category(housing: pd.DataFrame, show_plot=True):
    """
    Discretise the median incomes into bins/strata
    ref pg 54

    Instead of an array of indices [0:20640] and a corresponding median income
    (as given by housing["median_income"]), housing["income_cat"] returns the
    indices and a bin/stratum associated with the median income.
    i.e.  index i, median income 3.2
    ->    index i, stratum 3
    based on the bin ranges below
    """
    # create an income category attribute with five categories
    # where bin 0 contains:  0 <= median_income < 1.5
    # (1.5 == $15000)
    assert type(housing) == pd.DataFrame, "housing must be type(pd.DataFrame)"

    # TODO(Alex) Makes more sense to go to the max median income, rather than np.inf?
    housing["income_cat"] = pd.cut(housing["median_income"],
                                   bins=[0., 1.5, 3.0, 4.5, 6., np.inf],
                                   labels=[1, 2, 3, 4, 5])

    # TODO(Alex) Have x labels as the bin ranges
    if show_plot:
        housing["income_cat"].hist()
        plt.show()

    return housing


def split_training_test(housing: pd.DataFrame) -> tuple:
    """
    Training and testing data (not indices) applied to the stratised
    data for median income, rather than the raw median income

    Should avoid bias in random sampling of the total set to create
    test data
    """
    # Generate indices to split data into training and test set.
    # StratifiedShuffleSplit returns an object
    split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)

    # .split(Training data, The target variable for supervised learning problems)
    # My guess is this works by counting the number of instances for each housing["income_cat"]
    # (which has been constructed to vary from 1 to 5) and weighted appropriately

    for train_index, test_index in split.split(housing, housing["income_cat"]):
        strat_train_set = housing.loc[train_index]
        strat_test_set = housing.loc[test_index]

    return strat_train_set, strat_test_set


def inspect_training_test(strat_train_set, strat_test_set):
    # Ref pg 55
    # As you can see, the test set generated using stratified sampling has income
    # category proportions almost identical to those in the full dataset,
    # whereas the test set generated using purely random sampling is skewed.
    # TODO(Alex) Count what one gets from purely random sampling
    print(strat_test_set["income_cat"].value_counts() / len(strat_test_set))
    print(strat_train_set["income_cat"].value_counts() / len(strat_train_set))
    return