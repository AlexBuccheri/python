import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import pandas as pd
import scipy

from get_data import HOUSING_PATH
import load
import stratum
import visualise
import transformers


def print_info(housing):
    # Data dimensions and labels
    print(housing.info())

    # Preview pandas Dataframe object
    print(housing.head())

    # Get info on categories (number of entry types) for object "ocean_proximity"
    print(housing["ocean_proximity"].value_counts())

    # Generates descriptive statistics that summarize the central tendency,
    # dispersion and shape of a dataset's distribution, excluding `NaN` values.
    print(housing.describe())


def plot_info(housing):
    """ Plot histrogram for each of the initial 9 attributes (columns) in housing"""
    housing.hist(bins=50, figsize=(20, 15))
    plt.show()


def split_data(housing_data) -> tuple:
    """
    Splitting into training and test sets
    :return:
    """
    housing_data = stratum.stratify_median_income_category(housing_data, show_plot=False)
    strat_train_set, strat_test_set = stratum.split_training_test(housing_data)

    # Remove the income_cat attribute so the data is back to its original state
    # TODO(Alex) Iterating over a tuple?
    for set_ in (strat_train_set, strat_test_set):
        set_.drop("income_cat", axis=1, inplace=True)

    return strat_train_set, strat_test_set


def inspecting_data(strat_train_set):
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
        visualise.view_colour_map(housing)
        visualise.view_scatter_plots(housing)

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


def missing_entries_with_pandas(housing, method=3) -> tuple:
    """
    Most Machine Learning algorithms cannot work with missing features
    hence deal with them
    Should really pass attribute string too

    :param housing:
    :param method:
    :return:
    """
    assert type(method) == int

    if method == 1:
        # 1. Get rid of corresponding districts (entries with missing features)
        # drops entries that have N/A for their value
        housing.dropna(subset=["total_bedrooms"])
        return (housing)
    elif method == 2:
        # 2. Get rid of the whole attribute (total_bedrooms)
        # Drops whole column
        housing.drop("total_bedrooms", axis=1)
        return (housing)
    elif method == 3:
        # 3. Set missing values to some sensible value (i.e. 0, mean, median, etc)
        median = housing["total_bedrooms"].median()
        housing["total_bedrooms"].fillna(median, inplace=True)
        return (housing, median)
    else:
        quit("method choice is invalid:", method)


def fill_missing_entries(data: pd.DataFrame, strategy='median') -> tuple :
    """
    Replace each attribute’s missing values with the median of that attribute,
    using sklearn.

    NOTE: Performed per column. There might be more efficient pandas operations to do this.

    :param data: Input pd.DataFrame. Missing data MUST only be numerical
    :param strategy: Strategy for replacing missing (NaN) values in numerical data
    :return: pd.DataFrame data with NaN values replaced with fill values
             and the corresponding fill values
    """

    attributes_with_null_entries = []
    for key in data.keys():
        if data[key].isnull().any():
            attributes_with_null_entries.append(key)
    data_with_null_entries = data[attributes_with_null_entries].copy()

    for data_type in data_with_null_entries.dtypes.tolist():
        assert data_type in [np.int32, np.int64, np.float32, np.float64], \
            "data type of column with NaN entry/entries is not numerical"

    # Replaces NaN values with attribute median in 'data_with_null_entries'
    imputer = SimpleImputer(strategy=strategy)
    imputer.fit(data_with_null_entries)

    # Impute all values into X
    X = imputer.transform(data_with_null_entries)
    assert type(X) == np.ndarray
    assert X.shape == (16512, 1)

    # Replace columns (attributes) in data that contain ANY NaN values
    for i, key in enumerate(attributes_with_null_entries):
        data.loc[:, key] = X[:, i]

    return data, imputer.statistics_


def text_to_numerical(data: pd.Series):
    """
    Convert text attributes to numerical representations for a given
    pd.Series.

    one-hot encoding: Create one binary attribute per category.
    Represented in a sparse matrix of binary entries.
    Used when the text categories DO NOT map onto a linear numerical
    scale.

    View categories with onehot_encoder.categories_

    :param data: pd.Series containing a column of text-based categories (entries)
    :return: Numerical representations for the categories in scipy.sparse.csr.csr_matrix
    """

    # Not this function's responsibility to screen attributes
    # numerical_data_types = [np.int32, np.int64, np.float32, np.float64]
    # for key in data.keys():
    #     if data[key].dtype not in numerical_data_types:

    onehot_encoder = OneHotEncoder()
    data_encoded = onehot_encoder.fit_transform(data)
    assert type(data_encoded) == scipy.sparse.csr.csr_matrix
    return data_encoded



# --------------------------------------
# Main Routine
# --------------------------------------

housing_data = load.load_housing_data(HOUSING_PATH)
strat_train_set, strat_test_set = split_data(housing_data)
#inspecting_data(strat_train_set)

# Data cleaning
# Separate the predictors and the labels
# where labels are the target values
# Create a copy, with label name(s) and corresponding axis removed
housing = strat_train_set.drop("median_house_value", axis=1)
housing_labels = strat_train_set["median_house_value"].copy()

# ---------------------
# Fill missing values
# ---------------------
# Rather than follow the strategy in the book of dropping non-numerical attributes (columns)
# and filling all NaN values, I've operated only on columns with NaN values then replaced
# these columns in the full DataFrame. ref pg 63
housing, median_fill = fill_missing_entries(housing, strategy='median')

# We know that "ocean_proximity" is the only non-numerical attribute
# Because the data returned is 2D sparse array, this clearly cannot be put back into
# the housing DataFrame
op_encoded = text_to_numerical(housing[["ocean_proximity"]])

# ---------------------
# Custom Transformers
# ---------------------
# Remove the text-based attribute
# housing.drop(columns="ocean_proximity")

# add_bedrooms_per_rooms = hyper-parameter
# This hyper-parameter will allow you to easily find out whether adding this
# attribute helps the Machine Learning algorithms or not.
#
attr_adder = transformers.CombinedAttributesAdder(add_bedrooms_per_rooms=False)
housing_extra_attribs = attr_adder.transform(housing.to_numpy())

# ------------------------------------------------------------------------
# Create a pipeline
# We have a preprocessing pipeline that takes the full housing data
# and applies the appropriate transformations to each column.
# ------------------------------------------------------------------------

num_pipeline = Pipeline([
    ('imputer', SimpleImputer(strategy="median")),
    ('attribs_adder', transformers.CombinedAttributesAdder()),
    ('std_scaler', StandardScaler())
])

numerical_columns = list(housing)
numerical_columns.remove('ocean_proximity')

full_pipeline = ColumnTransformer([
    # Name,   transformer or pipeline of transformers, list of column labels or indices to apply the transformer to
    ("numerical", num_pipeline, numerical_columns),
    ("categories", OneHotEncoder(), ["ocean_proximity"]),
])

housing_prepared = full_pipeline.fit_transform(housing)

print(housing_prepared)


# In the text: housing_num_tr = num_pipeline.fit_transform(housing_num)
# In my code: housing, after removing "ocean_proximity"

# NOTE: ColumnTransformer can handle both categorical columns and the numerical column
# See page 71

# Train the model




# -----------------------------------
# Feature Scaling of the training data
# -----------------------------------
# Machine Learning algorithms don’t perform well when the input numerical
# attributes have very different scales.

# Approaches
# 1. Normalise the data (min-max scaling) to [0:1]
# Always required for neural networks
# Use MinMaxScaler in sklearn

# 2. Standardisation scaling
# Doesn't bound to range [0:1] but is much less affected by outliers
# Use StandardScaler in sklearn transformers




# For the next chapter, write the front end with a notebook in pycharm:
# https://www.jetbrains.com/help/pycharm/jupyter-notebook-support.html#ui







