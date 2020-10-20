import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt

# Pipelines
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer

# Transformers
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder

# My modules for this task
from get_data import HOUSING_PATH
import load
import stratum
import visualise
import transformers


# My functions


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


def transformer_fill_missing_values(housing):
    """
    Example data transformers that fills missing values for
    a) Numerical values
    b) Categorical values (which requires assigning some numerial value first)

    """
    # Numerical Data
    # Rather than follow the strategy in the book of dropping non-numerical attributes (columns)
    # and filling all NaN values, I've operated only on columns with NaN values then replaced
    # these columns in the full DataFrame. ref pg 63
    housing, median_fill = transformers.fill_missing_entries(housing, strategy='median')

    # Categorical Data
    # We know that "ocean_proximity" is the only non-numerical attribute
    # Because the data returned is 2D sparse array, this clearly cannot be put back into
    # the housing DataFrame
    op_encoded = transformers.text_to_numerical(housing[["ocean_proximity"]])

    return housing, median_fill, op_encoded


def transformer_combine_attributes(housing):
    """
    Custom transformer for combining attributes of the housing data

    :param housing: Housing training data
    :return: housing data with additional, combined attributes
    ALTHOUGH, our custom class is returning it in a numpy array rather than
    a pandas DataFrame
    """
    # Remove the text-based attribute
    # housing.drop(columns="ocean_proximity")

    # add_bedrooms_per_rooms = hyper-parameter
    # This hyper-parameter will allow you to easily find out whether adding this
    # attribute helps the Machine Learning algorithms or not.

    attr_adder = transformers.CombinedAttributesAdder(add_bedrooms_per_rooms=False)
    housing_extra_attribs = attr_adder.transform(housing.to_numpy())

    return housing_extra_attribs


# -----------------------------------
# Feature Scaling of the training data
# -----------------------------------
# Machine Learning algorithms don’t perform well when the input numerical
# attributes have very different scales.
#
# Approaches
# 1. Normalise the data (min-max scaling) to [0:1]
# Always required for neural networks
# Use MinMaxScaler in sklearn
#
# 2. Standardisation scaling
# Doesn't bound to range [0:1] but is much less affected by outliers
# Use StandardScaler in sklearn transformers


# ----------------------------------------------------------------------------
# Main Routine
# Load data, split into training and test sets (ensuring no biasing in
# how the test set is sampled (stratum), and set up a pipeline containing
# transformers that operate on the training data
# ----------------------------------------------------------------------------

housing_data = load.load_housing_data(HOUSING_PATH)
strat_train_set, strat_test_set = split_data(housing_data)
# inspecting_data(strat_train_set)

# Separate the predictors and the labels (target values)
# Create a copy, with label name(s) and corresponding axis removed
housing = strat_train_set.drop("median_house_value", axis=1)
housing_labels = strat_train_set["median_house_value"].copy()

# Create a pipeline
# i.e. a set of transformers that are sequentially applied to the data
#
# We have a pre-processing pipeline that takes the full housing data
# (numerical and categorised attributes)
# and applies the appropriate transformations to each column.

# TODO(Alex) Try combining imputer, attribs_adder and std_scaler
# into my own class, that contains the methods:
# fit() (returning self), transform(), and fit_transform().

num_pipeline = Pipeline([
    ('imputer', SimpleImputer(strategy="median")),
    ('attribs_adder', transformers.CombinedAttributesAdder()),
    ('std_scaler', StandardScaler())
])

numerical_columns = list(housing)
numerical_columns.remove('ocean_proximity')

# Would be useful to have a function that splits numerical and categorised data (I have that)
full_pipeline = ColumnTransformer([
    # Name,   transformer or pipeline of transformers, list of column labels or indices to apply the transformer to
    ("numerical", num_pipeline, numerical_columns),
    ("categories", OneHotEncoder(), ["ocean_proximity"]),
])

housing_prepared = full_pipeline.fit_transform(housing)

print(housing_prepared)






# Train the model









# For the next chapter, write the front end with a notebook in pycharm:
# https://www.jetbrains.com/help/pycharm/jupyter-notebook-support.html#ui







