import numpy as np
import scipy
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder


def missing_entries_with_pandas(housing, method=3) -> tuple:
    """
    Most Machine Learning algorithms cannot work with missing features
    hence deal with them
    Should really pass attribute string too

    :param housing: housing training data
    :param method: method used to replace missing (NaN) data
    :return: housing data with NaNs replaced
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



class CombinedAttributesAdder(BaseEstimator, TransformerMixin):
    """
    Scikit-Learn relies on duck typing (not inheritance).
    All you need to do is create a class and implement three methods:
    fit() (returning self), transform(), and fit_transform().
    """
    def __init__(self, add_bedrooms_per_rooms=True):
        self.add_bedrooms_per_room = add_bedrooms_per_rooms

    def fit(self, X, y=None):
        """ Dummy routine that is required by skilearn"""
        return self

    def transform(self, X, indices=None) -> np.ndarray:
        """
        Trivially manipulate the housing data to create more informative attributes
        TODO(Alex) Weird that the original attributes (columns) are also retained
        in the returned array
        """
        if indices is None:
            rooms_ix, bedrooms_ix, population_ix, households_ix = 3, 4, 5, 6
        else:
            assert len(indices) == 4, "Expect 4 indices "
            rooms_ix, bedrooms_ix, population_ix, households_ix = indices

        rooms_per_household = X[:, rooms_ix] / X[:, households_ix]
        population_per_household = X[:, population_ix] / X[:, households_ix]
        data = np.c_[X, rooms_per_household, population_per_household]

        if self.add_bedrooms_per_room:
            bedrooms_per_rooms = X[:, bedrooms_ix] / X[:, rooms_ix]
            data = np.c_[data, bedrooms_per_rooms]

        return data

