from sklearn.base import BaseEstimator, TransformerMixin
import numpy as np

# Scikit-Learn relies on duck typing (not inheritance).
# All you need to do is create a class and implement three methods:
# fit() (returning self), transform(), and fit_transform().

class CombinedAttributesAdder(BaseEstimator, TransformerMixin):
    def __init__(self, add_bedrooms_per_rooms=True):
        self.add_bedrooms_per_room = add_bedrooms_per_rooms

    # Looks like a dummy routine
    def fit(self, X, y=None):
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

