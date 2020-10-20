import numpy as np

# Doesn't work because shuffled indices will differ each time one runs the code
# => Could fix the random seed.
# If the random seed is set np.ran dom.seed(42), this will be valid until one fetches an updated dataset.
# I.e. if the reordering is random and the data set changes in size, then one might include a datum in the
# test set that was previously in the training set.

def bad_way_to_split_data_as_train_and_test(data, test_ratio=0.2) -> tuple:
    np.random.seed(42)
    shuffled_indices = np.random.permutation(len(data))
    test_set_size = int(test_ratio * len(data))
    test_indices = shuffled_indices[:test_set_size]
    train_indices = shuffled_indices[test_set_size:]
    return data.iloc[train_indices], data.iloc[test_indices]

# TODO(Alex)
# What is the identifier
# How is this creating a HASH value?
# How is this any different to just using the order in which the data is in the file ? (Maybe it
# it gives consistent results if the ordering changes)
# def test_set_check(identifier, test_ratio):
#     return crc32(np.int64(identifier)) & 0xffffffff < test_ratio * 2**32
#
# Can't use this as the housing dataset does not have an identifier column
# def split_train_test_by_id(data, test_ratio, id_column):
#     ids = data[id_column]
#     in_test_set = ids.apply(lambda id_: test_set_check(id_, test_ratio))
#     return data.loc[~in_test_set], data.loc[in_test_set]

# See page 53
# The simplest solution is to use the row index as the ID... which is what I was rationalising above
# This breaks if the order of the data changes
# If this is not possible, then you can try to use the most stable features to build a unique identifier
# i.e. longitude and latitude in the housing dataset
#
# TODO(Alex) Should write these functions above in a separate file and test

# Can use ski-kit learn's function
# Can take multiple DataFrames and split each w.r.t the same indices
# Using one DataFrame here
#train_set, test_set = train_test_split(housing_data, test_size=0.2, random_state=42)