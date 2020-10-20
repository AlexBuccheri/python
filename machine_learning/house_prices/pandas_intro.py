"""
10 minutes to Pandas
https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html
"""

import pandas as pd
import numpy as np

# Series. Pandas creates a default integer index
s = pd.Series([1, 3, 5, np.nan, 6, 8])
print(s)

# DataFrame for dates
dates = pd.date_range('20130101', periods=6)
#  data, index = row labels, columns = column labels
df = pd.DataFrame(np.random.randn(6, 4), index=dates, columns=list('ABCD'))
print(df)

# Initiliase DataFrame with a dictionary.
# keys = columns, values can be a whole range of things
df2 = pd.DataFrame({'A': 1.,
                    'B': pd.Timestamp('20130102'),
                    'C': pd.Series(1, index=list(range(4)), dtype='float32'),
                    'D': np.array([3] * 4, dtype='int32'),
                    'E': pd.Categorical(["test", "train", "test", "train"]),
                    'F': 'foo'})
print(df2)
print(df2.dtypes)

# Access column dtype
print(df2.A)

# Viewing data
# df.head()
# df.tail(3)

# Row and column labels
#print(df.index)
#print(df.columns)

# NumPy arrays have one dtype for the entire array, while pandas DataFrames have one dtype per column.
# When you call DataFrame.to_numpy(), pandas will find the NumPy dtype that can hold all of the dtypes
# in the DataFrame. This may end up being object, which requires casting every value to a Python object.

# Convert DataFrame to numpy
print(df.to_numpy())

print(df.at[dates[0], 'A'])
print(df.loc['20130101', 'A'])
#print(df['20130101', 'A'])
#print()


