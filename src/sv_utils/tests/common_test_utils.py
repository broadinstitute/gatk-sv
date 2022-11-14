import numpy
import pandas


def assert_indices_equal(i1: pandas.Index, i2: pandas.Index, context: str):
    assert numpy.array_equal(i1.isnull(), i2.isnull()), f"{context}: index null values not equal"
    assert numpy.array_equal(i1[~i1.isnull()], i2[~i2.isnull()]), f"{context}: index non-null values not equal"


def assert_series_equal(s1: pandas.Series, s2: pandas.Series, context: str, check_index: bool = True):
    if check_index:
        # noinspection PyTypeChecker
        assert_indices_equal(s1.index, s2.index, context=context)

    assert numpy.array_equal(s1.isnull().values, s2.isnull().values), f"{context}: null values not equal"
    assert numpy.array_equal(s1.loc[~s1.isnull().values], s2.loc[~s2.isnull().values]),\
        f"{context}: non-null values not equal"


def assert_dataframes_equal(
        df1: pandas.DataFrame,
        df2: pandas.DataFrame,
        context: str,
        check_index: bool = True,
        check_column_order: bool = True
):
    """
    For the purposes of testing we *DO* care about:
        - values in each column
        - names of each column
        - row order
    we *MAY* care about:
        - the index
        - column order
    we *DON'T* care about:
        - if columns are categorical are not (otherwise just use pandas.DataFrame.equals)
    Args:
        df1: pandas.DataFrame
            first DataFrame to test for equality
        df2:  pandas.DataFrame
            second DataFrame to test for equality
        context: str
            extra message to put in assertion if it fails
        check_index: bool
            if True, assert the indices are identical
        check_column_order: bool
            if True, assert the columns are in the same order
    """
    if check_column_order:
        if not df1.columns.equals(df2.columns):
            assert not set(df1.columns).symmetric_difference(df2.columns), f"{context}: extra/missing column values"
            assert df1.columns.equals(df2.columns), f"{context}: columns orders differ"
    else:
        assert not set(df1.columns).symmetric_difference(df2.columns), f"{context}: column values not equal"
    if check_index:
        assert_indices_equal(df1.index, df2.index, context)
    for column in df1.columns:
        assert_series_equal(df1[column], df2[column], context=f"{context}: {column}", check_index=False)
