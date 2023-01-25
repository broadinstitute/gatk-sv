from typing import Iterable
import numpy
import pandas


def assert_indices_equal(i1: pandas.Index, i2: pandas.Index, context: str):
    assert numpy.array_equal(i1.isnull(), i2.isnull()), f"{context}: index null values not equal"
    assert numpy.array_equal(i1[~i1.isnull()], i2[~i2.isnull()]), \
        f"{context}: index non-null values not equal"


def assert_series_equal(
        s1: pandas.Series,
        s2: pandas.Series,
        context: str,
        check_index: bool = True,
        dtypes_equal: bool = True,
        rtol: float | None = None,
        atol: float | None = None
):
    """
    Assert that two series are equal.
    If
    """
    if check_index:
        # noinspection PyTypeChecker
        assert_indices_equal(s1.index, s2.index, context=context)
    if dtypes_equal:
        if pandas.api.types.is_categorical_dtype(s1) and pandas.api.types.is_categorical_dtype(s2):
            # noinspection PyUnresolvedReferences
            assert s1.dtype == s2.dtype, (
                f"{context}: categorical dtypes not equal: {s1.dtype.categories} != "
                f"{s2.dtype.categories}"
            )
        else:
            # For reasons I can't understand, checking dtype equality sometimes fails for String
            # dtypes, but checking for dtype.__class__ equality works
            assert s1.dtype.__class__ == s2.dtype.__class__, (
                f"{context}: dtypes not equal: {s1.dtype.__class__.__name__} != "
                f"{s2.dtype.__class__.__name__}"
            )

    if pandas.api.types.is_float_dtype(s1) and pandas.api.types.is_float_dtype(s2):
        tol_overrides = {}
        if rtol is not None:
            tol_overrides["rtol"] = rtol
        if atol is not None:
            tol_overrides["atol"] = atol
        assert numpy.allclose(s1, s2, equal_nan=True, **tol_overrides), \
            f"{context}: values not equal"
    else:
        # could be an object or an extension dtype, so could have nulls
        assert numpy.array_equal(s1.isnull(), s2.isnull()), f"{context}: null values not equal"
        assert numpy.array_equal(s1[~s1.isnull()], s2[~s2.isnull()]), \
            f"{context}: non-null values not equal"


def assert_sets_equal(
        values1: Iterable,
        values2: Iterable,
        context: str
):
    set1 = set(values1)
    set2 = set(values2)
    if set1 == set2:
        return  # no problem
    set1_not_2 = ','.join(sorted(str(obj) for obj in set1.difference(set2)))
    set2_not_1 = ','.join(sorted(str(obj) for obj in set2.difference(set1)))
    if set1_not_2:
        if set2_not_1:
            delta_str = f"in 1st not second: {set1_not_2}\nin 2nd not first: {set2_not_1}"
        else:
            delta_str = f"in 1st not second: {set1_not_2}"
    else:
        # must be set2_not_1
        delta_str = f"in 2nd not first: {set2_not_1}"
    assert set1 == set2, f"{context}\n{delta_str}"


def assert_dataframes_equal(
        df1: pandas.DataFrame,
        df2: pandas.DataFrame,
        context: str,
        check_index: bool = True,
        check_column_order: bool = True,
        dtypes_equal: bool = True,
        rtol: float | None = None,
        atol: float | None = None
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
            assert_sets_equal(df1.columns, df2.columns, f"{context}: extra/missing columns")
            assert df1.columns.equals(df2.columns), f"{context}: columns orders differ"
    else:
        assert_sets_equal(df1.columns, df2.columns, f"{context}: extra/missing columns")
    if check_index:
        assert_indices_equal(df1.index, df2.index, context)
    for column in df1.columns:
        assert_series_equal(
            df1.loc[:, column], df2.loc[:, column],
            context=f"{context}: {column}", check_index=False, dtypes_equal=dtypes_equal,
            rtol=rtol, atol=atol
        )
