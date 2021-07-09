#!/usr/bin/env python

"""
Useful utilities for testing.
"""


def test_sets_equal(iterable_a, iterable_b, item_str="item", name_a="set A", name_b="set B"):
    set_a = set(iterable_a)
    set_b = set(iterable_b)
    a_not_in_b = set_a - set_b
    if len(a_not_in_b) > 0:
        raise ValueError("One or more %s(s) found in %s but not %s: %s" % (
            item_str, name_a, name_b, a_not_in_b))
    b_not_in_a = set_b - set_a
    if len(b_not_in_a) > 0:
        raise ValueError("One or more %s(s) found in %s but not %s: %s" % (
            item_str, name_b, name_a, b_not_in_a))


def test_iterable_sizes_equal(iterable_a, iterable_b, name_a="iterable A", name_b="iterable B"):
    if len(iterable_a) != len(iterable_b):
        raise ValueError("%s (%d) was not the same size as %s (%d)" % (
            name_a, len(iterable_a), name_b, len(iterable_b)))


def test_iterable_size(iterable, size):
    if len(iterable) != size:
        raise ValueError("Expected %d values but found %d in: %s" %
                         (size, len(iterable), str(iterable)))


def test_column_equals(columns, idx, val):
    if columns[idx] != val:
        raise ValueError("Expected column %d to equal %s but found %s" % (
            idx + 1, val, columns[idx]))


def test_column_in_iterable(columns, idx, iter, msg=None):
    if columns[idx] not in iter:
        if msg is None:
            raise ValueError("Expected column %d to be one of %s but found %s" % (
                idx + 1, str(iter), columns[idx]))
        else:
            raise ValueError(msg)


def test_is_not_empty(list, name):
    if len(list) == 0:
        raise ValueError("List %s was empty" % name)


def test_is_float(columns, idx):
    if not is_float(columns[idx]):
        raise ValueError("Column %d was not a float: %s" %
                         (idx + 1, str(columns)))


def test_is_int(columns, idx):
    if not is_int(columns[idx]):
        raise ValueError("Column %d was not an int: %s" %
                         (idx + 1, str(columns)))


def is_float(str):
    try:
        float(str)
        return True
    except ValueError:
        return False


def is_int(str):
    try:
        int(str)
        return True
    except ValueError:
        return False
