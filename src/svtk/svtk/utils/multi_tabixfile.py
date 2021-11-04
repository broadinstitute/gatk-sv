#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Access multiple TabixFiles from a single interface.

Fetching from the list of TabixFiles returns an iterator over merge-sorted
entries in the tabix files.
"""

import heapq
import pysam
from svtk.utils import is_smaller_chrom


class _TabixRow:
    def __init__(self, row):
        """
        Sortable wrapper around a row in a tabix file.

        Can be instantiated from a str or from a pysam TupleProxy.

        Parameters
        ----------
        row : str or pysam.TupleProxy

        Attributes
        ----------
        row : str
            Tab-delimited row entry
        tup : tup
            Tuple representation of the row
        base_type : type
            str or tuple
        """

        if isinstance(row, str):
            self.row = row
            self.tup = tuple(row.strip().split('\t'))
            self.base_type = str
        elif isinstance(row, pysam.libctabixproxies.TupleProxy):
            self.row = str(row)
            self.tup = tuple(row)
            self.base_type = tuple
        else:
            t = type(row).__name__
            raise Exception('Invalid type for tabix row: {0}'.format(t))

    def __lt__(self, other):
        if self.tup[0] == other.tup[0]:
            return int(self.tup[1]) < int(other.tup[1])
        else:
            return is_smaller_chrom(self.tup[0], other.tup[0])


class _SortableTabixIterator:
    def __init__(self, iterator):
        """
        Wrapper around the base pysam.TabixIterator that permits merge-sorting
        of the tabix records

        Parameters
        ----------
        iterator : pysam.TabixIterator
        """
        self.iterator = iterator

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        return _TabixRow(next(self.iterator))


class _TabixIterator:
    def __init__(self, iterator):
        """
        Helper class to convert _TabixRow objects back to original parsed type

        Parameters
        ----------
        iterator : _SortableTabixIterator
        """
        self.iterator = iterator

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        try:
            row = next(self.iterator)
        except StopIteration:
            raise StopIteration

        if row.base_type == str:
            return row.row
        else:
            return row.tup


class MultiTabixFile:
    def __init__(self, tabixfiles):
        """
        From multiple tabixfiles, simulate access into a single tabixfile

        Parameters
        ----------
        tabixfiles : list of pysam.TabixFile

        Attributes
        ----------
        tabixfiles : list of pysam.TabixFile
        """

        self.tabixfiles = tabixfiles

    def fetch(self, *args, **kwargs):
        """
        TabixFile.fetch(self, reference=None, start=None, end=None, region=None,
                        parser=None, multiple_iterators=False)
        """
        iterators = []
        for tbx in self.tabixfiles:
            iterators.append(_SortableTabixIterator(
                tbx.fetch(*args, **kwargs)))

        return _TabixIterator(heapq.merge(*iterators))

    def close(self):
        for tbx in self.tabixfiles:
            tbx.close()
