#!/usr/bin/env python

"""
Useful utilities for IO.
"""


def read_samples_list(path, fail_on_empty=True):
    with open(path, 'r') as f:
        samples = f.read().splitlines()
    if fail_on_empty and len(samples) == 0:
        raise ValueError("Samples list empty")
    return samples


def read_contig_list(path, fail_on_empty=True):
    with open(path, 'r') as f:
        contigs = f.read().splitlines()
    if fail_on_empty and len(contigs) == 0:
        raise ValueError("Contig list empty")
    return contigs
