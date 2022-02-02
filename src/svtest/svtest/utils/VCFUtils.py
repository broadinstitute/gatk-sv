#!/usr/bin/env python

"""
Useful utilities for working with pysam variant objects
"""


def get_info_field(record, name):
    if name not in record.info:
        if name == 'SVLEN':
            if record.info['SVTYPE'] in ['DEL', 'DUP', 'INV']:
                record.info['SVLEN'] = record.stop - record.pos
            else:
                record.info['SVLEN'] = -1
        else:
            raise ValueError("%s info field not found: %s" %
                             (name, record.info.keys()))
    return record.info[name]


def get_record_length(record):
    return get_info_field(record, "SVLEN")


def get_sv_type(record, expected_types):
    if "SVTYPE" not in record.info:
        raise ValueError("SVTYPE info field not found: %s" %
                         record.info.keys())
    type = record.info["SVTYPE"]
    if type not in expected_types:
        raise ValueError("Unexpected SVTYPE: %s" % type)
    return type


def get_evidence_types(record, expected_types):
    if "EVIDENCE" not in record.info:
        raise ValueError("EVIDENCE info field not found: %s" %
                         record.info.keys())
    types = record.info["EVIDENCE"]
    for type in types:
        if type not in expected_types:
            raise ValueError("Unexpected EVIDENCE: %s" % type)
    return types
