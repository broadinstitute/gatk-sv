#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from functools import reduce
import math

# Caches genotype carrier status for _is_non_ref()
_carrier_map = {}


def choose_best_genotype(sample, records):
    """
    Return record where sample has best non-reference genotype, if available

    Parameters
    ----------
    sample : str
    records : list of pysam.VariantRecord

    Returns
    -------
    best_record : pysam.VariantRecord
    """

    # Returns true if the GT is called and non-ref
    # TODO: does not handle multi-allelic CNVs
    def _is_non_ref(format_fields):
        gt = format_fields.get('GT', None)
        if gt in _carrier_map:
            return _carrier_map[gt]
        else:
            is_carrier = gt is not None and any(a is not None and a > 0 for a in gt)
            _carrier_map[gt] = is_carrier
            return is_carrier

    # find record with max GQ, prioritizing non-ref GTs if any exist
    def _record_key(record):
        format_field = record.samples[sample]
        return _is_non_ref(format_field), format_field.get('GQ', -math.inf)
    return max(records, key=_record_key)


def check_multiallelic(records):
    """
    Returns true if any record is multiallelic

    Parameters
    ----------
    records : list of pysam.VariantRecord
    """
    for record in records:
        if record.alts[0] in ['<CNV>', '<CN0>']:
            return True
    return False


def make_multiallelic_alts(records):
    """
    Sets simple symbolic alt for multi-allelic records
    """
    for record in records:
        record.alts = ('<CNV>',)


def update_best_genotypes(new_record, records, preserve_multiallelic=False):
    """
    For each sample in record, update GT and other formats with best genotype

    Parameters
    ----------
    new_record : pysam.VariantRecord
    records : list of SVRecord
    preserve_multiallelic : bool
        If any record is multiallelic, make all genotypes multiallelic
    """

    def _binary_or(x, y):
        return x | y

    def _str_to_tuple(x):
        return (x,) if type(x) is str else x

    if preserve_multiallelic:
        is_multiallelic = check_multiallelic(records)
    else:
        is_multiallelic = False

    if is_multiallelic:
        make_multiallelic_alts(records)

    new_record_formats = new_record.header.formats.keys()
    for sample in new_record.samples:
        best_record = choose_best_genotype(sample, records)
        for key in best_record.format.keys():
            if key in new_record_formats:
                # If any record is multiallelic, replace non-multiallelic
                # genotypes with multiallelic equivalents
                if key == 'GT':
                    gt = best_record.samples[sample][key]
                    if is_multiallelic:
                        # No genotypes for multi-allelic records since they can't always be determined
                        new_record.samples[sample][key] = (None, None)
                    else:
                        new_record.samples[sample][key] = tuple([x if x is None else min(x, 1) for x in gt])
                elif key == 'EV':
                    ev_values = [r.samples[sample][key] for r in records]
                    ev_values = [_str_to_tuple(v) for v in ev_values]
                    ev_types = list(set(type(v) for v in ev_values if v is not None))
                    if len(ev_types) == 0:
                        new_record.samples[sample][key] = None
                    elif len(ev_types) > 1:
                        raise ValueError(f"Found multiple EV types: {str(ev_types)}")
                    else:
                        ev_type = ev_types[0]
                        if ev_type is tuple:
                            new_record.samples[sample][key] = ",".join(
                                sorted(tuple(set(t for ev in ev_values for t in ev))))
                        elif ev_type is int:
                            new_record.samples[sample][key] = reduce(_binary_or, ev_values)
                        else:
                            raise ValueError(f"Unsupported EV type: {str(ev_type)}")
                else:
                    new_record.samples[sample][key] = best_record.samples[sample][key]
