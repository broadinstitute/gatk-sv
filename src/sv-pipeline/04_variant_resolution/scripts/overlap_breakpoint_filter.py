#!/bin/python

import pysam
import sys
from collections import defaultdict
from operator import attrgetter
from itertools import combinations

import svtk.utils as svu

# Inputs
VCF_PATH = sys.argv[1]
BOTHSIDE_PASS_PATH = sys.argv[2]
BACKGROUND_FAIL_PATH = sys.argv[3]
DROPPED_RECORD_OUTPUT_VCF_PATH = sys.argv[4]
if len(sys.argv) >= 6:
    DEBUG_OUTPUT_PATH = sys.argv[5]
else:
    DEBUG_OUTPUT_PATH = None

# Sorts list xs by specified attributes
def multisort(xs, specs):
    for key, reverse in reversed(specs):
        xs.sort(key=attrgetter(key), reverse=reverse)
    return xs


# Container for record data used for duplicate filtering
class RecordData:
    def __init__(self, record):
        self.id = record.id
        if 'EVIDENCE' in record.info:
            ev = set(record.info['EVIDENCE'])
        else:
            ev = set()
        if 'PE' in ev and 'SR' in ev and 'RD' in ev:
            self.level_of_support = 1
        elif 'PE' in ev and 'RD' in ev:
            self.level_of_support = 2
        elif 'PE' in ev and 'SR' in ev:
            self.level_of_support = 3
        elif 'RD' in ev and 'SR' in ev:
            self.level_of_support = 4
        elif 'PE' in ev:
            self.level_of_support = 5
        elif 'RD' in ev:
            self.level_of_support = 6
        elif 'SR' in ev:
            self.level_of_support = 7
        elif len(ev) == 0:
            self.level_of_support = 8
        else:
            raise ValueError("Uninterpretable evidence: {}".format(ev))
        if record.id in bothside_pass:
            self.both_end_support = bothside_pass[record.id]
        else:
            self.both_end_support = 0
        self.sr_fail = record.id in background_fail
        self.is_bnd = record.info['SVTYPE'] == 'BND'
        self.vargq = record.info['varGQ']
        self.called_samples = [sample_id_to_idx[s] for s in svu.get_called_samples(record)]
        self.freq = len(self.called_samples)
        self.length = record.info['SVLEN']
        self.gt_50bp = self.length >= 50
        self.is_mei = 'melt' in record.info['ALGORITHMS']

    def __str__(self):
        return ",".join(str(x) for x in
                        (self.is_bnd, self.level_of_support, self.is_mei, self.both_end_support,
                         self.sr_fail, self.vargq, self.freq, self.gt_50bp, self.length, self.id))


vcf = pysam.VariantFile(VCF_PATH)
sample_id_to_idx = {s: i for i, s in enumerate(vcf.header.samples)}

# Get breakend-to-record-id(s) map
sys.stderr.write("Gathering breakends...\n")
bnds_to_ids = defaultdict(list)
for record in vcf:
    if (record.info['SVTYPE'] == 'DEL' or record.info['SVTYPE'] == 'DUP') and record.info['SVLEN'] >= 5000:
        continue
    strands = record.info['STRANDS']
    bnd1 = "{}_{}_{}".format(record.chrom, record.pos, strands[0])
    bnd2 = "{}_{}_{}".format(record.info['CHR2'], record.stop, strands[1])
    bnds_to_ids[bnd1].append(record.id)
    bnds_to_ids[bnd2].append(record.id)

vcf.close()
# Filter down to breakends associated with at least 2 records
dup_bnds_to_ids = {bnd: bnds_to_ids[bnd] for bnd in bnds_to_ids if len(bnds_to_ids[bnd]) > 1}
sys.stderr.write("Found {} breakends.\n".format(len(dup_bnds_to_ids)))
del bnds_to_ids

# Reverse map of record-id-to-breakend(s)
ids_to_bnds = defaultdict(list)
for bnd, ids in dup_bnds_to_ids.items():
    for record_id in ids:
        ids_to_bnds[record_id].append(bnd)


# Read bothside-pass/background-fail records
sys.stderr.write("Reading SR files...\n")
with open(BOTHSIDE_PASS_PATH) as f:
    bothside_pass = {line.strip().split('\t')[-1]: float(line.strip().split('\t')[0]) for line in f}

with open(BACKGROUND_FAIL_PATH) as f:
    background_fail = set([line.strip().split('\t')[-1] for line in f])

# Store duplicate record data in memory
sys.stderr.write("Reading putative duplicate record data...\n")
vcf = pysam.VariantFile(VCF_PATH)
bnds_to_data = defaultdict(list)
num_samples = len(vcf.header.samples)
for record in vcf:
    if record.id in ids_to_bnds:
        for bnd in ids_to_bnds[record.id]:
            bnds_to_data[bnd].append(RecordData(record))
vcf.close()

# From the lists (each possibly containing more than 2 records), convert to all possible pairwise associations
pairwise_record_data = []
for data_list in bnds_to_data.values():
    pairwise_record_data.extend(combinations(data_list, 2))

sys.stderr.write("Found {} possible duplicate pairs.\n".format(len(pairwise_record_data)))
del bnds_to_data

# This is how we sort record pairs to determine which one gets filtered
sort_spec = [
    ('is_bnd', False),
    ('level_of_support', False),
    ('is_mei', True),
    ('both_end_support', True),
    ('sr_fail', False),
    ('vargq', True),
    ('freq', True),
    ('gt_50bp', False),
    ('length', False),
    ('id', False)
]

# Iterate through record pairs and generate list of record ids to filter out
ids_to_remove_dict = dict()
if DEBUG_OUTPUT_PATH is not None:
    debug = open(DEBUG_OUTPUT_PATH, 'w')
    debug.write("#record_kept\trecord_dropped\n")
for data_list in pairwise_record_data:
    # Check for 50% sample overlap
    sample_intersection = set(data_list[0].called_samples).intersection(data_list[1].called_samples)
    max_freq = max(data_list[0].freq, data_list[1].freq)
    if len(sample_intersection) < 0.50 * max_freq:
        continue
    # Determine which to filter
    sorted_data_list = multisort(list(data_list), sort_spec)
    ids_to_remove_dict[sorted_data_list[1].id] = sorted_data_list[0].id
    if DEBUG_OUTPUT_PATH is not None:
        debug.write("\t".join(str(x) for x in sorted_data_list) + "\n")
if DEBUG_OUTPUT_PATH is not None:
    debug.close()

# Perform filtering
sys.stderr.write("Filtering {} records\n".format(len(ids_to_remove_dict)))
vcf = pysam.VariantFile(VCF_PATH)
header = vcf.header
sys.stdout.write(str(header))

# Create
header.add_line(
    '##INFO=<ID=BPID,Number=.,Type=String,'
    'Description="ID of retained variant from breakpoint overlap filtering">')
dropped_record_vcf = pysam.VariantFile(DROPPED_RECORD_OUTPUT_VCF_PATH, 'w', header=header)

for record in vcf:
    if record.id in ids_to_remove_dict:
        record.info['BPID'] = ids_to_remove_dict[record.id]
        dropped_record_vcf.write(record)
    else:
        sys.stdout.write(str(record))
vcf.close()
dropped_record_vcf.close()
