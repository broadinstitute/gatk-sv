#!/bin/python

import sys
from collections import defaultdict


def count_vids(list_path):
	counts = defaultdict(lambda: 0)
	with open(list_path, 'r') as f_list:
		for path in f_list:
			with open(path.strip(), 'r') as f:
				for vid in f:
					counts[vid.strip()] += 1
	return counts


NON_REF_VIDS_LIST = sys.argv[1]
BOTHSIDE_PASS_LIST = sys.argv[2]

non_ref_counts = count_vids(NON_REF_VIDS_LIST)
bothside_pass_counts = count_vids(BOTHSIDE_PASS_LIST)

for vid, bothside_pass_count in bothside_pass_counts.items():
	if bothside_pass_count == 0:
		continue
	non_ref_count = non_ref_counts[vid]
	if non_ref_count == 0:
		continue
	fraction_support = min(1., bothside_pass_count / float(non_ref_count))
	sys.stdout.write("{}\t{}\n".format(fraction_support, vid))
