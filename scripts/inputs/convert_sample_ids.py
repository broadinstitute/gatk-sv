#!/usr/bin/env python

"""
Converts external sample IDs to safe GATK-SV IDs

"""

import argparse
import hashlib
import re
import sys


DEFAULT_HASH_LENGTH = 6
ID_FORMAT = "__{id_name:s}__{id_hash:s}"


def convert_id(sample_id, hash_input, hash_length):
    if sample_id is None or sample_id == "":
        raise ValueError("Encountered None-type or empty id")
    if "__" in sample_id:
        raise ValueError("Encountered double-underscore in sample id: {:s}".format(sample_id))
    hash_fn = hashlib.sha1()
    hash_fn.update(hash_input.encode(encoding='UTF-8', errors='strict'))
    id_hash = hash_fn.hexdigest()[:hash_length]
    id_no_special_chars = re.sub('[^0-9a-zA-Z]+', '_', sample_id).lower()
    return ID_FORMAT.format(id_name=id_no_special_chars, id_hash=id_hash)


def convert_ids_list(external_ids, hash_input, hash_length):
    if len(external_ids) != len(hash_input):
        raise ValueError("There were {:d} external ids but {:d} hash inputs".format(len(external_ids), len(hash_input)))
    zipped = zip(external_ids, hash_input)
    return [convert_id(sample_id=z[0], hash_input=z[1], hash_length=hash_length) for z in zipped]


def read_list_file(path):
    with open(path, 'r') as f:
        entries = f.read().splitlines()
    if len(entries) == 0:
        raise ValueError("List empty: {}".format(path))
    return entries


def test_ids(external_ids, converted_ids, hash_input, skip_substring):
    num_external_ids = len(external_ids)
    num_converted_ids = len(converted_ids)
    if num_external_ids != num_converted_ids:
        raise ValueError("Number of external ids was {:d} but there were {:d} converted ids"
                         .format(num_external_ids, num_converted_ids))

    num_hash_inputs = len(hash_input)
    num_unique_hash_inputs = len(set(hash_input))
    if num_hash_inputs != num_unique_hash_inputs:
        raise ValueError("{:d} hash inputs were provided but only {:d} were unique".format(num_hash_inputs, num_unique_hash_inputs))

    num_unique_converted_ids = len(set(converted_ids))
    if num_unique_converted_ids != num_converted_ids:
        raise ValueError("There are {:d} converted ids but only {:d} are unique".format(num_converted_ids, num_unique_converted_ids))

    # Check if any converted ID is a substring of another (slow naive approach)
    if not skip_substring:
        for i in range(num_converted_ids):
            other_ids = " ".join(converted_ids[:i] + converted_ids[(i + 1):])
            if converted_ids[i] in other_ids:
                for j in range(num_converted_ids):
                    if j != i and converted_ids[i] in converted_ids[j]:
                        raise ValueError("Conflicting substring in converted ids: \"{:s}\" \"{:s}\"".format(converted_ids[i], converted_ids[j]))
            if (i + 1) % 1000 == 0:
                sys.stderr.write("Checked for substring collisions in {} / {} converted ids\n".format(i + 1, num_converted_ids))


def write_ids(converted_ids):
    for i in converted_ids:
        print(i)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sample_list', type=str, help="Newline-delimited list of external sample IDs")
    parser.add_argument('--hash-input', type=str, required=False,
                        help="Newline-delimited list of corresponding unique strings to hash [default sample_list]")
    parser.add_argument('--hash-length', type=int, required=False, default=DEFAULT_HASH_LENGTH,
                        help="Appended hash length in characters [default {:d}]".format(DEFAULT_HASH_LENGTH))
    parser.add_argument('--skip-substring-check', action='store_true',
                        help="Skip converted id substring check, which is slow")
    args = parser.parse_args()

    external_ids = read_list_file(args.sample_list)
    if args.hash_input is not None:
        hash_input = read_list_file(args.hash_input)
    else:
        hash_input = external_ids
    converted_ids = convert_ids_list(external_ids=external_ids, hash_input=hash_input, hash_length=args.hash_length)
    test_ids(external_ids=external_ids, converted_ids=converted_ids, hash_input=hash_input, skip_substring=args.skip_substring_check)
    write_ids(converted_ids)


if __name__ == '__main__':
    main()
