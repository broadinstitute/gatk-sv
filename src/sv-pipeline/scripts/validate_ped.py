#!/bin/python

import argparse
import re


"""
Summary: Validates a PED file for use with GATK-SV. Performs some sample ID validations as well.

Usage: python validate_ped.py -p pedigree.ped -s samples.list

Outputs: The script will write to stdout "PED file passes validation!" if successful, and
    otherwise it will print an error message to stderr describing the PED file format violation.
"""


FIELD_NUMBER_ID = 1
FIELD_NUMBER_SEX = 4
ILLEGAL_ID_SUBSTRINGS = ["chr", "name", "DEL", "DUP", "CPX", "CHROM"]
ID_TYPE_SAMPLE, ID_TYPE_FAMILY, ID_TYPE_PARENT = "sample", "family", "parent"


def validate_id(identifier, id_set, id_type, source_file):
    # check for empty IDs
    if identifier is None or identifier == "":
        raise ValueError(f"Empty {id_type} ID in {source_file}")

    # check all characters are alphanumeric or underscore
    if not re.match(r'^\w+$', identifier):
        raise ValueError(f"Invalid {id_type} ID in {source_file}: '{identifier}'. " +
                         "IDs should only contain alphanumeric and underscore characters.")

    # check for all-numeric IDs
    # except: allow maternal and paternal ID to be 0
    # except: allow all-numeric family IDs (as in current ref panel PED file)
    if id_type != ID_TYPE_FAMILY and not (id_type == ID_TYPE_PARENT and identifier == "0") and identifier.isdigit():
        raise ValueError(f"Invalid {id_type} ID in {source_file}: {identifier}. " +
                         "IDs should not contain only numeric characters.")

    # check for illegal substrings
    for sub in ILLEGAL_ID_SUBSTRINGS:
        if sub in identifier:
            raise ValueError(f"Invalid {id_type} ID in {source_file}: {identifier}. " +
                             f"IDs cannot contain the following substrings: {', '.join(ILLEGAL_ID_SUBSTRINGS)}.")

    # check for duplicate IDs
    if id_set is not None:
        if identifier in id_set:
            raise ValueError(f"Duplicate {id_type} ID in {source_file}: {identifier}")


def get_samples(samples_file):
    # Note that this sample ID validation is incomplete:
    # * It does not check if sample IDs are substrings of other sample IDs
    # * It only checks the provided sample list for duplicates, which is not the full cohort in GatherBatchEvidence
    samples = set()
    with open(samples_file, 'r') as samp:
        for line in samp:
            sample = line.strip()
            validate_id(sample, samples, ID_TYPE_SAMPLE, "sample list")
            samples.add(sample)

    if len(samples) < 1:
        raise ValueError("Empty samples list provided")

    return samples


def validate_ped(ped_file, samples):
    seen_sex_1 = False
    seen_sex_2 = False
    samples_found = set()
    with open(ped_file, 'r') as ped:
        for line in ped:
            # specification allows commented lines, which should be removed by SubsetPedFile and CleanVcfPart1
            if line.startswith("#"):
                continue

            # we require tab-delimited PED files although the specification does not
            fields = line.strip().split("\t")
            if len(fields) != 6:
                raise ValueError("Invalid PED file. PED file must be tab-delimited and have 6 columns: " +
                                 "family_ID, sample_ID, paternal_ID, maternal_ID, sex, phenotype.")

            # validate IDs
            # don't check for duplicates here:
            # family and parent IDs may appear multiple times, and sample IDs checked elsewhere
            for identifier, id_type in zip(fields[:FIELD_NUMBER_SEX],
                                           [ID_TYPE_FAMILY, ID_TYPE_SAMPLE, ID_TYPE_PARENT, ID_TYPE_PARENT]):
                validate_id(identifier, None, id_type, "PED file")

            # check for at least one appearance each of 1 and 2 in sex column
            # this check is to ensure M/F are coded as 1 and 2 (not 0 and 1)
            # and both M and F are present (for sex-specific steps)
            sample_id = fields[FIELD_NUMBER_ID]
            sex = fields[FIELD_NUMBER_SEX]
            if sex == "1":
                seen_sex_1 = True
            elif sex == "2":
                seen_sex_2 = True

            # require sex = 0, 1, or 2
            elif sex != "0":
                raise ValueError(f"Sample {sample_id} has an invalid value for sex: {sex}. " +
                                 "PED file must use the following values for sex: 1=Male, 2=Female, 0=Unknown/Other.")

            # check all samples in samples list are present in PED file exactly once
            # no duplicate sample IDs
            if sample_id in samples_found:
                raise ValueError(f"Invalid PED file. PED file has duplicate entries for sample {sample_id}.")
            elif sample_id in samples:
                samples_found.add(sample_id)
            # if not in samples list, ignore: ok for PED file to contain samples not in sample list

    # check if all samples in sample list are present in PED file
    if len(samples_found) < len(samples):
        missing = samples - samples_found
        raise ValueError(f"Invalid PED file. PED file is missing sample(s): {','.join(missing)}.")

    if not (seen_sex_2 and seen_sex_1):
        raise ValueError("Invalid PED file. PED file must use the following values for sex: " +
                         "1=Male, 2=Female, 0=Unknown/Other. PED file must contain at least " +
                         "one sample with sex=1 and one with sex=2.")

    # passed validation tests
    print("PED file passes validation!")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--ped-file", required=True, help="PED file to validate")
    parser.add_argument("-s", "--samples-file", required=True,
                        help="File containing samples (one per line) that should be in PED file")
    args = parser.parse_args()

    samples = get_samples(args.samples_file)
    validate_ped(args.ped_file, samples)


if __name__ == "__main__":
    main()
