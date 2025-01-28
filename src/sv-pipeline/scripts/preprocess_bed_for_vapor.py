import gzip
import argparse
import logging


"""
Reformats an SV BED file for VaPoR.
* Outputs headerless BED file containing columns chrom, start, end, name, svtype
* Extracts contig of interest
* Keeps DEL, DUP, INS, INV, and CTX but drops BND, CPX, and CNV types
* Appends SVLEN to INS
* Extracts sites in sample of interest

Input BED file:
* Must be sorted by contig
* First 5 columns: chrom, start, end, name, svtype
* Header with SVLEN column required to reformat INS with SVLEN
* If sample provided, uses header to find samples column.
  If header is absent, assumes samples are in 6th column if present.
"""


REMOVE_TYPES = {"BND", "CPX", "CNV"}


def is_gzipped(path):
    return path.endswith(".gz")


def handle_header(line, columns, fields, default_num_columns, sample_to_extract):
    if line.startswith("#"):
        # get column names beyond first default ones from header if available
        for i, name in enumerate(fields[default_num_columns:]):
            columns[name] = default_num_columns + i
        if "SVLEN" not in columns:
            raise ValueError("SVLEN column not found in header")
    else:
        raise ValueError("Header not found. Header must exist and start with #")
        if len(fields) >= default_num_columns:
            columns['samples'] = default_num_columns  # if no header but extra fields, assume samples is next column
    if sample_to_extract is not None and "samples" not in columns:
        raise ValueError("Sample to extract provided but no samples column found")


def reformat(bed_in, bed_out, contig, sample_to_extract):
    open_fn = gzip.open if is_gzipped(bed_in) else open
    open_mode = 'rt' if is_gzipped(bed_in) else 'r'

    default_columns = "chrom start end name svtype".split()
    default_num_columns = len(default_columns)
    columns = dict(zip(default_columns, range(default_num_columns)))
    first = True
    found_contig = False
    found_sample = False
    with open_fn(bed_in, open_mode) as inp, open(bed_out, 'w') as out:
        for line in inp:
            fields = line.lstrip("#").rstrip('\n').split('\t')
            if first:
                handle_header(line, columns, fields, default_num_columns, sample_to_extract)
            first = False
            if fields[columns["chrom"]] != contig:
                if found_contig:
                    break  # stop if gone past the contig of interest (optimization for early contigs)
                continue  # extract only contig of interest. also drops header if exists
            found_contig = True
            if fields[columns["svtype"]] in REMOVE_TYPES:
                continue  # drop BND, CNV, CPX
            if sample_to_extract is not None and "samples" in columns and \
                    sample_to_extract not in fields[columns["samples"]]:
                continue  # extract events in sample of interest if provided
            found_sample = True
            svtype_write = fields[columns["svtype"]]
            # for INS, format svtype as INS_SVLEN for VaPoR if SVLEN column in input BED
            if ("INS" in svtype_write or "MEI" in svtype_write) and "SVLEN" in columns:
                svtype_write = f"INS_{fields[columns['SVLEN']]}"
            # write chrom, pos, end, SVID, svtype/description for vapor
            out.write("\t".join([fields[columns["chrom"]], fields[columns["start"]], fields[columns["end"]],
                                 fields[columns["name"]], svtype_write]) + "\n")

    if not found_contig:
        logging.warning(f"Could not find any records for contig {contig} in the provided BED file")
    if sample_to_extract is not None and "samples" in columns and not found_sample:
        logging.warning(f"Could not find any records for sample {sample_to_extract} in the provided BED file")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--contig", help="Contig to extract")
    parser.add_argument("--bed-in", help="Input SV BED file")
    parser.add_argument("--bed-out", help="Output Vapor-formatted SV BED file")
    parser.add_argument("-s", "--sample",
                        help="Sample to extract (if input BED file is multi-sample)",
                        required=False)
    parser.add_argument("-l", "--log-level", required=False, default="INFO",
                        help="Specify level of logging information, ie. info, warning, error (not case-sensitive). "
                             "Default: INFO")
    args = parser.parse_args()

    # Set logging level from -l input
    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')

    reformat(args.bed_in, args.bed_out, args.contig, args.sample)


if __name__ == "__main__":
    main()
