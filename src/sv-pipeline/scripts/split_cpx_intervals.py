import pysam
import sys
import argparse


INS_TYPE_CPX = {'dDUP', 'dDUP_iDEL', 'INS_iDEL'}


def copy_record(rec, out, counter):
    record = out.new_record()
    record.info['ORIGINAL_VID'] = rec.id
    record.ref = rec.ref
    record.alts = rec.alts
    record.qual = rec.qual
    for key in rec.info:
        record.info[key] = rec.info[key]

    if rec.info['CPX_TYPE'] in INS_TYPE_CPX:
        record.info['SINK_CHROM'] = rec.chrom
        record.info['SINK_POS'] = rec.start
        record.info['SINK_END'] = rec.stop

    for filt in rec.filter:
        record.filter.add(filt)

    record.id = f"{rec.id}_{counter}"

    return record



def process(vcf_in, vcf_out):
    with pysam.VariantFile(vcf_in) as vcf:
        header = vcf.header
        header.add_line('##INFO=<ID=ORIGINAL_VID,Number=1,Type=String,Description="Original VID for split CPX variants">')
        header.add_line('##INFO=<ID=SINK_CHROM,Number=1,Type=String,Description="Chromosome for sink for split dDUP variants">')
        header.add_line('##INFO=<ID=SINK_POS,Number=1,Type=Integer,Description="Position for sink for split dDUP variants">')
        header.add_line('##INFO=<ID=SINK_END,Number=1,Type=Integer,Description="Position for sink for split dDUP variants">')

        with pysam.VariantFile(vcf_out, 'w', header=vcf.header) as out:
            for rec in vcf:
                if rec.info['SVTYPE'] == "CPX":
                    counter = 1

                    # output a record for the sink
                    if rec.info['CPX_TYPE'] in INS_TYPE_CPX:
                        record = copy_record(rec, out, counter)
                        record.chrom = rec.chrom
                        record.start = rec.start
                        record.stop = rec.stop
                        record.info['SVTYPE'] = 'INS'
                        record.alts = ("<INS>",)

                        if record.info['CPX_TYPE'] in ['dDUP_iDEL', 'INS_iDEL']:
                            svlen = record.info['SVLEN']
                            deletion = [x for x in record.info['CPX_INTERVALS'] if x.startswith("DEL_")][0].split(":")[1].split("-")
                            del_len = int(deletion[1]) - int(deletion[0])
                            record.info['SVLEN'] = svlen - del_len  # SVLEN for iDELs is source + deletion length
                        out.write(record)
                        counter += 1

                    # output a record for every interval
                    for interval in rec.info['CPX_INTERVALS']:
                        record = copy_record(rec, out, counter)
                        svtype, coords = interval.split("_")
                        chrom, breakpoints = coords.split(":")
                        pos, end = [int(x) for x in breakpoints.split("-")]
                        if svtype == "INS":
                            svtype = "DUP"  # for INS intervals in INS_iDEL, allow match with DUP sources
                        record.info['SVTYPE'] = svtype
                        record.alts = (f"<{svtype}>",)
                        record.chrom = chrom
                        record.start = pos  # copy the record in order to avoid nastiness with start/end
                        record.stop = end  # note: had to set symbolic alt in order for END to show up
                        record.info['SVLEN'] = end - pos
                        out.write(record)
                        counter += 1


def _parse_arguments(argv):
    parser = argparse.ArgumentParser(
        description="Create VID match and score table for input to SVFederate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True,
                        help="VCF containing CPX records to extract and split")
    parser.add_argument("--out", type=str, required=True,
                        help="Output table")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main():
    args = _parse_arguments(sys.argv)

    process(args.vcf, args.out)


if __name__ == "__main__":
    main()
