import argparse
from pysam import VariantFile


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_vcf')
    parser.add_argument('shard_prefix')
    parser.add_argument('records_per_shard')
    args = parser.parse_args()

    records_per_shard = int(args.records_per_shard)

    record_idx = 0
    current_shard = -1
    out_vcf = None
    input_vcf = VariantFile(args.input_vcf, 'r')
    for record in input_vcf:
        shard = int(record_idx / records_per_shard)
        if shard != current_shard:
            if out_vcf is not None:
                out_vcf.close()
            out_vcf = VariantFile("{}.{:05d}.vcf.gz".format(args.shard_prefix, shard), 'w', header=input_vcf.header)
        current_shard = shard
        out_vcf.write(record)
        record_idx = record_idx + 1
    out_vcf.close()
    input_vcf.close()


if __name__ == '__main__':
    main()