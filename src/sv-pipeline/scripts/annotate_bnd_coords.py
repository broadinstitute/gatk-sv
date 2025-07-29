import pysam
import argparse


def annotate_bnd_coords(input_vcf_path, output_vcf_path):
    with pysam.VariantFile(input_vcf_path, 'r') as vcf_in:
        header = vcf_in.header

        with pysam.VariantFile(output_vcf_path, 'w', header=header) as vcf_out:
            for record in vcf_in:
                if record.info['SVTYPE'] == 'BND' and 'END2' not in record.info:
                    record.info['END2'] = record.stop
                    record.stop = record.pos + 1
                vcf_out.write(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update coordinates for BND variants in a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file path.")
    parser.add_argument("output_vcf", help="Output VCF file path.")
    args = parser.parse_args()

    annotate_bnd_coords(args.input_vcf, args.output_vcf)
