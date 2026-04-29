import pysam
import gzip
import argparse
import sys


CONSEQUENCES = 'PREDICTED_LOF PREDICTED_COPY_GAIN PREDICTED_INTRAGENIC_EXON_DUP PREDICTED_PARTIAL_EXON_DUP PREDICTED_TSS_DUP PREDICTED_DUP_PARTIAL PREDICTED_INV_SPAN PREDICTED_UTR PREDICTED_INTRONIC PREDICTED_PROMOTER'.split()
FREQUENCIES = 'rare ultra_rare singleton rare_large common'.split()


def process_vcf(vcf_path):
    gene_dict = dict()
    with pysam.VariantFile(vcf_path, 'r') as vcf:
        for record in vcf:
            if (len([x for x in record.filter if x!= "PASS"]) > 0) or record.info['SVTYPE'] == 'CNV':
                continue  # PASS only
            for key in CONSEQUENCES:
                if key in record.info:
                    for gene in record.info[key]:
                        if gene not in gene_dict:
                            gene_dict[gene] = dict()
                        if key not in gene_dict[gene]:
                            gene_dict[gene][key] = {x: 0 for x in FREQUENCIES}
                        af = record.info['AF'][0]
                        if af < 0.01:
                            gene_dict[gene][key]['rare'] += 1
                            if 'SVLEN' in record.info and record.info['SVLEN'] > 100000:
                                gene_dict[gene][key]['rare_large'] += 1
                            if af < 0.001:
                                gene_dict[gene][key]['ultra_rare'] += 1
                            if record.info['AC'] == 1:
                                gene_dict[gene][key]['singleton'] += 1
                        else:
                            gene_dict[gene][key]['common'] += 1
    return gene_dict


def write_gene_data(gene_info_path, out_path, gene_dict):
    with gzip.open(gene_info_path, 'rt') as inp, gzip.open(out_path, 'wt') as out:
        first = True
        for line in inp:
            fields = line.strip('\n').split('\t')
            if first:
                header = dict([(x,i) for i,x in enumerate(fields)])
                out_header = fields
                for conseq in CONSEQUENCES:
                    for freq in FREQUENCIES:
                        out_header.append(f"{conseq}.{freq}")
                out.write("\t".join(out_header) + "\n")
                first = False
                continue
            gene_name = fields[header['gene_name']]
            gene_data = fields
            if gene_name not in gene_dict:
                gene_data.extend(["0"]*(len(CONSEQUENCES)*len(FREQUENCIES)))
            else:
                for conseq in CONSEQUENCES:
                    if conseq in gene_dict[gene_name]:
                        for freq in FREQUENCIES:
                            gene_data.append(str(gene_dict[gene_name][conseq][freq]))
                    else:
                        gene_data.extend(["0"]*len(FREQUENCIES))
            out.write("\t".join(gene_data) + "\n")


def _parse_arguments(argv):
    parser = argparse.ArgumentParser(
        description="Create VID match and score table for input to SVFederate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True,
                        help="VCF with allele frequency information and SVAnnotate annotations")
    parser.add_argument("--gene-info", type=str, required=True,
                        help="Gene info table (gzipped)")
    parser.add_argument("--out", type=str, required=True,
                        help="Output table name")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main():
    args = _parse_arguments(sys.argv)

    gene_dict = process_vcf(args.vcf)
    write_gene_data(args.gene_info, args.out, gene_dict)


if __name__ == "__main__":
    main()

