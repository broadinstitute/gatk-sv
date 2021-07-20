import hail as hl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('out_bucket')
parser.add_argument('cluster_name')
args = parser.parse_args()

files = [f.rstrip() for f in open("files.list", "r").readlines()]

all_datasets = [hl.import_vcf(f, reference_genome='GRCh38', force_bgz=True) for f in files]

merged = hl.MatrixTable.union_rows(*all_datasets)

hl.export_vcf(merged, "gs://{}/{}/merged.vcf.bgz".format(args.out_bucket, args.cluster_name), metadata=hl.get_vcf_metadata(files[0]))

