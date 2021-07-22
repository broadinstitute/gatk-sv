import hail as hl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('out_bucket')
parser.add_argument('cluster_name')
args = parser.parse_args()

files = [f.rstrip() for f in open("files.list", "r").readlines()]

all_datasets = [hl.import_vcf(f, reference_genome='GRCh38', force_bgz=True) for f in files]

mt = hl.MatrixTable.union_rows(*all_datasets)
# rest the qual to missing because hail by default populates it with -1.00e+01
merged_reset_qual = mt.annotate_rows(qual=hl.missing('float64'))

hl.export_vcf(merged_reset_qual,
              "gs://{}/{}/merged.vcf.bgz".format(args.out_bucket, args.cluster_name),
              metadata=hl.get_vcf_metadata(files[0]))
