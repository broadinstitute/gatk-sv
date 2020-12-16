import hail as hl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('out_bucket')
parser.add_argument('cluster_name')
args = parser.parse_args()

files = [f.rstrip() for f in open("files.list", "r").readlines()]

# Define custom reference with only primary contigs, otherwise Hail adds all GRCh38 contigs,
# which may be problematic downstream

contigs = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY"
]

lengths = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415
}

ref = hl.ReferenceGenome(name="hg38", contigs=contigs, lengths=lengths, x_contigs="chrX", y_contigs="chrY")
all_datasets = hl.import_vcf(files, reference_genome=ref, force_bgz=True)

# union_rows approach causes ClassTooLargeException
#mt = hl.MatrixTable.union_rows(*all_datasets)
mt = all_datasets
# rest the qual to missing because hail by default populates it with -1.00e+01
merged_reset_qual = mt.annotate_rows(qual=hl.missing('float64'))

hl.export_vcf(merged_reset_qual,
              "gs://{}/{}/merged.vcf.bgz".format(args.out_bucket, args.cluster_name),
              metadata=hl.get_vcf_metadata(files[0]))
