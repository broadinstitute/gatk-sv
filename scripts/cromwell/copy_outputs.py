#!/bin/python

import sys
import json
import os.path
import argparse
from google.cloud import storage

# Synopsis:
#  Copies workflow outputs needed for downstream processing to a destination bucket.
#
# Author: Mark Walker (markw@broadinstitute.org)

# Output file definitions: FILENAME_MAP[workflow_name][output_variable] = destination_file_suffix
FILENAME_MAP = {
    'GATKSVPipelinePhase1': {
        'filtered_depth_vcf': 'filtered_depth_vcf.vcf.gz',
        'filtered_pesr_vcf': 'filtered_pesr_vcf.vcf.gz',
        'cutoffs': 'rf_cutoffs.tsv',
        'outlier_samples_excluded_file': 'outliers.list',
        'batch_samples_postOutlierExclusion_file': 'filtered_samples.list',
        'ped_file_postOutlierExclusion': 'filtered.ped',
        'merged_SR': 'SR.txt.gz',
        'merged_SR_index': 'SR.txt.gz.tbi',
        'merged_PE': 'PE.txt.gz',
        'merged_PE_index': 'PE.txt.gz.tbi',
        'merged_bincov': 'RD.txt.gz',
        'merged_bincov_index': 'RD.txt.gz.tbi',
        'median_cov': 'median_cov.bed'
    },
    'MergeCohortVcfs': {
        'cohort_pesr_vcf': 'cohort_pesr.vcf.gz',
        'cohort_depth_vcf': 'cohort_depth.vcf.gz',
        'cohort_combined': 'cohort.combined.bed',
        'lookup': 'master_cluster_dups.bed',
        'cohort_sort': 'cohort.sort.bed',
        'cluster_combined': 'cluster.combined.bed'
    },
    'Module04': {
        'sr_bothside_pass': 'sr_bothside_pass.txt',
        'sr_background_fail': 'sr_background_fail.txt',
        'trained_PE_metrics': 'trained_PE_metrics.txt',
        'trained_SR_metrics': 'trained_SR_metrics.txt',
        'trained_genotype_pesr_pesr_sepcutoff': 'trained_genotype_pesr_pesr_sepcutoff.txt',
        'trained_genotype_pesr_depth_sepcutoff': 'trained_genotype_pesr_depth_sepcutoff.txt',
        'trained_genotype_depth_pesr_sepcutoff': 'trained_genotype_depth_pesr_sepcutoff.txt',
        'trained_genotype_depth_depth_sepcutoff': 'trained_genotype_depth_depth_sepcutoff.txt',
        'genotyped_depth_vcf': 'genotyped.depth.vcf.gz',
        'genotyped_depth_vcf_index': 'genotyped.depth.vcf.gz.tbi',
        'genotyped_pesr_vcf': 'genotyped.pesr.vcf.gz',
        'genotyped_pesr_vcf_index': 'genotyped.pesr.vcf.gz.tbi',
        'regeno_depth': 'regeno_depth.bed'
    }
}


def get_uris(metadata, output_name, dest_prefix):
    if 'workflowName' not in metadata:
        raise ValueError("Workflow name not found. Check metadata file.")
    workflow = metadata['workflowName']
    if workflow not in FILENAME_MAP:
        raise ValueError(f"Unknown workflow {workflow}")
    outputs = metadata['outputs']
    for var in FILENAME_MAP[workflow]:
        key = f"{workflow}.{var}"
        if outputs[key] is not None:
            source_uri = outputs[key]
            dest_filename = f"{output_name}.{FILENAME_MAP[workflow][var]}"
            dest_uri = os.path.join(dest_prefix, dest_filename)
            yield (source_uri, dest_uri)


def copy_blob(storage_client, bucket_name, blob_name, destination_bucket_name, destination_blob_name):
    source_bucket = storage_client.bucket(bucket_name)
    source_blob = source_bucket.blob(blob_name)
    destination_bucket = storage_client.bucket(destination_bucket_name)
    destination_blob = destination_bucket.blob(destination_blob_name)
    source_uri = f"gs://{source_bucket.name}/{source_blob.name}"
    destination_uri = f"gs://{destination_bucket.name}/{destination_blob_name}"
    if destination_blob.exists():
        sys.stderr.write(
            f"Target {destination_uri} exists, cautiously refusing to overwrite. Aborting...\n")
        sys.exit(1)
    sys.stderr.write(f"Copying {source_uri}...")
    (token, bytes_rewritten, total_bytes) = destination_blob.rewrite(source=source_blob)
    while token is not None:
        (token, bytes_rewritten, total_bytes) = destination_blob.rewrite(
            source=source_blob, token=token)
    size_kb = int(bytes_rewritten / 1024)
    sys.stderr.write(f"done ({size_kb} KB)\n")


def copy_uri(source_uri, dest_uri, storage_client):
    def _parse_uri(uri):
        tokens = uri.split('/')
        bucket_name = tokens[2]
        bucket_object = '/'.join(tokens[3:])
        return bucket_name, bucket_object
    source_bucket_name, source_blob_name = _parse_uri(source_uri)
    dest_bucket_name, dest_blob_name = _parse_uri(dest_uri)
    copy_blob(storage_client, source_bucket_name,
              source_blob_name, dest_bucket_name, dest_blob_name)


# Main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", help="Batch or cohort name", required=True)
    parser.add_argument(
        "--metadata", help="Workflow metadata JSON file", required=True)
    parser.add_argument(
        "--dest", help="Destination GCS URI (e.g. \"gs://my-bucket/output\")", required=True)
    args = parser.parse_args()
    metadata = json.load(open(args.metadata, 'r'))
    output_uris = get_uris(metadata, args.name, args.dest)
    client = storage.Client()
    for source_uri, dest_uri in output_uris:
        copy_uri(source_uri, dest_uri, client)


if __name__ == "__main__":
    main()
