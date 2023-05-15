#!/bin/python

import sys
import json
import os.path
import argparse
from google.cloud import storage

# Synopsis:
#  Copies workflow outputs needed for downstream processing to a destination bucket.
# Author: Mark Walker (markw@broadinstitute.org)


def get_uris(metadata, dest_prefix):
    if 'workflowName' not in metadata:
        raise ValueError("Workflow name not found. Check metadata file.")
    outputs = metadata['outputs']
    for source_uri in outputs.values():
        if source_uri is not None and source_uri.startswith("gs://"):
            dest_filename = os.path.basename(source_uri)
            dest_uri = os.path.join(dest_prefix, dest_filename)
            yield source_uri, dest_uri


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
    parser.add_argument("--metadata", help="Workflow metadata JSON file", required=True)
    parser.add_argument("--dest", help="Destination GCS URI (e.g. \"gs://my-bucket/output\")", required=True)
    args = parser.parse_args()
    metadata = json.load(open(args.metadata, 'r'))
    output_uris = get_uris(metadata, args.dest)
    client = storage.Client()
    for source_uri, dest_uri in output_uris:
        copy_uri(source_uri, dest_uri, client)


if __name__ == "__main__":
    main()
