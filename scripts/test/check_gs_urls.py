#!/bin/python

from google.cloud import storage
from google.api_core import exceptions
import json
import argparse
from urllib.parse import urlparse

# Synopsis:
#   This script checks if all GCS URLs in a Cromwell input JSON file exist. URLs that do not exist are printed to stdout.
#
# Requirements:
#   Python >= 3.5
#   Google Cloud Storage Python API
#     - Install with "pip install google-cloud-storage"
#
# Usage:
#   python check_gs_urls.py inputs.json
#
# Parameters:
#   inputs.json : workflow input file
#
# Author: Mark Walker (markw@broadinstitute.org)

# URI scheme for Cloud Storage.
GOOGLE_STORAGE = 'gs'

# Checks if the string is a Google bucket URL


def is_gcs_url(str):
    return urlparse(str).scheme == GOOGLE_STORAGE

# Checks if the object exists in GCS


def check_gcs_url(source_uri, client, project_id):
    def _parse_uri(uri):
        parsed = urlparse(uri)
        bucket_name = parsed.netloc
        bucket_object = parsed.path[1:]
        return bucket_name, bucket_object
    source_bucket_name, source_blob_name = _parse_uri(source_uri)
    source_bucket = client.bucket(source_bucket_name, user_project=project_id)
    source_blob = source_bucket.blob(source_blob_name)
    try:
        if not source_blob.exists():
            print(f"{source_uri} not found")
            return False
    except exceptions.BadRequest as e:
        print(e)
    return True

# Main function


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs_json")
    parser.add_argument("--project-id", required=False,
                        help="Project ID to charge for requester pays buckets")
    args = parser.parse_args()

    with open(args.inputs_json, 'r') as f:
        client = storage.Client()
        inputs = json.load(f)
        for x in inputs:
            if isinstance(inputs[x], str) and is_gcs_url(inputs[x]):
                check_gcs_url(inputs[x], client, args.project_id)
            elif isinstance(inputs[x], list):
                for y in inputs[x]:
                    if is_gcs_url(y):
                        check_gcs_url(y, client, args.project_id)


if __name__ == "__main__":
    main()
