#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Run before 04_integrate_batches
"""

import argparse
import firecloud.api as fapi
from firecloud import fccore
from firecloud import errors as ferrors
from google.cloud import storage


fcconfig = fccore.config_parse()


def get_bucket_name(namespace, workspace):
    response = fapi.get_workspace(namespace, workspace)
    fapi._check_response_code(response, 200)
    return response.json()['workspace']['bucketName']


def get_sample_sets(namespace, workspace, batches):
    response = fapi.get_entities(namespace, workspace, 'sample_set')
    fapi._check_response_code(response, 200)

    return [entity for entity in response.json() if entity['name'] in batches]


def make_list(sample_sets, attribute, fname):
    fs = [s['attributes'][attribute] for s in sample_sets]

    fout = open(fname, 'w')
    for f in fs:
        fout.write(f + '\n')
    fout.close()


def make_vcf_lists(namespace, workspace, batches):
    sample_sets = get_sample_sets(namespace, workspace, batches)

    make_list(sample_sets, 'genotyped_pesr_vcf', 'cohort.pesr_vcfs.list')
    make_list(sample_sets, 'genotyped_depth_vcf', 'cohort.depth_vcfs.list')
    make_list(sample_sets, 'PE_file', 'cohort.discfile.list')
    make_list(sample_sets, 'PE_file_idx', 'cohort.discfile_idx.list')

    return [('cohort_genotyped_pesr_vcf_list', 'cohort.pesr_vcfs.list'),
            ('cohort_genotyped_depth_vcf_list', 'cohort.depth_vcfs.list'),
            ('cohort_discfile_list', 'cohort.discfile.list'),
            ('cohort_discfile_idx_list', 'cohort.discfile_idx.list')]


def upload_vcf_lists(namespace, workspace, data):
    """

    Arguments
    ---------
    list of (str, str)
    """

    bucket_name = get_bucket_name(namespace, workspace)
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    for attr, fname in data:
        blob = bucket.blob('cohort_lists/{0}'.format(fname))
        blob.upload_from_filename(fname)


def update_metadata(namespace, workspace, sample_set, data):
    fpath = 'gs://{0}/cohort_lists/{1}'
    bucket_name = get_bucket_name(namespace, workspace)

    for attr, fname in data:
        list_path = fpath.format(bucket_name, fname)
        update = fapi._attr_set(attr, list_path)
        try:
            r = fapi.update_entity(namespace, workspace,
                                   'sample_set', sample_set, [update])
            fapi._check_response_code(r, 200)
        except ferrors.FireCloudServerError:
            pass


def main():
    proj_required = not bool(fcconfig.project)
    workspace_required = not bool(fcconfig.workspace)

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('batch_list', type=argparse.FileType('r'))
    parser.add_argument('master_sample_set')
    parser.add_argument('-w', '--workspace', default=fcconfig.workspace,
                        required=workspace_required,
                        help='Workspace name. Required if no default '
                        'workspace configured.')
    parser.add_argument('-p', '--project', default=fcconfig.project,
                        required=proj_required,
                        help='Project (workspace namespace). Required if no '
                        'default project configured.')
    args = parser.parse_args()

    batches = [l.strip() for l in args.batch_list.readlines()]

    data = make_vcf_lists(args.project, args.workspace, batches)

    upload_vcf_lists(args.project, args.workspace, data)

    update_metadata(args.project, args.workspace, args.master_sample_set, data)

    #  os.remove(pesr_vcflist)
    #  os.remove(depth_vcflist)


if __name__ == '__main__':
    main()
