#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import os
import firecloud.api as fapi
from firecloud import fccore
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


def make_other_vcfs_list(namespace, workspace, batches):
    sample_sets = get_sample_sets(namespace, workspace, batches)

    for batch in batches:
        others = [e for e in sample_sets if e['name'] != batch]
        other_pesr_vcfs = [o['attributes']['filtered_pesr_vcf']
                           for o in others]
        other_depth_vcfs = [o['attributes']
                            ['filtered_depth_vcf'] for o in others]

        pesr_fname = '{0}.other_pesr_vcfs.list'.format(batch)
        depth_fname = '{0}.other_depth_vcfs.list'.format(batch)

        fout = open(pesr_fname, 'w')
        for vcf in other_pesr_vcfs:
            fout.write(vcf + '\n')
        fout.close()

        fout = open(depth_fname, 'w')
        for vcf in other_depth_vcfs:
            fout.write(vcf + '\n')
        fout.close()

        yield batch, pesr_fname, depth_fname


def upload_vcf_lists(namespace, workspace, vcf_lists):
    """

    Arguments
    ---------
    list of (str, str)
    """

    bucket_name = get_bucket_name(namespace, workspace)
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    for batch, pesr_vcflist, depth_vcflist in vcf_lists:
        blob = bucket.blob('other_vcf_lists/{0}'.format(pesr_vcflist))
        blob.upload_from_filename(pesr_vcflist)
        blob = bucket.blob('other_vcf_lists/{0}'.format(depth_vcflist))
        blob.upload_from_filename(depth_vcflist)


def update_metadata(namespace, workspace, vcf_lists):
    fpath = 'gs://{0}/other_vcf_lists/{1}'
    bucket_name = get_bucket_name(namespace, workspace)

    for batch, pesr_vcflist, depth_vcflist in vcf_lists:
        list_path = fpath.format(bucket_name, pesr_vcflist)
        update = fapi._attr_set(
            'other_batches_filtered_pesr_vcf_list', list_path)
        r = fapi.update_entity(namespace, workspace,
                               'sample_set', batch, [update])
        fapi._check_response_code(r, 200)

        list_path = fpath.format(bucket_name, depth_vcflist)
        update = fapi._attr_set(
            'other_batches_filtered_depth_vcf_list', list_path)
        r = fapi.update_entity(namespace, workspace,
                               'sample_set', batch, [update])
        fapi._check_response_code(r, 200)


def main():
    proj_required = not bool(fcconfig.project)
    workspace_required = not bool(fcconfig.workspace)

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('batch_list', type=argparse.FileType('r'))
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

    vcf_lists = list(make_other_vcfs_list(
        args.project, args.workspace, batches))

    upload_vcf_lists(args.project, args.workspace, vcf_lists)

    update_metadata(args.project, args.workspace, vcf_lists)

    for batch, pesr_vcflist, depth_vcflist in vcf_lists:
        os.remove(pesr_vcflist)
        os.remove(depth_vcflist)


if __name__ == '__main__':
    main()
