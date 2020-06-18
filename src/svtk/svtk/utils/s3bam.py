#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 ec2-user <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""
Load S3-hosted bam into pysam.AlignmentFile
"""

import os
import boto3
import pysam


def load_s3bam(bam_path, index_dir=None):
    if not bam_path.startswith('s3://'):
        raise Exception('Bam {0} is not a valid S3 path'.format(bam_path))

    # Pysam doesn't accept explicit path to index file, expects index to be
    # present in working directory. If a local copy of the index is available,
    # move to its directory to use it.
    # Otherwise, index is downloaded automatically
    if index_dir is not None:
        os.chdir(index_dir)
    else:
        msg = ('Local index directory not specified for {0}. Downloading '
               'remote copy of index to working directory.')
        raise Warning(msg.format(bam_path))

    # Parse bucket and key from filepath
    s3path = bam_path[5:]
    bucket = s3path.split('/')[0]
    bam_path = '/'.join(s3path.split('/')[1:])

    # Create S3 client and get presigned URL
    # Necessary to take advantage of pysam's https support until the library
    # supports S3 paths directly
    s3 = boto3.client('s3')
    url = s3.generate_presigned_url(
            ClientMethod='get_object',
            Params={'Bucket': bucket, 'Key': bam_path},
            ExpiresIn=86400)

    return pysam.AlignmentFile(url)
