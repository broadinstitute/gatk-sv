#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Wrap bgzip output
"""

import sys
import subprocess


class BgzipFile:
    def __init__(self, filename, bgzip=True):
        self.filename = filename
        self.bgzip = bgzip

    def __enter__(self):
        # Open output file
        if self.is_stdout:
            self.file = sys.stdout.buffer
        else:
            self.file = open(self.filename, 'wb')

        if self.bgzip:
            if not self.filename.endswith('gz'):
                msg = 'Bgzipped filename "{0}" does not end with .gz'
                raise Exception(msg.format(self.filename))

            self.pipe = subprocess.Popen(['bgzip', '-c'],
                                         stdin=subprocess.PIPE,
                                         stdout=self.file)
            self.file = self.pipe.stdin

        return self.file

    def __exit__(self, *args):
        if self.bgzip:
            stdout, stderr = self.pipe.communicate()

            if not self.is_stdout:
                tabix = 'tabix -f -s1 -b2 -e2 %s' % self.filename
                subprocess.call(tabix.split())

        else:
            self.file.close()

    @property
    def is_stdout(self):
        return self.filename in '- stdout'.split()
