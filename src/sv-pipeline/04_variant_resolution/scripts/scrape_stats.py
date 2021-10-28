#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import itertools
from collections import defaultdict
import pysam
import svtk.utils as svu
from svtk.famfile import parse_famfile


def get_inh(called, fam):
    """
    Get list of children per inheritance status

    Returns
    -------
    inh_status : dict of {str: list of str}
        Maps inheritance status to list of corresponding samples
    """

    inh_status = dict(denovo=[], maternal=[], paternal=[], biparental=[])

    for family, samples in itertools.groupby(called, key=lambda s: fam.samples[s].family):
        samples = list(samples)

        for sample in samples:
            if fam.samples[sample].has_parents:
                if fam.samples[sample].mother in samples and fam.samples[sample].father in samples:
                    status = 'biparental'
                elif fam.samples[sample].mother in samples:
                    status = 'maternal'
                elif fam.samples[sample].father in samples:
                    status = 'paternal'
                else:
                    status = 'denovo'

                inh_status[status].append(sample)

    return inh_status


def scrape_record_stats(record, fam):
    """
    record : pysam.VariantRecord
    sample_keys : dict of {str: list of str}
        {batch: [samples]}
    """
    name = record.id
    svtype = record.info['SVTYPE']

    #  svsize = record.stop - record.pos
    svsize = record.info['SVLEN']
    algorithms = ','.join(record.info['ALGORITHMS'])

    called = svu.get_called_samples(record)
    parents = [s for s in called if not fam.samples[s].has_parents]
    children = [s for s in called if fam.samples[s].has_parents]

    n_called = len(called)
    n_parents = len(parents)
    n_children = len(children)

    homs = [s for s in called if record.samples[s]['GT'] == (1, 1)]
    hom_parents = [s for s in homs if s in parents]
    hets = [s for s in called if record.samples[s]['GT'] == (0, 1)]

    n_homs = len(homs)
    n_hets = len(hets)

    inh_status = get_inh(called, fam)
    n_denovo = len(inh_status['denovo'])
    n_maternal = len(inh_status['maternal'])
    n_paternal = len(inh_status['paternal'])
    n_biparental = len(inh_status['biparental'])

    chrom, start, end = record.chrom, record.pos, record.stop

    statline = ('{chrom}\t{start}\t{end}\t'
                '{name}\t{svtype}\t{svsize}\t{algorithms}\t'
                '{n_called}\t{n_parents}\t{n_children}\t'
                '{n_homs}\t{n_hets}\t'
                '{n_denovo}\t{n_maternal}\t{n_paternal}\t{n_biparental}')

    statline = statline.format(**locals())
    return statline


class StatsScraper:
    def __init__(self, vcf, fam, var_fout, obs_fout):
        """
        vcf : pysam.VariantRecord
        sample_keys : dict of {str: list of str}
            {batch : [samples]}
        var_fout : writable File object
            Per-variant statistics
        obs_fout : writable File object
            Per-sample statistics
        """

        self.vcf = vcf
        self.fam = fam

        self.samples = list(vcf.header.samples)

        self.var_fout = var_fout
        self.obs_fout = obs_fout

        self.var_fout.write(self.var_header)
        self.obs_fout.write(self.obs_header)

        self.sample_stats = {s: defaultdict(int) for s in self.samples}

    def scrape(self):
        for record in self.vcf:
            var_statline = scrape_record_stats(record, self.fam)
            self.var_fout.write(var_statline + '\n')

            self.scrape_sample_stats(record)

    def scrape_sample_stats(self, record):
        name = record.id
        chrom = record.chrom
        svtype = record.info['SVTYPE']
        if record.info['ALGORITHMS'] == ('depth',):
            source = 'depth-only'
        elif 'depth' in record.info['ALGORITHMS']:
            source = 'pesr+depth'
        else:
            source = 'pesr-only'

        called = svu.get_called_samples(record)
        inh_status = get_inh(called, self.fam)
        inh_map = {}
        for status, samples in inh_status.items():
            for s in samples:
                inh_map[s] = status

        fmt = ('{sample}\t{name}\t{chrom}\t{svtype}\t{source}\t{inh}\n')

        for sample in called:
            if self.fam.samples[sample].has_parents:
                inh = inh_map[sample]
            else:
                inh = 'parent'

            self.obs_fout.write(fmt.format(**locals()))

    @property
    def var_header(self):
        # TODO: add variable batches
        header = ('chrom\tstart\tend\t'
                  'name\tsvtype\tsvsize\talgorithms\t'
                  'n_called\tn_parents\tn_children\t'
                  'n_homs\tn_hets\t'
                  'n_denovo\tn_maternal\tn_paternal\tn_biparental\n')
        return header

    @property
    def obs_header(self):
        header = ('sample\t'
                  'name\tchrom\tsvtype\tsource\tinh\n')
        return header


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('famfile', type=argparse.FileType('r'))
    parser.add_argument('var_fout', type=argparse.FileType('w'))
    parser.add_argument('obs_fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    fam = parse_famfile(args.famfile)

    scraper = StatsScraper(vcf, fam, args.var_fout, args.obs_fout)
    scraper.scrape()


if __name__ == '__main__':
    main()
