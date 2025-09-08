#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""


import argparse
import numpy as np
import pysam
import svtk.utils as svu
import pandas as pd
from svtk.pesr import PESRTestRunner, SRTest, PETest
from svtk.famfile import parse_famfile


def get_denovo_candidates(record, fam, max_parents=20):
    """
    Obtain list of samples which are putatively called de novo
    """
    called = svu.get_called_samples(record)

    parents = [s for s in called if fam.samples[s].is_parent]

    if len(parents) > max_parents:
        return []

    denovo = []

    for ID in called:
        sample = fam.samples[ID]
        if sample.has_parents:
            if sample.mother not in called and sample.father not in called:
                denovo.append(sample.ID)

    return denovo


class DenovoTestRunner(PESRTestRunner):
    def __init__(self, vcf, fam, countfile, discfile, pe_fout, sr_fout,
                 n_background=160, max_parents=10):

        super().__init__(vcf, n_background)

        self.srtest = SRTest(countfile)
        self.sr_fout = sr_fout

        self.petest = PETest(discfile)
        self.pe_fout = pe_fout

        self.fam = fam

        self.max_parents = max_parents

    def run(self):
        for record in self.vcf:
            # Skip records without any Mendelian violations
            candidates = get_denovo_candidates(
                record, self.fam, self.max_parents)
            if len(candidates) == 0:
                continue

            # Restrict to rare (parental VF<0.1%) variants
            called = svu.get_called_samples(record)
            parents = [s for s in called if self.fam.samples[s].is_parent]
            if len(parents) > self.max_parents:
                continue

            # Skip non-stranded (wham)
            if record.info['STRANDS'] not in '+- -+ ++ --'.split():
                continue

            self.test_record(record, candidates)

    def test_record(self, record, candidates):
        for child in candidates:
            called = [child]
            background = self.choose_background(record, candidates)

            self.sr_test(record, called, background)
            self.pe_test(record, called, background)

    def sr_test(self, record, called, background):
        sr_results = self.srtest.test_record(record, called, background)

        child = called[0]
        sr_results['sample'] = child

        # Check parents at predicted coordinates in child
        posA, posB = sr_results.set_index('coord')['pos'][['posA', 'posB']]
        parents = self.sr_test_parents(record, child, background,
                                       posA, posB)

        # Merge parents and children
        results = pd.concat([sr_results, parents], ignore_index=True)

        # sneak posB-posA distance into output
        dist = posB - posA
        results.loc[results.coord == 'sum', 'pos'] = dist

        cols = 'name sample coord pos log_pval called background'.split()
        results = results[cols]
        results['pos'] = results.pos.astype(int)
        results = results.rename(columns={'called': 'called_median',
                                          'background': 'bg_median'})

        results.to_csv(self.sr_fout, index=False, header=False,
                       sep='\t', na_rep='NA', float_format='%.20f')

    def pe_test(self, record, called, background):
        pe_results = self.petest.test_record(record, called, background)

        child = called[0]
        pe_results['sample'] = child

        p_results = []
        for parent in [self.fam.samples[child].father, self.fam.samples[child].mother]:
            result = self.petest.test_record(record, [parent], background)
            result['sample'] = parent
            p_results.append(result)

        p_results = pd.concat(p_results)
        results = pd.concat([pe_results, p_results], ignore_index=True)

        cols = 'name sample log_pval called background'.split()
        results[cols].to_csv(self.pe_fout, index=False, header=False,
                             sep='\t', na_rep='NA', float_format='%.20f')

    def choose_background(self, record, candidates):
        # Exclude called samples and all candidate families from background
        related = [c for c in candidates]
        for s in candidates:
            if self.fam.samples[s].has_parents:
                related.append(self.fam.samples[s].father)
                related.append(self.fam.samples[s].mother)

                related += self.fam.samples[self.fam.samples[s].father].children

        called = set(svu.get_called_samples(record))
        blacklist = called.union(related)

        background = [s for s in self.samples if s not in blacklist]

        if len(background) >= self.n_background:
            background = np.random.choice(background, self.n_background,
                                          replace=False).tolist()

        return background

    def sr_test_parents(self, record, child, background, posA, posB):
        results = []
        for parent in [self.fam.samples[child].father, self.fam.samples[child].mother]:
            p_results = []
            for pos, strand in zip([posA, posB], record.info['STRANDS']):
                result = self.srtest.test(record.chrom, pos, strand,
                                          [parent], background)
                result = result.to_frame().transpose()
                result['coord'] = 'posA' if pos == posA else 'posB'
                result['pos'] = pos

                p_results.append(result)

            p_results = pd.concat(p_results, ignore_index=True)
            p_total = self.srtest._test_total(p_results)

            p_results = pd.concat([p_results, p_total], ignore_index=True)
            p_results['sample'] = parent
            results.append(p_results)

        results = pd.concat(results, ignore_index=True)
        results['name'] = record.id

        return results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('famfile', type=argparse.FileType('r'))
    parser.add_argument('-c', '--countfile', required=True)
    parser.add_argument('-d', '--discfile', required=True)
    parser.add_argument('--discfile-index')
    parser.add_argument('--countfile-index')
    parser.add_argument('--background', type=int, default=160)
    parser.add_argument('--max-parents', type=float, default=10)
    parser.add_argument('petest', type=argparse.FileType('w'), help='fout')
    parser.add_argument('srtest', type=argparse.FileType('w'), help='fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    fam = parse_famfile(args.famfile)

    if args.discfile_index is None:
        discfile = pysam.TabixFile(args.discfile)
    else:
        discfile = pysam.TabixFile(args.discfile, index=args.discfile_index)

    if args.countfile_index is None:
        countfile = pysam.TabixFile(args.countfile)
    else:
        countfile = pysam.TabixFile(args.countfile, index=args.countfile_index)

    header = 'name sample log_pval called_median bg_median'.split()
    args.petest.write('\t'.join(header) + '\n')

    header = 'name sample coord pos log_pval called_median bg_median'.split()
    args.srtest.write('\t'.join(header) + '\n')

    runner = DenovoTestRunner(vcf, fam, countfile, discfile,
                              args.petest, args.srtest,
                              args.background, args.max_parents)
    runner.run()


if __name__ == '__main__':
    main()
