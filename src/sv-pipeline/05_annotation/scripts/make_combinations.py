#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pandas as pd


def make_combinations(genic, genelists, noncoding, functional):
    # Assign all variants ANY effect and ANY genelist
    any_variants = genic.copy()
    any_variants['effect'] = 'Any'
    any_variants['gene'] = 'Any'
    any_variants = any_variants.rename(columns=dict(gene='genelist'))
    any_variants = any_variants.drop_duplicates()

    # Make combinations of ANY effect with specific gene list
    any_effects = genelists.copy()
    any_effects['effect'] = 'Genic'
    any_effects = any_effects.drop_duplicates()

    # Make combinations of specific effect with ANY gene list
    any_genelists = genelists.copy()
    any_genelists['genelist'] = 'Any'
    any_genelists = any_genelists.drop_duplicates()

    any_effect_genelist = genelists.copy()
    any_effect_genelist['effect'] = 'Genic'
    any_effect_genelist['genelist'] = 'Any'
    any_effect_genelist = any_effect_genelist.drop_duplicates()

    # Get all combos of genelists (restricted to genic variation)
    genelists = pd.concat([genelists, any_effects, any_genelists,
                           any_effect_genelist])
    genelists = genelists.drop_duplicates()

    # Combine noncoding and intergenic variation
    # Only requires noncoding effect plus ANY gene list
    noncoding = noncoding.rename(columns=dict(gene='genelist'))
    noncoding['genelist'] = 'Any'

    intergenic = genic.loc[(genic.effect == 'INTERGENIC') &
                           (~genic.name.isin(noncoding.name))].copy()
    intergenic = intergenic.rename(columns=dict(gene='genelist'))
    intergenic['genelist'] = 'Any'
    noncoding = pd.concat([noncoding, intergenic])

    # Merge all combos of effect and genelist
    genic = pd.concat([any_variants, genelists, noncoding]).drop_duplicates()

    # Any functional variation
    any_functional = genic['name effect'.split()].copy()
    any_functional = any_functional.rename(columns={'effect': 'functional'})
    any_functional['functional'] = 'Any'
    functional = functional.rename(columns={'NONCODING': 'functional'})
    functional = pd.concat([functional, any_functional]).drop_duplicates()

    # Take outer join to make all combinations of effect, genelist, and
    # functional element
    combos = pd.merge(genic, functional, on='name', how='outer')

    combos['combo'] = (combos.effect + '_' +
                       combos.genelist + '_' +
                       combos.functional)

    return combos


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('genic')
    parser.add_argument('genelists')
    parser.add_argument('noncoding')
    parser.add_argument('functional')
    parser.add_argument('fout')
    args = parser.parse_args()

    genic = pd.read_table(args.genic)
    genelists = pd.read_table(args.genelists)
    noncoding = pd.read_table(args.noncoding)
    functional = pd.read_table(args.functional)

    combos = make_combinations(genic, genelists, noncoding, functional)

    combos.to_csv(args.fout, index=False, sep='\t')


if __name__ == '__main__':
    main()
