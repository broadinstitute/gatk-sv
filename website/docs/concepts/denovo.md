---
title: de-novo SVs
description: Downstream filtering (work in progress)
sidebar_position: 1
slug: denovo
---


The de-novo workflow identifies a variant as de-novo based the criteria defined 
on site and sample specific metrics of a variant.
These criteria are discussed in details in the following sections. 


## Site-specific filters

- Exclude multi-allelic copy number variants (mCNV) and breakends (BND).
- Exclude variants with gnomAD `AF > 0.01`, cohort `AF > 0.05` or parents `AF > 0.05` (except for GDs).
- Remove Wham-only calls with `GQ = 1`.

## Sample-specific filters

- De novo SV (call absent in parents)
- CNVs filters:
  - Large CNVs (>5,000bp) in autosomal chromosomes: 
    - Remove if read depth copy number (RD_CN) of parents is equal to RD_CN of proband.
    - Remove if the variant overlaps with a CNV in the parents, with minimum overlap of 0.5 (using bedtools coverage).
    - Check that variant has raw depth support (uses bedtools coverage for >5000bp, and intersect for 1000-5000bp).
    - Remove if variant overlaps with parental raw depth CNV (uses bedtools coverage for >5000bp, and intersect for 1000-5000bp).
  - Small CNVs (<=5000bp):
    - If the CNV is small, ignore RD evidence if “RD,SR” (and use only SR)
    - Check that variant has raw evidence (bedtools intersect)
    - Remove if variant overlaps with parental raw evidence (bedtools intersect)
    - Remove if variant is SR only but doesn’t have “bothsides support”
  - Remove depth only calls that are <10000bp
  - Remove deletions that are >500bp and RD_CN=2 and PE only evidence
  - Flag if parents GQ is smaller than min_gq (does not apply, we have been using 0)
  - Flag INS that are either manta or melt and SR only and have high SR background (filter removed)
  - Flag if median coverage in parents is <= 10
- Blacklist: remove variants that overlap >0.5 with “blacklist” (repetitive and bad regions)
- INS filters:
  - Requires raw evidence
  - Remove variant if it’s in parental raw evidence file
- No filtering is applied on complex SVs, translocations and inversions
- Reformatting is applied at the end to remove duplicated CPX SVs.