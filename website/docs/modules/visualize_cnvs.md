---
title: VisualizeCnvs
description: Visualize CNVs
sidebar_position: 23
slug: viz
---

Generates plots of read depth across all samples. This is useful for visualizing large deletions and duplications 
(at least 5kbp).

This module is not a part of the core pipeline but may be used to investigate variants of interest.

### Inputs

:::info
All array inputs of batch data must match in order. In particular, the ordering of `median_files` and `rd_files` must 
be the same.
:::

#### `vcf_or_bed`
VCF or bed file containing variants to plot. All variants will be automatically subsetted to DEL and DUP types subject 
to the [min_size](#min_size) constraint. VCF files must end in `.vcf.gz` and bed files must end in either `.bed` or 
`.bed.gz`. Bed files must contain columns: `chrom,start,end,name,svtype,samples`.

#### `prefix`
Output prefix, such as cohort name. May be alphanumeric with underscores.

#### `median_files`
Array of median coverage files for all batches in the input variants, generated in [GatherBatchEvidence](./gbe#median_cov).

#### `rd_files`
Array of RD evidence files for all batches in the input variants, generated in [GatherBatchEvidence](./gbe#merged_bincov).

#### `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format).

#### `min_size`
Minimum size in bases of variants to plot.

#### `flags`
Additional flags to pass to the [RdTest plotting script](https://github.com/broadinstitute/gatk-sv/blob/main/src/RdTest/RdTest.R).

:::warning
Due to a bug, the `flags` parameter must contain `-s 999999999` in order to properly plot variants over 1 Mb. 
:::

### Outputs

#### `rdtest_plots`
Tarball containing output plots.

