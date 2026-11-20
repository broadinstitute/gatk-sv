---
title: PlotCnvReadDepth
description: Plot CNV read depth
sidebar_position: 25
slug: viz
---

Generates plots of read depth across all samples. This is useful for visualizing large deletions and duplications 
(at least 5kbp).

This module is not a part of the core pipeline but may be used to investigate variants of interest.

### Inputs

:::info
All array inputs of batch data must match in order. In particular, the ordering of `batch_name_list`, 
`batch_sample_lists`, `medians_files`, `bincov_files`, and `bincov_index_files` must be the same.
:::

#### `vcfs`
Array of VCFs containing variants to plot. Variants will be automatically subsetted to `PASS` DEL and DUP types subject 
to the [min_size](#min_size) constraint. Ignored if [cnvs](#cnvs) is provided.

#### `cnvs`
Optional. Array of tab-delimited files of CNVs to plot, used instead of [vcfs](#vcfs). Files must not have a header and 
must contain columns: `chrom,start,end,name,svtype,samples`, where `samples` is comma-delimited.

#### `variant_ids`
Optional. File of variant IDs to plot from [vcfs](#vcfs), one per line. Variants will still be restricted to DEL and 
DUP types.

#### `plot_prefix`
Output prefix, such as cohort name. May be alphanumeric with underscores.

#### `min_size`
Minimum size in bases of variants to plot.

#### `one_sample_per_plot`
Optional. If `true`, generates a separate plot for each carrier sample. Default: `false`.

#### `padding`
Optional. Padding around each CNV when plotting, as a fraction of CNV length. Default: `0.5`.

#### `variants_per_shard`
Optional. Number of variants to plot per shard. Default: `40`.

#### `batch_name_list`
Array of batch IDs for all batches in the input variants.

#### `batch_sample_lists`
Array of sample ID lists, one per batch, generated in [FilterBatch](./fb#filtered_batch_samples_file).

#### `medians_files`
Array of median coverage files for all batches in the input variants, generated in [GatherBatchEvidence](./gbe#median_cov).

#### `bincov_files`
Array of RD evidence files for all batches in the input variants, generated in [GatherBatchEvidence](./gbe#merged_bincov), 
and their indexes `bincov_index_files`.

### Outputs

#### `cnv_plots`
Tarball containing output plots.
