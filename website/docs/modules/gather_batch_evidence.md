---
title: GatherBatchEvidence
description: Gather Batch Evidence
sidebar_position: 4
slug: gbe
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/GatherBatchEvidence.wdl)

Runs CNV callers ([cn.MOPS](https://academic.oup.com/nar/article/40/9/e69/1136601), [GATK-gCNV](https://github.com/broadinstitute/gatk)) 
and combines single-sample raw evidence into batched files.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gbe: GatherBatchEvidence
  t: TrainGCNV
  cb: ClusterBatch
  t --> gbe
  gbe --> cb
  
  class gbe thisModule
  class t inModules
  class cb outModules
```

## Inputs

:::info
All array inputs of sample data must match in order. For example, the order of the `samples` array should match that of
`counts`, `PE_files`, etc.
:::

#### `batch`
An identifier for the batch; may only be alphanumeric with underscores.

#### `samples`
Sample IDs. Must match the sample IDs used in [GatherSampleEvidence](./gse#sample_id) unless `rename_samples` is enabled, in 
which case sample IDs will be overwritten. See [sample ID requirements](/docs/gs/inputs#sampleids) for specifications 
of allowable sample IDs.

#### `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format).

#### `counts`
Binned read count files (`*.rd.txt.gz`) generated in [GatherSampleEvidence](./gse#coverage-counts).

#### `PE_files`
Discordant pair evidence files (`*.pe.txt.gz`) generated in [GatherSampleEvidence](./gse#pesr-disc).

#### `SR_files`
Split read evidence files (`*.sr.txt.gz`) generated in [GatherSampleEvidence](./gse#pesr-split).

#### `SD_files`
Site depth files (`*.sd.txt.gz`) generated in [GatherSampleEvidence](./gse#pesr-sd).

#### `*_vcfs`
Raw caller VCFs generated in [GatherSampleEvidence](./gse#outputs). Callers may be omitted if they were not run.

#### `run_matrix_qc`
Enables running QC tasks.

#### `contig_ploidy_model_tar`
Contig ploidy model tarball generated in [TrainGCNV](./gcnv#cohort_contig_ploidy_model_tar).

#### `gcnv_model_tars`
CNV model tarball generated in [TrainGCNV](./gcnv#cohort_gcnv_model_tars).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `rename_samples`
Default: `false`. Overwrite sample IDs with the [samples](#samples) input.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `run_ploidy`
Default: `false`. Runs ploidy estimation. Note this calls the same method used in [EvidenceQc](./eqc).

## Outputs

#### `merged_BAF`
Batch B-allele frequencies file (`.baf.txt.gz`) derived from site depth evidence.

#### `merged_SR`
Batch split read evidence file (`.sr.txt.gz`).

#### `merged_PE`
Batch paired-end evidence file (`.pe.txt.gz`).

#### `merged_bincov`
Batch binned read counts file (`.rd.txt.gz`).

#### `merged_dels`, `merged_dups`
Batch CNV calls (`.bed.gz`).

#### `median_cov`
Median coverage table.

#### `std_*_vcf_tar`
Tarballs containing per-sample raw caller VCFs in standardized formats. This will be ommitted for any callers not 
provided in the inputs.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `batch_ploidy_*`
Ploidy analysis files. Enabled with [run_ploidy](#optional-run_ploidy).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `*_stats`, `Matrix_QC_plot`
QC files. Enabled with [run_matrix_qc](#run_matrix_qc).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `manta_tloc`
Supplemental evidence for translocation variants. These records are hard filtered from the main call set but may be of 
interest to users investigating reciprocal translocations and other complex events.



