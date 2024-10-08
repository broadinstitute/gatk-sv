---
title: GatherSampleEvidence 
description: Gather Sample Evidence
sidebar_position: 1
slug: gse
---

Runs raw evidence collection on each sample with the following SV callers: 
Manta, Wham, Scramble, and/or MELT, and collects raw SV evidence (PE/SR/RD/SD). For guidance on pre-filtering prior 
to GatherSampleEvidence, refer to the Sample Exclusion section.

The following diagram illustrates the downstream workflows of the `GatherSampleEvidence` workflow 
in the recommended invocation order. You may refer to 
[this diagram](https://github.com/broadinstitute/gatk-sv/blob/main/terra_pipeline_diagram.jpg) 
for the overall recommended invocation order.


```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#00509d,color:#caf0f8
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gse: GatherSampleEvidence
  eqc: EvidenceQC
  gse --> eqc
  
  class gse thisModule
  class eqc outModules
```


### Required inputs

#### `bam_or_cram_file`
An indexed BAM or CRAM file aligned to hg38. See [Input data requirements](/docs/gs/inputs).

#### `sample_id`
Refer to the [sample ID requirements](/docs/gs/inputs#sampleids) for specifications of allowable sample IDs. 
IDs that do not meet these requirements may lead to errors.

### Optional inputs

#### `collect_coverage`
Default: `true`. Collect read depth.

#### `collect_pesr`
Default: `true`. Collect paired-end and split-read evidence.

#### `manta_docker`
If provided, runs the Manta tool.

#### `melt_docker`
If provided, runs the MELT tool.

#### `scramble_docker`
If provided, runs the Scramble tool.

#### `wham_docker`
If provided, runs the Wham tool.

#### `run_localize_reads`
Default: `false`. Copy input alignment files to the execution bucket before localizing to subsequent tasks. This 
may be desirable when BAM/CRAM files are stored in a requester-pays bucket or in another region to avoid egress charges.

:::warning
Enabling `run_localize_reads` can incur high storage costs. If using, make sure to clean up execution directories from 
failed runs.
:::

#### `run_module_metrics`
Default: `true`. Calculate QC metrics for the sample. If true, `primary_contigs_fai` must also be provided, and 
optionally the `baseline_*_vcf` inputs to run comparisons. 

#### `run_localize_reads`
Default: `false`. Copy input alignment files to the execution bucket before localizing to subsequent tasks.

### Outputs

#### `manta_vcf` {#manta-vcf}
VCF containing variants called by Manta

#### `melt_vcf` {#melt-vcf}
VCF containing variants called by MELT

#### `scramble_vcf` {#scramble-vcf}
VCF containing variants called by Scramble

#### `wham_vcf` {#wham-vcf}
VCF containing variants called by Wham

#### `coverage_counts` {#coverage-counts}
Binned read counts collected by `GATK-CollectReadCounts` (`.counts.tsv.gz`)

#### `pesr_disc` {#pesr-disc}
Discordant read pairs collected by `GATK-CollectSVEvidence` (`.pe.txt.gz`)

#### `pesr_split` {#pesr-split}
Split read positions collected by `GATK-CollectSVEvidence` (`.sr.txt.gz`)

#### `pesr_sd` {#pesr-sd}
Site depth counts collected by `GATK-CollectSVEvidence` (`.sd.txt.gz`)

#### `sample_metrics_files`
(Optional) Sample metrics for QC. Enable by setting `run_module_metrics = true`.