---
title: TrainGCNV
description: Train gCNV
sidebar_position: 3
slug: gcnv
---

import { Highlight, HighlightOptionalArg } from "../../src/components/highlight.js"

[GATK-gCNV](https://www.nature.com/articles/s41588-023-01449-0)
is a method for detecting rare germline copy number variants (CNVs)
from short-read sequencing read-depth information.
The [TrainGCNV](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/TrainGCNV.wdl)
module trains a gCNV model for use in the [GatherBatchEvidence](./gbe) workflow. 
The upstream and downstream dependencies of the TrainGCNV module are illustrated in the following diagram.


The samples used for training should be homogeneous and similar 
to the samples on which the model will be applied in terms of sample type, 
library preparation protocol, sequencer, sequencing center, and etc.


For small, relatively homogeneous cohorts, a single gCNV model is usually sufficient. 
However, for larger cohorts, especially those with multiple data sources, 
we recommend training a separate model for each batch or group of batches (see 
[batching section](/docs/run/joint#batching) for details).
The model can be trained on all or a subset of the samples to which it will be applied. 
A subset of 100 randomly selected samples from the batch is a reasonable
input size for training the model; when the `n_samples_subsample` input is provided, 
the `TrainGCNV` workflow can automatically perform this random selection.

The following diagram illustrates the upstream and downstream workflows of the `TrainGCNV` workflow 
in the recommended invocation order. You may refer to 
[this diagram](https://github.com/broadinstitute/gatk-sv/blob/main/terra_pipeline_diagram.jpg) 
for the overall recommended invocation order.

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#00509d,color:#caf0f8
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  batching: Batching, sample QC, and sex assignment
  t: TrainGCNV
  gbe: GatherBatchEvidence
  
  batching --> t
  t --> gbe 
  
  class t thisModule
  class batching inModules
  class gbe outModules
```

## Inputs

This section provides a brief description on the _required_ inputs of the TrainGCNV workflow.
For a description on the _optional_ inputs and their default values, you may refer to the 
[source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/TrainGCNV.wdl) of the TrainGCNV workflow.
Additionally, the majority of the optional inputs of the workflow map to the optional arguments of the 
tool the workflow uses, `GATK GermlineCNVCaller`; hence, you may refer to the 
[documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360040097712-GermlineCNVCaller) 
of the tool for a description on these optional inputs. 

#### `samples`
A list of sample IDs. 
The order of IDs in this list should match the order of files in `count_files`.

#### `count_files`
A list of per-sample coverage counts generated in the [GatherSampleEvidence](./gse#outputs) workflow.

#### `contig_ploidy_priors`
A tabular file with ploidy prior probability per contig.
You may find the link to this input from 
[this reference](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json)
and a description to the file format 
[here](https://gatk.broadinstitute.org/hc/en-us/articles/360037224772-DetermineGermlineContigPloidy).


#### `reference_fasta`
`reference_fasta`, `reference_index`, `reference_dict` are respectively the 
reference genome sequence in the FASTA format, its index file, and a corresponding 
[dictionary file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format).
You may find links to these files from 
[this reference](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json).


## Outputs

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `annotated_intervals` {#annotated-intervals}

The count files from [GatherSampleEvidence](./gse) with adjacent intervals combined into 
locus-sorted `DepthEvidence` files using `GATK CondenseDepthEvidence` tool, which are
annotated with GC content, mappability, and segmental-duplication content using 
[`GATK AnnotateIntervals`](https://gatk.broadinstitute.org/hc/en-us/articles/360041416652-AnnotateIntervals)
tool. This output is generated if the optional input `do_explicit_gc_correction` is set to `True`.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `filtered_intervals_cnv` {#filtered-intervals-cnv}

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `cohort_contig_ploidy_model_tar` {#contig-ploidy-model-tarball}

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `cohort_gcnv_model_tars` {#model-tarballs}
