---
title: Calling modes
description: Description of single-sample and joint calling
sidebar_position: 1
---

# Calling Modes

GATK-SV offers three different modes for SV calling. 
This section guides you to determine
which mode is appropriate for your use case.

Please refer to the table below for a quick comparison; 
each mode is covered in more detail later in this page.

|  | [Cohort Mode](#join-calling) | [Single-sample (Terra)](#single-sample-terra) | [Single-sample (SV-Shell)](#single-sample-sv-shell) |
|--|------------------------------|-----------------------------------------------|-----------------------------------------------------|
| **Designed for**  | Biobank-level joint-calling SVs                               | Clinical use-cases                                      | Clinical use-cases                                      |
| **Use when**      | All the data is available when starting to run the pipeline   | Data is available in smaller batches on a rolling basis | Data is available in smaller batches on a rolling basis |
| **Scalability**   | Hundreds of thousands of WGS                                  | One WGS at a time                                       | One WGS at a time                                       |
| **Availability**  | Exclusively Terra platform running on GCP                     | Exclusively Terra platform running on GCP               | Platform-agnostic and available on Terra                |


## Single-Sample Mode (GCP Only) {#single-sample-terra}

GATK-SV can perform SV calling on individual samples. In this mode, a sample is jointly called against a fixed reference 
panel of [156 high-quality samples from the 1000 Genomes Project](https://app.terra.bio/#workspaces/anvil-datastorage/1000G-high-coverage-2019). 
Single-sample mode is a good option for the following 
use cases:

- Studies involving fewer 100 samples
- Studies with rolling data delivery, i.e. in small batches over time

Users should also consider that the single-sample mode is provided as a single workflow and is therefore considerably 
simpler to run than joint calling. However, it also has higher compute costs on a per-sample basis and will not be as sensitive 
as joint calling with larger cohorts. Additionally, SV quality will be best when the case sample closely resembles the samples
in the reference panel in terms of sequencing depth, sample quality, and library preparation.

Please refer to [this page](/docs/execution/single) for more details on this mode.

## Single-Sample Mode (Platform-Agnostic) {#single-sample-sv-shell}

Single-sample mode runs exclusively on the Terra platform, which supports only GCP.
To address this limitation, we developed SV-Shell, 
a platform-agnostic reimplementation of the single-sample mode.
It offers the same functionality as single-sample mode, 
and additionally supports optional DRAGEN-generated SV and CNV calls.

You may want to use SV-Shell for the following use cases:

- Running on platforms other than GCP/Terra, such as AWS, Azure, or institutional HPC
- Studies involving fewer than 100 samples
- Studies with rolling data delivery, i.e. in small batches over time
- Using DRAGEN-generated SV and CNV calls

Please refer to [this page](/docs/execution/svshell) for more details on this mode.

## Joint Calling Mode {#join-calling}

GATK-SV can also perform joint calling on a set of samples. Users may opt for this mode in the following use cases:

- Studies involving at least 100 samples
- When maximum sensitivity is desired
- Data sets that are technically heterogeneous, i.e. with strong batch effects, or are very different from the single-sample mode reference panel

Joint calling has the advantage of increasing SV discovery sensitivity and providing allele frequency estimates, and there are 
some features, such as genotype recalibration and filtering and in-depth QC plotting, that are only available in joint calling mode. 
However, this pipeline is considerably more complex to execute than the single-sample mode, requiring sample batching and the execution of 
several individual modules.

Please refer to [this page](/docs/execution/joint) for more details on this mode.
