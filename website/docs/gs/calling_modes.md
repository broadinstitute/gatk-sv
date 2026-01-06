---
title: Calling modes
description: Description of single-sample and joint calling
sidebar_position: 1
---

# Calling modes

GATK-SV offers two different modes for SV calling. Users should carefully review the following sections to determine
which mode is appropriate for their use case.

## Single-sample mode

GATK-SV can perform SV calling on individual samples. In this mode, a sample is jointly called against a fixed reference 
panel of [156 high-quality samples from the 1000 Genomes Project](https://app.terra.bio/#workspaces/anvil-datastorage/1000G-high-coverage-2019). Single-sample mode is a good option for the following 
use cases:

- Studies involving fewer 100 samples
- Studies with rolling data delivery, i.e. in small batches over time

Users should also consider that the single-sample mode is provided as a single workflow and is therefore considerably 
simpler to run than joint calling. However, it also has higher compute costs on a per-sample basis and will not be as sensitive 
as joint calling with larger cohorts. Additionally, SV quality will be best when the case sample closely resembles the samples
in the reference panel in terms of sequencing depth, sample quality, and library preparation.

## Joint calling mode

GATK-SV can also perform joint calling on a set of samples. Users may opt for this mode in the following use cases:

- Studies involving at least 100 samples
- When maximum sensitivity is desired
- Data sets that are technically heterogeneous, i.e. with strong batch effects, or are very different from the single-sample mode reference panel

Joint calling has the advantage of increasing SV discovery sensitivity and providing allele frequency estimates, and there are 
some features, such as genotype recalibration and filtering and in-depth QC plotting, that are only available in joint calling mode. 
However, this pipeline is considerably more complex to execute than the single-sample mode, requiring sample batching and the execution of 
several individual modules.

## Related content

More information on single-sample and joint calling can be found in the [Execution](/docs/execution/overview) section.

