---
title: Overview
description: Overview
sidebar_position: 1
slug: overview
---

GATK-SV uses [Manta](https://github.com/Illumina/manta), [WHAM](https://github.com/zeeev/wham), 
[Scramble](https://github.com/GeneDx/scramble), [GATK gCNV](https://github.com/broadinstitute/gatk), 
and [cn.MOPS](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)
as raw calling algorithms, and then integrates, filters, refines, and annotates
the calls from these tools to produce a final call set. 

:::info
Please note that most large published joint call sets produced by
GATK-SV, including gnomAD-SV, included the tool MELT, a state-of-the-art mobile element insertion (MEI) detector, as 
part of the pipeline. Due to licensing restrictions, we cannot provide a public docker image or reference panel VCFs 
for this algorithm. The version of the pipeline configured in this workspace does not run MELT or include MELT calls 
for the reference panel. However, the Scramble tool has replaced MELT as an MEI caller and should provide similar 
performance.
:::

This section provides technical documentation for the running single-sample and joint calling modes in Terra.

Please make sure to review the [Getting Started](/docs/gs/overview) material before proceeding to run the pipeline.
