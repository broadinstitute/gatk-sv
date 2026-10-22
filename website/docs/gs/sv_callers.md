---
title: SV/CNV callers
description: Summary of SV discovery tools
sidebar_position: 3
---

GATK-SV uses an ensemble of SV discovery tools to produce a raw call set and then integrates, filters, refines, 
and annotates the calls from these tools to produce a final call set.

The SV calling tools, sometimes referred to as "PE/SR" tools, include:
- [Manta](https://github.com/Illumina/manta)
- [Wham](https://github.com/zeeev/wham)
- [Scramble](https://github.com/GeneDx/scramble)
- [DRAGEN-SV](https://help.dragen.illumina.com/product-guides/dragen-v4.3/dragen-dna-pipeline/sv-calling) (not yet fully supported)

Depth-based calling of copy number variants (CNVs) is performed by two tools:
- [GATK-gCNV](https://github.com/broadinstitute/gatk)
- [cn.MOPS](https://bioconductor.org/packages/release/bioc/html/cn.mops.html)

While it is possible to omit individual tools from the pipeline, it is strongly recommended to use the default 
configuration that runs all of them.

:::note
As of 2024, most published joint call sets generated with GATK-SV used the tool MELT, a state-of-the-art mobile element 
insertion (MEI) detector, instead of Scramble. Due to licensing restrictions, we cannot provide a public docker image 
or reference panel VCFs for this algorithm. The version of the pipeline configured in the Terra workspaces does not run 
MELT or include MELT calls for the reference panel. However, the Scramble tool has replaced MELT as an MEI caller and 
should provide similar performance.
:::
