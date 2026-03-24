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
public configuration, which runs Manta, Wham, Scramble, GATK-gCNV, and cn.MOPS.

:::note
Historically, many published GATK-SV call sets used MELT, a state-of-the-art mobile element 
insertion (MEI) detector, instead of Scramble. Due to licensing restrictions, we cannot provide the public Docker image, 
reference-panel assets, or Terra workspace defaults needed to run MELT out of the box. The repository still contains 
optional MELT plumbing for advanced users who can supply their own MELT resources, but the public Terra workspaces and 
distributed reference panel use Scramble as the supported MEI caller.
:::
