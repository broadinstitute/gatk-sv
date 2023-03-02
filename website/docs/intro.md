---
title: About
description: GATK-SV
sidebar_position: 1
---

GATK-SV is a comprehensive, cloud-based ensemble pipeline to capture and annotate all 
classes of structural variants (SV) from whole genome sequencing (WGS). It can detect 
deletions, duplications, multi-allelic copy number variants, balanced inversions, 
insertions, translocations, and a diverse spectrum of complex SV. Briefly, GATK-SV 
maximizes the sensitivity of SV discovery by harmonizing output from five tools 
(Manta, Wham, cnMOPS, GATK-gCNV, MELT). In order to reduce false positives, raw SV 
are adjudicated and re-genotyped from read evidence considering all potential 
sequencing evidence including anomalous paired-end (PE) reads, split reads (SR) crossing 
a breakpoint, normalized read-depth (RD), and B-allele frequencies (BAF). It also 
fully resolves 11 classes of complex SV based on their breakpoint signature. GATK-SV 
has been designed to be deployed in the Google Cloud Platform via the cromwell 
execution engine, which allows massively parallel scaling. Further details about 
GATK--SV can be found in [Collins et al. 2020, Nature](https://www.nature.com/articles/s41586-020-2287-8).


A high-level description of GATK-SV is available [here](https://gatk.broadinstitute.org/hc/en-us/articles/9022487952155-Structural-variant-SV-discovery).

### Citation 

Please cite the following publication:

- [Collins, Brand, et al. 2020. "A structural variation reference for medical and population genetics." Nature 581, 444-451.](https://doi.org/10.1038/s41586-020-2287-8)

Additional references: 

- [Werling et al. 2018. "An analytical framework for whole-genome sequence association studies and its implications for autism spectrum disorder." Nature genetics 50.5, 727-736.](https://doi.org/10.1038/s41588-018-0107-y)
