---
title: About
description: GATK-SV
sidebar_position: 1
---

GATK-SV is a comprehensive, cloud-based ensemble pipeline for discovering and annotating all 
classes of structural variants (SV) from short-read whole genome sequencing (WGS) data. It can detect 
deletions, duplications, multi-allelic copy number variants, balanced inversions, 
insertions, translocations, and a diverse spectrum of complex SV. Briefly, GATK-SV 
maximizes the sensitivity of SV discovery by harmonizing output from five tools: 
Manta, Wham, Scramble, cn.MOPS, and GATK-gCNV. To minimize false positives, raw SVs 
are adjudicated and re-genotyped, considering all potential 
sequencing evidence including anomalous paired-end (PE) reads, split reads (SR), 
read-depth (RD), and B-allele frequencies (BAF). It also fully resolves 11 classes of complex 
SVs composed of multiple breakpoints. GATK-SV is intended for use on the [Terra](https://app.terra.bio/) 
 platform.

### Methods

Further details about GATK-SV methods can be found in [Collins et al. 2020](https://www.nature.com/articles/s41586-020-2287-8).

### GATK Best Practices

Additional guidance on running GATK-SV is also available [here](https://gatk.broadinstitute.org/hc/en-us/articles/9022653744283-GATK-Best-Practices-for-Structural-Variation-Discovery-on-Single-Samples).

### Where to go from here

This documentation includes instructions for running the pipeline, technical implementation details, troubleshooting 
information, and guides for advanced users who wish to work with the source code or rebuild the project.

We recommend new users continue to the [Getting Started overview](/docs/gs/overview.md).
