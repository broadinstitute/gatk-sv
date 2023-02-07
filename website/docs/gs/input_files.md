---
title: Input Data
description: Supported input and reference data.
sidebar_position: 3
slug: ./inputs
---

GATK-SV supports input data in the following formats.

- Illumina short-read whole-genome CRAMs or BAMs, aligned to hg38 with [bwa-mem](https://github.com/lh3/bwa). 
  BAMs must also be indexed.

- Family structure definitions file in 
  [PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).  
  
  :::info
  Sex aneuploidies (detected in EvidenceQC) should be entered as sex = 0.
  :::
