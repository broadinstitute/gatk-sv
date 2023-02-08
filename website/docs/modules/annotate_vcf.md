---
title: AnnotateVCF
description: Annotate VCF (work in progress)
sidebar_position: 13
slug: av
---

Add annotations, such as the inferred function and 
allele frequencies of variants, to final VCF.

Annotations methods include:

- Functional annotation - annotate SVs with inferred functional 
  consequence on protein-coding regions, regulatory regions 
  such as UTR and promoters, and other non-coding elements.

- Allele Frequency annotation - annotate SVs with their allele 
  frequencies across all samples, and samples of specific sex, 
  as well as specific sub-populations.

- Allele Frequency annotation with external callset - annotate 
  SVs with the allele frequencies of their overlapping SVs in 
  another callset, eg. gnomad SV callset.
