---
title: Downstream Filtering
description: Downstream filtering (work in progress)
sidebar_position: 12
slug: df
---

Apply downstream filtering steps to the cleaned VCF to further 
control the false discovery rate; all steps are optional 
and users should decide based on the specific purpose of 
their projects.

Filtering methods include:

- minGQ - remove variants based on the genotype quality across 
  populations. Note: Trio families are required to build the minGQ 
  filtering model in this step. We provide tables pre-trained with 
  the 1000 genomes samples at different FDR thresholds for projects 
  that lack family structures, and they can be found at the paths below. 
  These tables assume that GQ has a scale of [0,999], so they will 
  not work with newer VCFs where GQ has a scale of [0,99].

  ```shell
  gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.10perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
  gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.1perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
  gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.5perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
  ```
  
- BatchEffect - remove variants that show significant discrepancies 
  in allele frequencies across batches

- FilterOutlierSamplesPostMinGQ - remove outlier samples with unusually 
  high or low number of SVs

- FilterCleanupQualRecalibration - sanitize filter columns and 
  recalibrate variant QUAL scores for easier interpretation
