---
title: Overview
description: Overview
sidebar_position: 1
slug: overview
---

There are two factors to consider when deciding how to run GATK-SV. 


1. **Variant calling modes: single-sample and cohort-based calling.**
  GATK-SV offers two distinct pipeline configurations for detecting 
  structural variations (SVs), each tailored for different research needs:

    - **Single-sample analysis:** 
      This configuration is ideal for examining SVs in individual samples, 
      focusing exclusively on data from that single sample. Running this mode is less complex, 
      involving just one workflow per sample.

    - **Joint calling:** 
      This configuration is designed for more extensive studies, such as those 
      involving population genetics or disease association studies. 
      It analyzes SVs across a cohort by collectively assessing data from all samples. 
      However, this comes with increased complexity compared to the single-sample mode, 
      requiring the execution of multiple workflows and involves data preparation steps 
      (e.g., batching files from the cohort).
   

2. **Which platform you would like to use for running GATK-SV?** 
   You may run GATK-SV on the following platforms. 

   - [Terra.bio](https://terra.bio): A user-friendly cloud-native platform for scalable data analysis. 
     The primary focus of this documentation is on supporting the execution of GATK-SV within the Terra platform.
   
   - [Cromwell](https://github.com/broadinstitute/cromwell): 
     You may run GATK-SV on a self-hosted and managed cromwell instance, which is ideal for 
     power-users and developers. We provide guidelines for this option in the 
     [_advanced guides_](/docs/advanced/development/cromwell) section.


Your decision regarding the execution modes and platform should be guided by 
the objectives of your study, the size of your cohort, data access needs, 
and the trade-off between a straightforward interface (Terra) 
and more detailed customization options (self-managed Cromwell server). 
Please refer to the following documentation on running GATK-SV within the Terra platform.

- [Single-sample on Terra](single.md);
- [Joint calling on Terra](joint.md).
