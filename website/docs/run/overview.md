---
title: Overview
description: Overview
sidebar_position: 1
slug: overview
---

There are two factors to consider when deciding how to run GATK-SV. 

1. **Which platform you would like to use for running GATK-SV?** 
   You may run GATK-SV on the following platforms. 
   - [Terra.bio](https://terra.bio): A user-friendly cloud-native platform for scalable data analysis. 
     The primary focus of this documentation is on supporting the execution of GATK-SV within the Terra platform.
    
   - [Cromwell](https://github.com/broadinstitute/cromwell): 
     You may run GATK-SV on a self-hosted and managed cromwell instance, which is ideal for 
     power-users and developers. We provide guidelines for this option in the 
     [_advanced guides_](/docs/advanced/cromwell) section.

2. **How many samples are you planning to process?**
    We have developed optimized flavors of GATK-SV allowing you to select the one that best suits 
    your input size for minimal cost and shorter execution time.

    | Execution Mode       | Ideal input size   | Per-sample execution cost | Per-sample execution time |
    |----------------------|--------------------|---------------------------| ------------------------- |
    | Singe-sample         | <100 samples       | $7.00-$10.00              | 24-36 hours               | 
    | Single-batch         | 100-500 samples    | $2.00-$2.50               || 
    | Multi-batch (cohort) | >200 samples       |                           ||





## GATK-SV on Terra

Please refer to the following documentation on running the GATK-SV pipeline within the Terra platform. 

- **Single-sample on Terra**: We have developed a 
  [Terra workspace](https://app.terra.bio/#workspaces/help-gatk/GATK-Structural-Variants-Single-Sample)
  for running GATK-SV on a single sample, which contains the related documentation. 
- [Single-batch on Terra](batch.md);
- [Multiple-batch (cohort) on Terra](cohort.md).

## GATK-SV on Cromwell

- [Single-sample on Cromwell](/docs/advanced/cromwell/single);
- [Single-batch on Cromwell](/docs/advanced/cromwell/batch);
- [Multiple-batch (cohort) on Cromwell](/docs/advanced/cromwell/cohort).
