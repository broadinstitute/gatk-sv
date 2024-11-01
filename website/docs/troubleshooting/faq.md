---
title: FAQ
slug: faq
sidebar_position: 0
---

Please consult the following resources for additional troubleshooting guides:

- [Troubleshooting GATK-SV (article)](https://gatk.broadinstitute.org/hc/en-us/articles/5334566940699-Troubleshooting-GATK-SV)
- [Troubleshooting GATK-SV Error Messages on Terra (video)](https://www.youtube.com/watch?v=3UVV03H9p1w)

### VM runs out of memory or disk

- Default pipeline settings are tuned for batches of 150 samples. 
  Larger batches or cohorts may require additional VM resources. 
  Most runtime attributes can be modified through 
  the `RuntimeAttr` inputs. These can be formatted like this in an input json:

  ```json
  "MyWorkflow.runtime_attr_override": {
    "disk_gb": 100,
    "mem_gb": 16
  },
  ```
  
  Note that a subset of the struct attributes can be specified. 
  See [wdl/Structs.wdl](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/Structs.wdl) for available attributes.

### Calculated read length causes error in MELT workflow

Example error message from `GatherSampleEvidence.MELT.GetWgsMetrics`:

> Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 
> The requested index 701766 is out of counter bounds. 
> Possible cause of exception can be wrong READ_LENGTH 
> parameter (much smaller than actual read length)
 

This error message was observed for a sample with an average 
read length of 117, but for which half the reads were of length 
90 and half were of length 151. As a workaround, override the 
calculated read length by providing a `read_length` input of 151 
(or the expected read length for the sample in question) to `GatherSampleEvidence`.