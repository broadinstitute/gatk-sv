---
title: Overview
description: Overview
sidebar_position: 1
slug: overview
---

# Overview

SV-Shell is a platform-agnostic implementation of the single-sample pipeline: 
it runs end-to-end in a single Docker image, 
so it works on any platform that can run Docker images, 
such as Azure, AWS, GCP, and on-premises infrastructure.

It genotypes structural variants in one case sample against a pre-built reference panel: 
collecting per-sample evidence (read depth, PE/SR, and raw SV calls from multiple callers), 
merging that evidence with the batch/panel, 
clustering and genotyping the calls, 
resolving complex multi-breakpoint SVs, 
scoring and filtering genotypes with an ML-based recalibrator, 
and annotating the result. 
It supports two evidence configurations, 
using either a combination of callers or DRAGEN's built-in calls, 
and produces a single genotyped, filtered, 
and annotated per-sample VCF along with QC metrics.

SV-Shell is implemented in bash and runs without a WDL execution engine, 
without any Terra-platform dependency, 
or any GCP dependency. 


## Evidence Source Configurations

| | Single-Sample Pipeline without DRAGEN | Single-Sample Pipeline with DRAGEN |
|---|---------------------------------------|------------------------------------|
| SV Callers | Manta, Scramble, Wham                 | DRAGEN SV, Scramble, Wham          |
| CNV Callers | cn.MOPS, gCNV                         | DRAGEN CNV                         |


## Invocation

`single_sample_pipeline.sh` is the entrypoint and drives all other scripts in order. Every script shares one calling convention:

```
bash <script>.sh <inputs.json> [outputs.json] [output_dir]
```

- **inputs.json**: the task's parameters (file paths, flags, tuning values).
- **outputs.json**: written by the script with paths to its key output files (mirrors a WDL task's outputs); defaults to `output_dir/output.json`.
- **output_dir**: holds outputs and logs; defaults to a fresh temp directory under `SV_SHELL_BASE_DIR`.

`SV_SHELL_BASE_DIR` and `TMPDIR` must be set in the environment before running anything. Because every module shares this contract, any script can be run and tested on its own using the matching file in `sample_inputs/`.