## Overview

This directory is a pure bash reimplementation of the GATK-SV WDL workflows, 
running the **single-sample ("case") mode** of GATK-SV outside of Cromwell. 

Please refer to [this page](https://broadinstitute.github.io/gatk-sv/docs/category/execution/svshell)
for a detailed documentation.

Each `*.sh` script mirrors one WDL workflow or task one-to-one. 
`single_sample_pipeline.sh` is the top-level driver: 
it chains the module scripts together in the same order as GATK-SV's `GATKSVPipelineSingleSample.wdl`, 
calling a case sample against a pre-built reference panel 
(bincov matrix, PE/SR/SD files, gCNV models, 
standardized Manta/Wham/Scramble/Dragen VCF tars, depth BED files, etc.) 
to produce a genotyped, filtered, and annotated per-sample SV VCF.

## Running a module

Every script, including the top-level pipeline, follows the same calling convention:

```
bash <script>.sh <inputs.json> [outputs.json] [output_dir]
```

* `inputs.json` supplies the task's input parameters; 
* `outputs.json` is written with the paths to the task's output files (mirroring WDL task outputs), 
defaulting to `output_dir/output.json`; 
* `output_dir` holds outputs and logs, defaulting to a fresh `mktemp -d` under `SV_SHELL_BASE_DIR`. 
 
Two environment variables must be set before running: 
* `SV_SHELL_BASE_DIR` (scratch/output root with sufficient disk space)
* `TMPDIR`.


Because every module shares this contract, 
each can be run and tested independently of the full pipeline. 
See `sample_inputs/` for a working `inputs.json` per script.


## Directory layout

- `*.sh`: one script per WDL workflow/task, each following the 
  `inputs.json` / `outputs.json` / `output_dir` contract above. 
  `single_sample_pipeline.sh` is the entrypoint; the rest are modules it invokes, 
  except `sv_cluster.sh`, which is a shared clustering helper called from several modules.

- `sample_inputs/`: a working example `inputs.json` for every script.

- `sample_outputs/`: example `outputs.json` for a subset of the larger modules.
