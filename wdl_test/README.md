This directory contains [Carrot](https://github.com/broadinstitute/carrot)
tests for the GAKT-SV pipeline's WDLs; the tests are organized in folders
containing Carrot resources (e.g., evaluation WDL, default/test inputs).
Additionally, a utility script, `carrot_help.py`, is provided that 
automates defining tests to Carrot, running, and checking their execution 
status. Generally, with tests organized in a particular folder hierarchy,
with a single call to the utility script 
(`python carrot_helper.py test run ./*`) every step from defining and 
running tests are automatically executed, ideally simplifying 
defining/running Carrot tests without requiring domain-specific expertise.

## Organize Tests in Directories

Consider a WDL that takes a set of BAM files as input, and it first
computes short-tandem repeat (STR) profiles for each of the BAM files, and 
then performs a comparative analysis on the determined STR profiles.
Suppose we want to test computing STR profiles and performing comparative
analysis separately using two cohorts: a small simulated cohort, and a 
large real data cohort. In order to use `carrot_helper` to setup Carrot
tests for this scenario, we need to create a directory structure as 
the following: 

```shell
├── wdl
│   └── STRAnalyzer.wdl
└── wdl_test
    └── STRAnalyzer
        ├── comparative
        │   ├── real_cohort
        │   │   ├── eval_input.json
        │   │   └── test_input.json
        │   ├── simulated_cohort
        │   │   ├── eval_input.json
        │   │   └── test_input.json
        │   ├── eval.wdl
        │   ├── eval_input_defaults.json
        │   └── test_input_defaults.json
        └── compute_profiles
            ├── real_cohort
            │   ├── eval_input.json
            │   └── test_input.json
            ├── simulated_cohort
            │   ├── eval_input.json
            │   └── test_input.json
            ├── eval.wdl
            ├── eval_input_defaults.json
            └── test_input_defaults.json
```

Accordingly: 

- Create a folder with the _same name_ as the WDL file containing the
workflow you want to test (`STRAnalyzer` for `STRAnalyzer.wdl` in 
this example).


- Create a separate folder for every evaluation you want to perform (e.g., 
in this example, the `compute_profiles` and `comparative` folders to evaluate
`STRAnalyzer`'s STR profile computation and their comparative analysis, 
respectively). While you can make all assertions as part of a single
evaluation, it is generally a good practice to break assertions into smaller
atomic evaluations. 


- Inside every evaluation directory, create three files: `eval.wdl`, 
`eval_input_defaults.json`, and `test_input_defaults.json`. 
The `eval.wdl` WDL receives outputs of the workflow you're testing and 
asserts their values. The JSON files provide default inputs to the test
(`STRAnalyzer.wdl` in the STR analysis scenario) and `eval.wdl` WDLs.


- An evaluation can be performed using different set of inputs for the 
workflow we're testing and the workflow making the assertions. For instance,
in the STR analysis scenario, we want to pass `sample1.bam` and `sample2.bam`
to the `STRAnalyzer` workflow, and pass `sample1_expected_str_profiles.json`
and `sample2_expected_str_profiles.json` along with the `STRAnalyzer`'s 
produced STR profiles to the `eval.wdl` to assert the equality of the 
produced and expected STR profiles. We group different combinations of inputs
to the test and evaluation workflows under different folders. In the 
STR analysis scenario, the `simulated_cohort` and `real_cohort`. Inside
this folder create two files: `eval_input.json` and `test_input.json`, 
containing, the inputs to the eval and test workflows respectively.


## Carrot Helper

The `carrot_helper` utility script automates few routine task for 
running and updating Carrot tests. This script is not a replacement
for `carrot_cli` or Carrot's API that have more expressive power,
wider functionality, and generalization than `carrot_helper`. 

### Setup
- Install and configure 
[`carrot_cli`](https://github.com/broadinstitute/carrot_cli). Make sure to 
install the `dev` version since `carrot_helper` leverages features not 
included in the current latest release.

    ```shell
    git clone https://github.com/broadinstitute/carrot_cli/
    pip install -r dev-requirements.txt
    pip install -e .
    ```
    Make sure to configure `carrot_cli` and set your email address.
    `carrot_cli` needs to be configured to access a
    [Carrot server](https://github.com/broadinstitute/carrot).

- Install latest version of 
[`womtool`](https://github.com/broadinstitute/cromwell/releases).


### Run Carrot Helper

```shell
cd wdl_test
python carrot_helper.py test run ./*
```
_The script should be invoked from the `wdl_test` directory._

This above command will define every test (in the above-discussed directory
structure) to Carrot, and will run them all. The information of the created 
and executed tests are persisted in `.carrot_pipelines.json` and `.runs.json`
files. 

You can specify a single test to run; for instance:

```shell
python carrot_helper.py test run STRAnalyzer/comparative/real_cohort
```

Or you may use wildcards to specify particular tests to run. For instance:

```shell
python carrot_helper.py test run STRAnalyzer/*/real_cohort
```

To check for the status of the runs, you use the following command.

```shell
python carrot_helper.py test update_status
```

### Carrot Report
Carrot can pass the output of an evaluation workflow to a Jupyter notebook,
which enables more in-depth evaluations/assertions and visualizations. 
visualization. In general, this requires defining a template notebook 
(ideally separate notebooks for each test to have test-specific visualization), 
defining a `report` in carrot and mapping a template to the report.
Please refer to [Carrot documentation for details.](https://github.com/broadinstitute/carrot/blob/48c58446d4fb044cabbdafe8962b67ee511b483a/UserGuide.md#-2-define-a-report-in-carrot)
The `carrot_helper` does not currently support defining `report`.
