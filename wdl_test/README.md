This directory contains [carrot](https://github.com/broadinstitute/carrot)
tests for the GAKT-SV pipeline's WDLs; the tests are organized in folders
containing `carrot` resources (e.g., evaluation WDL, default/test inputs).
Additionally, a utility script, `carrot_help.py`, is provided that 
automates defining tests to Carrot, and running and checking their execution 
status. Generally, with tests organized in a particular folder hierarchy,
with a single call to the utility script 
(`python carrot_helper.py test run ./*`) every step from defining and 
running tests are automatically executed, ideally simplifying 
defining/running Carrot tests without requiring domain-specific expertise.

## Organize Tests in Directories

Test cases for a WDL need to be organized in directories where each 
contain separate subdirectories for test cohorts each containing 
`carrot` resources. For instance, the `ExpansionHunterDenovo.wdl`
performs case-control analysis and outlier detection on a set 
of BAM files based on their short-tandem repeat (STR) profiles.
In order to test this WDL for case-control analysis using a cohort
of simulated data, `carrot` resources need to be organized as
the following for the `carrot_helper` utility to automatically 
setup, run, and track `carrot` tests. 

```shell
├── wdl
│   └── ExpansionHunterDenovo.wdl
└── wdl_test
   └── ExpansionHunterDenovo
       └── casecontrol
           ├── simulated_data
           │   ├── eval_input.json
           │   └── test_input.json
           ├── eval.wdl
           ├── eval_input_defaults.json
           └── test_input_defaults.json
```

Accordingly: 

- Create a folder with the _same name_ as the WDL file containing the
workflow you want to test (`ExpansionHunterDenovo` for 
`ExpansionHunterDenovo.wdl` in this example).


- Create a separate folder for every evaluation you want to perform (e.g., 
the `casecontrol` folder to evaluate `ExpansionHunterDenovo.wdl`'s
case-control analysis on the STR profiles of input BAM files). 
While all assertions can be part of a single evaluation, it is generally 
a good practice to break assertions into smaller atomic evaluations. 


- Inside every evaluation directory, create three files: `eval.wdl`, 
`eval_input_defaults.json`, and `test_input_defaults.json`. 
The `eval.wdl` WDL receives outputs of the workflow you're testing and 
asserts their values. The JSON files provide default inputs to the test
(`ExpansionHunterDenovo.wdl`) and `eval.wdl` WDLs. For instance, if the 
majority of the tests are running `eval.wdl` on a common docker image, 
the image name can be set in the `eval_input_defaults.json`, which can be
overridden in the tests that execute `eval.wdl` on a different docker image.


- An evaluation can be performed using different set of inputs for the 
test and evaluation workflows. For instance, in the STR analysis scenario, 
we pass 
[seven BAM files](https://github.com/VJalili/gatk-sv/blob/89e67350ea7fec8edc687011ac7308e3e1db17ff/wdl_test/ExpansionHunterDenovo/casecontrol/simulated_data/test_input.json#L4-L12)
to the `ExpansionHunterDenovo.wdl`, run the WDL, and pass 
[its output](https://github.com/VJalili/gatk-sv/blob/89e67350ea7fec8edc687011ac7308e3e1db17ff/wdl_test/ExpansionHunterDenovo/casecontrol/simulated_data/eval_input.json#L2)
along with the 
[expected output](https://github.com/VJalili/gatk-sv/blob/89e67350ea7fec8edc687011ac7308e3e1db17ff/wdl_test/ExpansionHunterDenovo/casecontrol/simulated_data/eval_input.json#L3)
to the evaluation WDL. Different combinations of inputs the test and 
evaluation workflows are grouped under separate subdirectories (e.g., 
the `simulated_data` subdirectory for `casecontrol` assertion of 
`ExpansionHunterDenovo.wdl`). The inputs for test and evaluation
WDLs are specified using two JSON files, `test_input.json` and 
`eval_input.json`, containing inputs for the test and evaluation 
WDLs respectively. The files should be located in the subdirectory of the 
test cohort (e.g, `casecontrol/simulated_data/test_input.json`).


- In order to pass any file to the WDLs via the JSON files, the files 
need to be stored on a publicly accessible Google storage bucket.


- In order to pass the output of test WDL as input to the evaluation WDL, 
the value of the key should be prefixed with `test_output:` (see `carrot`'s 
[documentation](https://github.com/broadinstitute/carrot/blob/0f616c0a9933a44bb92bc9ddbc90b81b0b532de6/UserGuide.md#-mapping-test-outputs-to-eval-inputs)).
For instance:

  ```json
  "EvalCaseControlLocus.multisample_profile": "test_output:EHdnSTRAnalysis.multisample_profile",
  ```


## Carrot Helper

The `carrot_helper` utility script automates few routine task for 
running and updating Carrot tests. This script is not a replacement
for `carrot_cli` or Carrot's API that have more expressive power,
wider functionality, and generalization than `carrot_helper`. 

### Setup

1. Install `carrot_cli`:
    - Install the `dev` version of [`carrot_cli`](https://github.com/broadinstitute/carrot_cli) 
   as the following. We install the `dev` since `carrot_helper` leverages 
   unreleased feature of `carrot_cli`. 
    
        ```shell
        git clone https://github.com/broadinstitute/carrot_cli/
        pip install -r dev-requirements.txt
        pip install -e .
        ```
      
    - [Configure `carrot_cli`]((https://github.com/broadinstitute/carrot/blob/master/UserGuide.md#-carrot-cli)):
   configure it to access a [Carrot server](https://github.com/broadinstitute/carrot) 
   and set your email address.


2. Install latest version of 
[`womtool`](https://github.com/broadinstitute/cromwell/releases).


3. Setup `carrot_helper.py` by executing the following command providing
values for its prompts:

    ```shell
    $ cd gatk-sv/wdl_test
    $ python carrot_helper.py config
    ```
   
    Carrot fetches the test and evaluation WDLs for every test from 
a publicly accessible GitHub repository. Therefore, in order to define/update
tests, `carrot_helper` requires to know the GitHub repository and the git 
branch where the test and evaluation WDLs are available. If you want to run 
existing tests, you may use `https://github.com/broadinstitute/gatk-sv` and
`master` for repository and branch respectively. If you are developing 
a carrot test for a WDL, then you may set the repository to your fork
of `github.com/broadinstitute/gatk-sv` and set the branch to your feature
branch.


### Run Carrot Helper

```shell
cd wdl_test
python carrot_helper.py test run ./*
```
_Note that the script should be invoked from the `wdl_test` directory._

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

### Reusable Resources
The `carrot_helper.py` persists any metadata about the carrot resources it
creates (e.g., 
[pipeline](https://github.com/broadinstitute/carrot/blob/master/UserGuide.md#-pipeline),
[template](https://github.com/broadinstitute/carrot/blob/master/UserGuide.md#-template), 
[test](https://github.com/broadinstitute/carrot/blob/master/UserGuide.md#-test), 
[result](https://github.com/broadinstitute/carrot/blob/master/UserGuide.md#-result)
and any necessary mapping between them) in the `.carrot_pipelines.json`. 

The `.carrot_pipelines.json` file tracked on git contains metadata belonging
to the `carrot` resources defined for tests and WDLs available from the 
`master` branch of the 
[`github.com/broadinstitute/gatk-sv`](https://github.com/broadinstitute/gatk-sv) 
repository on a `carrot` server maintained for internal use at the Broad 
institute. You may use this file to run and updated (read the following)
tests if you have access to Broad's VPN. Otherwise, you may remove or rename
the `.carrot_pipelines.json` file, **without tracking the changes on git**,
and let the `carrot_helper.py` create resources on the `carrot` server for
the repository and branch [you have configured](#setup-carrot-helper). 

`carrot_helper.py` automatically initializes and updates the 
`.carrot_pipelines.json`. When `carrot_helper.py test run` is invoked, 
the script traverses the `wdl_test` and initializes/updates `carrot` 
resources if any of the test or evaluations WDLs or their inputs are 
changed. Carrot reads test and evaluation WDLs from github; therefore, 
make sure you commit and push changes to your branch when updating 
test and evaluation WDLs. 


### Carrot Report
Carrot can pass the output of an evaluation workflow to a Jupyter notebook,
which enables more in-depth evaluations/assertions and visualizations. 
visualization. In general, this requires defining a template notebook 
(ideally separate notebooks for each test to have test-specific visualization), 
defining a `report` in carrot and mapping a template to the report.
Please refer to [Carrot documentation for details.](https://github.com/broadinstitute/carrot/blob/48c58446d4fb044cabbdafe8962b67ee511b483a/UserGuide.md#-2-define-a-report-in-carrot)
The `carrot_helper` does not currently support defining `report`.


### Current limitations

Carrot is under active development and new functionalities emerge
as new versions are released. There are a few functionalities that are
under development and not yet released that impact the workflows that
can be tested using `carrot`. Specifically, Carrot does not currently 
support relative imports in WDL files (i.e., importing workflow via 
a WDL file is provided via the `--imports` argument of `cromwell`). 
A workaround to is to host required imports on a Google cloud storage 
bucket and import using the object's URL. However, this would require 
modifying all the WDLs of the GATK-SV pipeline. The carrot team is working 
on supporting an `--imports`-like functionality in carrot. 

Additionally, `carrot` do not currently support `Array` type outputs 
(e.g., `Array[File]`). In other words, the array type outputs of a 
test WDL cannot be passed to evaluation WDLs for assertions. A workaround 
is to encapsulate array output in a zip archive, hence the test WDL outputs
a single file, and extract the content of zip in the eval WDL. This workaround
would require a significant modification to GATK-SV pipeline workflows, hence
we currently do not assert array type outputs. 
