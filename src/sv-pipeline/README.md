# SV detection pipeline
This repository contains the full script developed by Talkowski Lab to detect structural variants.

## Installation and Usage
As a snakemake workflow, the pipeline is intended to be cloned for each project, e.g.
```
git clone https://github.com/talkowski-lab/SV-Adjudicator.git MySVDiscoveryProject
```
Create and activate a new Anaconda environment with all the necessary dependencies.

```
cd MySVDiscoveryProject
conda env create -f environment.yaml
source activate sv_pipeline
```
## Module configuration and input
After cloning the pipeline, edit `config.yaml` to update the configuration as
necessary for the project, then link or copy raw data into the `data/` or
`ref/` directories. (More detail below or to come. For current testing
purposes, symlink the `data/` and `ref/` directories in
`/data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline-devel/`).

Optionally, run a dry run of `snakemake` to test the configuration, then run
the pipeline with `snakemake`.

```
$ vim config.yaml
$ ln -s ${raw_SV_calls} data/
$ cp ${reference_data} ref/
$ snakemake -np
$ snakemake
```

The pipeline can remove all files it has generated without affecting
configuration or data files.

```
$ snakemake clean
```

## Install svtk
`svtk` is required for the process of `SV-Adjudicator` and is included in this repository. To install `svtk`,
```
cd svtk
python setup.py install --user
```
For more details about `svtk`, refer to the subfolder.

## Explicit Dependencies

If you would prefer to use your own python environment, the following packages
are necessary to run the pipeline.

### Python 3

The pipeline requires Python 3. If you already have an established Python 3
environment, you can ignore this subsection. Otherwise, the recommended way to
create a Python 3 environment with the required packages is with
[Anaconda](https://www.continuum.io/downloads).

```
$ conda create -n $environment -c bioconda python=3.5 numpy scipy pysam snakemake
$ source activate $environment
```
### Snakemake
The pipeline is built with `snakemake`, which can be installed through `pip` or
`conda` with one of the two following commands.

```
$ pip install snakemake
$ conda install -c bioconda snakemake
```

`snakemake` is an excellent and recommended tool for building bioinformatics
pipelines, but a comprehensive understanding of the tool is not necessary to
run this pipeline. If interested in learning more, extended `snakemake`
documentation can be found on their [Read the Docs
page](https://snakemake.readthedocs.io/en/stable/). A
[tutorial](https://snakemake.bitbucket.io/snakemake-tutorial.html) and
[demonstrative slides](http://slides.com/johanneskoester/deck-1#/) are also
available.

### SVtools
The pipeline requires the `svtk` Python package, which is currently
available only on github.

```
$ git clone git@github.com:talkowski-lab/svtk.git
$ cd svtk
$ pip install -e .
```

### Bedtools
The pipeline requires bedtools 2.26 or later. Earlier versions may throw an
error when `bedtools merge` is passed an empty file.

### pybedtools
In order to perform per-chromosome parallelization, the master branch of
`pybedtools` is required (or at least commit `b1e0ce0`).

```
$ pip install git+git://github.com/daler/pybedtools.git@master
$ pip install git+git://github.com/daler/pybedtools.git@b1e0ce0
```

# Running the pipeline

The pipeline consists of multiple independent modules. Documentation of each
module is provided in the respective subdirectory.

0. [Preprocessing](preprocessing/README.md)
   (optional, not intended for all users)
1. [Algorithm integration](algorithm_integration/README.md)
2. [RD-test](rdtest/README.md)
3. SR-test

## Rolling your own preprocessing

The preprocessing module here is provided for reproducibility and as an
example implementation of SV algorithm standardization, but is not intended to
be generalizable to all use cases.

If you would like to implement your own preprocessing before beginning
integration and filtering, the pipeline can be bootstrapped to begin with the
integration module by providing the following input files:
