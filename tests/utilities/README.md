
# Installation

1. Create a virtual environment:

   ```shell
   virtualenv .venv
   ```
   
2. Activate the virtual environment:

   ```shell
   source .venv/bin/activate
   ```
   
3. Install requirements:

   ```shell
   pip install -r requirements.txt
   ```
# Downsample single-sample pipeline inputs

Run these from the repository root.

Install `bwa` and `bgzip`, which are used to index the downsampled reference.

```shell
brew install bwa htslib
```

Set the directory to download the raw inputs to.

```shell
INPUTS=svshell_downsampling_inputs
```

Create the directory.

```shell
mkdir -p $INPUTS
```

Set the directory to write the downsampled data to.

```shell
OUTPUT=$PWD/downsampled_data
```

Download the raw inputs (~49 GB).

```shell
gcloud storage rsync -r gs://broad-dsde-methods-vj/sv-shell-downsampling-inputs $INPUTS
```

Generate the downsampled inputs in `$OUTPUT` and their input JSON (~10 min).

```shell
python tests/utilities/generate_test_data.py \
  $INPUTS/single_sample_pipeline.json \
  $INPUTS/picard.jar \
  --workflow-name single_sample_pipeline \
  --working-dir $OUTPUT \
  --target-regions tests/utilities/sv_shell_downsampling_regions.bed \
  --reference-fasta $INPUTS/reference/Homo_sapiens_assembly38.fasta \
  --reference-fasta-index $INPUTS/reference/Homo_sapiens_assembly38.fasta.fai \
  --reference-dict $INPUTS/reference/Homo_sapiens_assembly38.dict \
  --cnmops-contigs chr16 chr20 \
  --gcnv-contigs chr16 chr20 chrX chrY \
  --ref-panel-samples tests/utilities/sv_shell_ref_panel_samples.list \
  --output-root-in-json /test_inputs \
  --output-workflow-json src/sv_shell/sample_inputs/single_sample_pipeline_downsampled.json
```

Run the pipeline in the sv-shell container with `$OUTPUT` mounted at `/test_inputs`; see [VM setup](https://broadinstitute.github.io/gatk-sv/docs/svshell/vm-setup) for more details.

```shell
docker run --rm -t -d --entrypoint /bin/bash --platform linux/amd64 \
  -v "$OUTPUT:/test_inputs/" \
  -v "$PWD/src/sv_shell:/opt/sv_shell/" \
  us.gcr.io/broad-dsde-methods/gatk-sv/sv-shell:2026-10-08-v1.1.1-a55c498f
```

```shell
docker exec -it {CONTAINER_ID} /bin/bash
```

```shell
cd /opt/sv_shell
export SV_SHELL_BASE_DIR="/wd"
export TMPDIR="/wd/tmp"
mkdir -p /wd/tmp
```


```shell
{ time bash single_sample_pipeline.sh sample_inputs/single_sample_pipeline_downsampled.json; } \
  > >(tee "${SV_SHELL_BASE_DIR}/stdout.log") \
  2> >(tee "${SV_SHELL_BASE_DIR}/stderr.log" >&2)
```
