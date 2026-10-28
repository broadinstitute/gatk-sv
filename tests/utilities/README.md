
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

Download the raw inputs (~49 GB).

```shell
gcloud storage rsync -r gs://broad-dsde-methods-vj/sv-shell-downsampling-inputs $INPUTS
```

Generate the downsampled inputs in `downsampled_data/` and their input JSON (~10 min).

```shell
python tests/utilities/generate_test_data.py \
  $INPUTS/single_sample_pipeline.json \
  $INPUTS/picard.jar \
  --workflow-name single_sample_pipeline \
  --working-dir downsampled_data \
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

Run the pipeline in the sv-shell container with `downsampled_data/` mounted at `/test_inputs`.

```shell
bash /opt/sv_shell/single_sample_pipeline.sh /opt/sv_shell/sample_inputs/single_sample_pipeline_downsampled.json output.json out
```
