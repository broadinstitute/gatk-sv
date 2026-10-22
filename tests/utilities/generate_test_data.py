import argparse
import transformers
import json
import logging
import os
import shutil

from pathlib import Path
from dataclasses import dataclass
from google.cloud import storage
from typing import Callable, List, Type, Union


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


@dataclass
class Region:
    chr: str
    start: int
    end: int


@dataclass
class Handler:
    transformer: Union[transformers.BaseTransformer, Type[transformers.BaseTransformer]]
    callback: Callable[[str, ...], dict]


SUBJECT_WORKFLOW_INPUTS = {
    "GatherSampleEvidence": {
        "bam_or_cram_file": Handler(
            transformers.CramDownsampler,
            lambda cram, index: {"bam_or_cram_file": cram, "bam_or_cram_index": index}
        ),
        "preprocessed_intervals": Handler(
            transformers.BedToIntervalListConverter,
            lambda x: {"preprocessed_intervals": x, "melt_metrics_intervals": x}
        ),
        "sd_locs_vcf": Handler(
            transformers.VcfDownsampler,
            lambda x: {"sd_locs_vcf": x}
        ),
        "primary_contigs_list": Handler(
            transformers.PrimaryContigsDownsampler,
            lambda x, _: {"primary_contigs_list": x}
        ),
        "primary_contigs_fai": Handler(
            transformers.PrimaryContigsDownsampler,
            lambda x, _: {"primary_contigs_fai": x}
        ),
        "wham_include_list_bed_file": Handler(
            transformers.BedDownsampler,
            lambda x: {"wham_include_list_bed_file": x}
        ),
        "manta_region_bed": Handler(
            transformers.PrimaryContigsDownsampler,
            lambda bed_gz, bed_gz_tbi: {"manta_region_bed": bed_gz}
        )
    },
    # Inputs of sv-shell's `src/sv_shell/gather_sample_evidence.sh`.
    # sv-shell input JSON keys are not prefixed with a workflow name,
    # hence, the workflow name needs to be provided via `--workflow-name`.
    # Unless `--skip-reference-downsampling` is provided, the reference is downsampled
    # (see `transformers.MiniReference`) and the coordinates in all the downsampled files
    # are lifted over to the downsampled reference.
    "gather_sample_evidence": {
        "reference_fasta": Handler(
            transformers.ReferenceDownsampler,
            lambda fasta, fai, dict, alt, amb, ann, bwt, pac, sa: {
                "reference_fasta": fasta, "reference_index": fai, "reference_dict": dict,
                "reference_bwa_alt": alt, "reference_bwa_amb": amb, "reference_bwa_ann": ann,
                "reference_bwa_bwt": bwt, "reference_bwa_pac": pac, "reference_bwa_sa": sa}
        ),
        "bam_or_cram_file": Handler(
            transformers.CramDownsampler,
            lambda cram, index: {"bam_or_cram_file": cram, "bam_or_cram_index": index}
        ),
        "preprocessed_intervals": Handler(
            transformers.IntervalListDownsampler,
            lambda x: {"preprocessed_intervals": x}
        ),
        "sd_locs_vcf": Handler(
            transformers.VcfDownsampler,
            lambda x: {"sd_locs_vcf": x}
        ),
        "primary_contigs_list": Handler(
            transformers.PrimaryContigsDownsampler,
            lambda x, _: {"primary_contigs_list": x}
        ),
        "primary_contigs_fai": Handler(
            transformers.FaiDownsampler,
            lambda x: {"primary_contigs_fai": x}
        ),
        "wham_include_list_bed_file": Handler(
            transformers.BedDownsampler,
            lambda x: {"wham_include_list_bed_file": x}
        ),
        "manta_region_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"manta_region_bed": bed_gz, "manta_region_bed_index": bed_gz_tbi}
        ),
        "mei_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"mei_bed": bed_gz}
        )
    },
    # Inputs of sv-shell's `src/sv_shell/evidence_qc.sh`; see `gather_sample_evidence` for details.
    # Note that `count_files` are outputs of `gather_sample_evidence`, hence, you may want to write
    # their downsampled files to a separate directory using `--output-subdir`.
    "evidence_qc": {
        "count_files": Handler(
            transformers.IntervalListDownsampler,
            lambda tsv_gz, tsv_gz_tbi: {"count_files": tsv_gz}
        ),
        "bincov_matrix": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"bincov_matrix": bed_gz}
        ),
        "reference_dict": Handler(
            transformers.SequenceDictDownsampler,
            lambda x: {"reference_dict": x}
        ),
        "wgd_scoring_mask": Handler(
            transformers.BedDownsampler,
            lambda x: {"wgd_scoring_mask": x}
        ),
        "genome_file": Handler(
            transformers.FaiDownsampler,
            lambda x: {"genome_file": x}
        )
    }
}


WORKFLOW_INPUTS_TO_DROP = {
    "GatherSampleEvidence": [
        "melt_docker",
        "melt_metrics_intervals",
        "melt_standard_vcf_header"
    ]
}


# The inputs that are the reference or its index files. When the reference is not downsampled
# (see `--skip-reference-downsampling`), these inputs are left unchanged (i.e., neither downsampled nor copied).
REFERENCE_INPUTS = {
    "gather_sample_evidence": [
        "reference_fasta",
        "reference_index",
        "reference_dict",
        "reference_bwa_alt",
        "reference_bwa_amb",
        "reference_bwa_ann",
        "reference_bwa_bwt",
        "reference_bwa_pac",
        "reference_bwa_sa"
    ],
    "evidence_qc": [
        "reference_dict"
    ]
}


def parse_target_regions(input_filename):
    regions = []
    with open(input_filename, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) != 3:
                raise ValueError(
                    f"Invalid line in {input_filename}. Expected a line with three columns, "
                    f"chr, start, and stop positions, found: {repr(line.strip())}")
            regions.append(Region(str(cols[0]), int(cols[1]), int(cols[2])))
    return regions


def localize_file(input_filename, output_filename):
    """
    Returns the filename of a local copy of the given input file;
    files that are already on the local filesystem are used in place.
    """
    if os.path.isfile(input_filename):
        logging.info(f"File {input_filename} is a local file, skipping localization.")
        return str(input_filename)
    if os.path.isfile(output_filename):
        logging.info(f"File {input_filename} exists locally, skipping localization.")
        return str(output_filename)
    if input_filename.startswith("gs://"):
        logging.info(f"Localizing from GCP; blob: {input_filename} ...")
        download_blob_from_gs(input_filename, output_filename)
        logging.info(f"Finished localizing blob {input_filename}")
        return str(output_filename)
    else:
        raise NotImplementedError()


def initialize_transformers(
        working_dir: str,
        reference_fasta: str,
        reference_index: str,
        sequence_dict_filename: str,
        picard_path: str,
        input_root: str = None,
        mini_reference: transformers.MiniReference = None,
        output_subdirs: dict = None):
    output_subdirs = output_subdirs or {}
    for _, inputs in SUBJECT_WORKFLOW_INPUTS.items():
        for input_var, handler in inputs.items():
            handler.transformer = handler.transformer(
                working_dir=working_dir,
                callback=handler.callback,
                reference_fasta=reference_fasta,
                reference_index=reference_index,
                sequence_dict_filename=sequence_dict_filename,
                picard_path=picard_path,
                input_root=input_root,
                mini_reference=mini_reference,
                output_subdir=output_subdirs.get(input_var)
            )


def split_input_key(key: str, workflow_name: str = None):
    """
    Returns the workflow name, the input variable name, and the prefix of the given input key.
    """
    if "." in key:
        # Example of the following split:
        # key="a.b.c" --> workflow_name="a.b", input_var="c"
        key_workflow_name, input_var = key.rsplit(".", maxsplit=1)
        return key_workflow_name, input_var, f"{key_workflow_name}."
    # Keys of sv-shell input JSON files are not prefixed with a workflow name.
    return workflow_name, key, ""


def get_local_files(workflow_inputs: dict) -> List[str]:
    values = []
    for v in workflow_inputs.values():
        values.extend(v if isinstance(v, list) else [v])
    return [v for v in values if isinstance(v, str) and os.path.isfile(v)]


def uses_mini_reference(workflow_inputs: dict, workflow_name: str = None) -> bool:
    """
    Returns True if any of the inputs is a reference input (see `REFERENCE_INPUTS`) with a handler,
    hence, the reference should be downsampled to a mini reference.
    """
    for k in workflow_inputs.keys():
        key_workflow_name, input_var, _ = split_input_key(k, workflow_name)
        if input_var in SUBJECT_WORKFLOW_INPUTS.get(key_workflow_name, {}) \
                and is_reference_input(key_workflow_name, input_var):
            return True
    return False


def parse_output_subdir(value: str):
    input_var, subdir = value.split("=", maxsplit=1)
    return input_var, subdir


def copy_local_file(working_dir: str, input_root: str, filename: str) -> str:
    """
    Copies the given file to the working directory following the same directory hierarchy
    as the file has relative to the `input_root` directory.
    """
    relative_path = Path(filename).resolve().relative_to(Path(input_root).resolve())
    output_filename = Path(working_dir).resolve().joinpath(relative_path)
    output_filename.parent.mkdir(parents=True, exist_ok=True)
    logging.info(f"Copying {filename} to {output_filename}.")
    shutil.copyfile(filename, output_filename)
    return str(output_filename)


def is_reference_input(workflow_name: str, input_var: str) -> bool:
    return input_var in REFERENCE_INPUTS.get(workflow_name, [])


def update_workflow_json(
        working_dir: str, input_filename: str, output_filename: str, output_filename_prefix: str,
        regions: List[Region], bucket_name: str = None, blob_name: str = None,
        workflow_name: str = None, input_root: str = None, skip_reference_downsampling: bool = False):
    with open(input_filename, "r") as f:
        workflow_inputs = json.load(f)

    updated_workflow_inputs = {}

    # Inputs are processed in two passes: first the inputs with a handler, then the rest.
    # Hence, the inputs set by the handler of another input (e.g., the index of a downsampled CRAM file)
    # are not processed as unchanged inputs.
    for k, v in workflow_inputs.items():
        k_workflow_name, input_var, key_prefix = split_input_key(k, workflow_name)
        handler = SUBJECT_WORKFLOW_INPUTS.get(k_workflow_name, {}).get(input_var)
        if handler is None or (skip_reference_downsampling and is_reference_input(k_workflow_name, input_var)):
            continue

        logging.info(f"Processing input {k}.")
        # An input can be a list of files (e.g., `count_files`), where each file is downsampled separately.
        updated_files = []
        for x in (v if isinstance(v, list) else [v]):
            workflow_input_local_filename = localize_file(x, Path(working_dir).joinpath(Path(x).name))
            updated_files.append(handler.transformer.transform(
                input_filename=workflow_input_local_filename,
                output_prefix=output_filename_prefix,
                regions=regions
            ))

        for varname in updated_files[0].keys():
            filenames = [x[varname] for x in updated_files]
            if bucket_name is not None and blob_name is not None:
                for i, filename in enumerate(filenames):
                    logging.info(f"Uploading downsampled file {filename} to bucket {bucket_name}.")
                    filenames[i] = upload_to_gs_blob(filename, bucket_name, blob_name)
                    logging.info(f"Finished uploading {filename}.")
            updated_workflow_inputs[f"{key_prefix}{varname}"] = filenames if isinstance(v, list) else filenames[0]

    for k, v in workflow_inputs.items():
        k_workflow_name, input_var, _ = split_input_key(k, workflow_name)
        if k in updated_workflow_inputs:
            continue
        if k_workflow_name in WORKFLOW_INPUTS_TO_DROP and input_var in WORKFLOW_INPUTS_TO_DROP[k_workflow_name]:
            logging.info(f"Dropping {k}.")
        elif skip_reference_downsampling and is_reference_input(k_workflow_name, input_var):
            updated_workflow_inputs[k] = v
            logging.info(f"Leaving reference input {k} unchanged.")
        elif input_root is not None and isinstance(v, str) and os.path.isfile(v):
            # Local files without a handler are copied, so the output directory contains all the input files.
            updated_workflow_inputs[k] = copy_local_file(working_dir, input_root, v)
        else:
            updated_workflow_inputs[k] = v
            logging.info(f"Leaving {k} unchanged.")

    # Keep the order of the inputs as in the input JSON.
    ordered_keys = [k for k in workflow_inputs.keys() if k in updated_workflow_inputs]
    ordered_keys += [k for k in updated_workflow_inputs.keys() if k not in workflow_inputs]
    updated_workflow_inputs = {k: updated_workflow_inputs[k] for k in ordered_keys}

    logging.info(f"Creating output JSON {output_filename}.")
    with open(output_filename, "w") as f:
        json.dump(updated_workflow_inputs, f, indent=4)
    logging.info(f"Finished creating output JSON {output_filename}.")


def upload_to_gs_blob(filename, bucket_name, blob_name):
    path = Path(filename)
    blob_name = blob_name + "/" + path.stem + "".join(path.suffix)
    blob = storage.Client().bucket(bucket_name).blob(blob_name)
    blob.upload_from_filename(filename)
    return "gs://" + bucket_name + "/" + blob_name


def download_blob_from_gs(gs_link, local_filename):
    bucket_name, blob_name = gs_link.split("/", maxsplit=3)[2:]
    blob = storage.Client().bucket(bucket_name).blob(blob_name)
    blob.download_to_filename(local_filename)


def main():
    parser = argparse.ArgumentParser(
        description="This is a utility script to downsample the inputs of a workflow "
                    "in order to run the workflow faster or prepare data for unit testing."
                    "In addition to other inputs, the script takes a JSON file containing the inputs to a workflow, "
                    "downsamples the inputs according to the defined rules (see `SUBJECT_WORKFLOW_INPUTS`), "
                    "pushes the downsampled files to a given cloud storage, and creates a new JSON "
                    "file with the updated downsampled inputs."
                    "This script needs samtools version 1.19.2 or newer installed and added to PATH; "
                    "downsampling a reference (e.g., for sv-shell's gather_sample_evidence) "
                    "also needs bwa installed and added to PATH.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-w", "--working-dir",
        default=os.getcwd(),
        help="Sets the working directory where the downsampled files will be stored before pushed to a cloud storage."
    )

    parser.add_argument(
        "input_workflow_json",
        help="Sets a JSON filename containing the inputs to a workflow."
    )

    parser.add_argument(
        "--output-workflow-json",
        help="Sets a JSON filename containing the updated input arguments from the input JSON."
             "The default value is a JSON file created in the working directory with the same "
             "name as the input JSON with an added prefix."
    )

    parser.add_argument(
        "picard_path",
        help="Sets the absolute path to `picard.jar`."
             "You may download picard.jar from `https://github.com/broadinstitute/picard/releases`."
    )

    parser.add_argument(
        "--output-filename-prefix",
        default="downsampled_",
        help="Sets a prefix to be added to all the output files generated."
    )

    parser.add_argument(
        "--bucket-name",
        help="Sets the cloud bucket name where the downsampled files will be pushed. "
             "The script skips uploading to cloud if a value for this argument is not provided."
    )

    parser.add_argument(
        "--blob-name",
        help="Sets the cloud blob name where the downsampled files will be pushed."
             "The script skips uploading to cloud if a value for this argument is not provided."
    )

    parser.add_argument(
        "--workflow-name",
        choices=[k for k in SUBJECT_WORKFLOW_INPUTS.keys() if "." not in k],
        help="Sets the workflow name for an input JSON whose keys are not prefixed with "
             "a workflow name, e.g., the input JSON of an sv-shell script "
             "(`gather_sample_evidence` for `src/sv_shell/gather_sample_evidence.sh`)."
    )

    parser.add_argument(
        "--skip-reference-downsampling",
        action="store_true",
        help="If provided, the reference is not downsampled; the reference inputs (e.g., the reference "
             "FASTA, its index and dictionary, and its BWA index files) are left unchanged, and the "
             "coordinates of the downsampled files are the same as those of the input files. "
             "Otherwise, for the workflows that take a reference as input (e.g., sv-shell's "
             "gather_sample_evidence), the reference is downsampled to the target regions "
             "and the coordinates in all the downsampled files are lifted over to the downsampled reference."
    )

    parser.add_argument(
        "--input-root",
        help="Sets a directory containing the local input files referenced in the input JSON. "
             "The output files are written to the working directory following the same directory "
             "hierarchy as the input files have relative to this directory; the local input files "
             "that are not downsampled are copied. The default value is the deepest common "
             "directory of the local input files' directories."
    )

    parser.add_argument(
        "--output-subdir",
        action="append",
        default=[],
        type=parse_output_subdir,
        help="Sets the subdirectory of the working directory where the downsampled files of an input "
             "are written, formatted as `input_name=subdirectory`; e.g., "
             "`count_files=intermediates_HG00512/gather_sample_evidence` for an input that is "
             "the output of an upstream task. Can be provided multiple times. The inputs without "
             "a subdirectory set follow the directory hierarchy of the input files (see `--input-root`)."
    )

    this_script_folder = os.path.dirname(os.path.abspath(__file__))
    gatk_sv_path = os.path.dirname(os.path.dirname(this_script_folder))

    parser.add_argument(
        "--target-regions",
        default=os.path.join(gatk_sv_path, "tests", "utilities", "default_downsampling_regions.bed"),
        help="Sets a BED filename containing target regions for downsampling, "
             "such that the downsampled files contains data from the input files overlapping these regions."
    )

    parser.add_argument(
        "--reference-fasta",
        default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
        help="Set the path to reference fasta file. "
    )

    parser.add_argument(
        "--reference-fasta-index",
        default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        help="Set the path to index file of the reference fasta file. "
    )

    parser.add_argument(
        "--reference-dict",
        default="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
        help="Set the path to the reference dictionary file."
    )

    args = parser.parse_args()

    regions = parse_target_regions(args.target_regions)
    logging.info(f"Found {len(regions)} target regions for downsampling.")

    output_workflow_json = args.output_workflow_json
    if not output_workflow_json:
        output_workflow_json = os.path.join(
            args.working_dir,
            f"{args.output_filename_prefix}{Path(args.input_workflow_json).name}"
        )

    Path(args.working_dir).mkdir(parents=True, exist_ok=True)

    with open(args.input_workflow_json, "r") as f:
        workflow_inputs = json.load(f)

    input_root = args.input_root
    local_files = get_local_files(workflow_inputs)
    if input_root is None and len(local_files) > 0:
        input_root = os.path.commonpath([os.path.dirname(os.path.abspath(x)) for x in local_files])

    mini_reference = None
    if not args.skip_reference_downsampling and uses_mini_reference(workflow_inputs, args.workflow_name):
        logging.info(f"Creating a mini reference from the target regions of {args.reference_fasta}.")
        reference_fasta_local_filename = localize_file(
            args.reference_fasta, Path(args.working_dir).joinpath(Path(args.reference_fasta).name))
        mini_reference = transformers.MiniReference(reference_fasta_local_filename, regions, args.working_dir)

    sequence_dict_local_filename = localize_file(
        args.reference_dict, Path(args.working_dir).joinpath(Path(args.reference_dict).name))
    initialize_transformers(args.working_dir, args.reference_fasta, args.reference_fasta_index,
                            sequence_dict_local_filename, args.picard_path, input_root, mini_reference,
                            dict(args.output_subdir))

    update_workflow_json(
        working_dir=args.working_dir,
        input_filename=args.input_workflow_json,
        output_filename=output_workflow_json,
        output_filename_prefix=args.output_filename_prefix,
        regions=regions,
        bucket_name=args.bucket_name,
        blob_name=args.blob_name,
        workflow_name=args.workflow_name,
        input_root=input_root,
        skip_reference_downsampling=args.skip_reference_downsampling)

    if mini_reference is not None:
        mini_reference.cleanup()

    logging.info("All process finished successfully.")


if __name__ == '__main__':
    main()
