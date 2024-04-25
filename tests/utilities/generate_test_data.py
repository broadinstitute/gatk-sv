import argparse
import transformers
import json
import logging
import os

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
            lambda x: {"primary_contigs_list": x}
        ),
        "primary_contigs_fai": Handler(
            transformers.PrimaryContigsDownsampler,
            lambda x: {"primary_contigs_fai": x}
        ),
        "wham_include_list_bed_file": Handler(
            transformers.BedDownsampler,
            lambda x: {"wham_include_list_bed_file": x}
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
    if os.path.isfile(output_filename):
        logging.info(f"File {input_filename} exists locally, skipping localization.")
        return
    if input_filename.startswith("gs://"):
        logging.info(f"Localizing from GCP; blob: {input_filename} ...")
        download_blob_from_gs(input_filename, output_filename)
        logging.info(f"Finished localizing blob {input_filename}")
    else:
        raise NotImplementedError()


def initialize_transformers(
        working_dir: str,
        reference_fasta: str,
        reference_index: str,
        sequence_dict_filename: str,
        picard_path: str):
    for _, inputs in SUBJECT_WORKFLOW_INPUTS.items():
        for _, handler in inputs.items():
            handler.transformer = handler.transformer(
                working_dir=working_dir,
                callback=handler.callback,
                reference_fasta=reference_fasta,
                reference_index=reference_index,
                sequence_dict_filename=sequence_dict_filename,
                picard_path=picard_path
            )


def update_workflow_json(
        working_dir: str, input_filename: str, output_filename: str, output_filename_prefix: str,
        regions: List[Region], bucket_name: str = None, blob_name: str = None):
    with open(input_filename, "r") as f:
        workflow_inputs = json.load(f)

    updated_workflow_inputs = {}

    for k, v in workflow_inputs.items():
        # Example of the following split:
        # k="a.b.c" --> workflow_name="a.b", input_var="c"
        workflow_name, input_var = k.rsplit(".", maxsplit=1)

        try:
            handler = SUBJECT_WORKFLOW_INPUTS[workflow_name][input_var]
        except KeyError:
            if workflow_name in WORKFLOW_INPUTS_TO_DROP and input_var in WORKFLOW_INPUTS_TO_DROP[workflow_name]:
                logging.info(f"Dropping {k}.")
            else:
                updated_workflow_inputs[k] = v
                logging.info(f"Leaving {k} unchanged.")
            continue

        logging.info(f"Processing input {k}.")
        workflow_input_local_filename = Path(working_dir).joinpath(Path(v).name)
        localize_file(v, workflow_input_local_filename)
        updated_files = handler.transformer.transform(
            input_filename=workflow_input_local_filename,
            output_prefix=output_filename_prefix,
            regions=regions
        )

        for varname, filename in updated_files.items():
            input_key = f"{workflow_name}.{varname}"
            updated_workflow_inputs[input_key] = filename
            if bucket_name is not None and blob_name is not None:
                logging.info(f"Uploading downsampled file {filename} to bucket {bucket_name}.")
                blob = upload_to_gs_blob(filename, bucket_name, blob_name)
                logging.info(f"Finished uploading {filename}.")
                updated_workflow_inputs[input_key] = blob
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
                    "This script needs samtools version 1.19.2 or newer installed and added to PATH.",
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

    sequence_dict_local_filename = Path(args.working_dir).joinpath(Path(args.reference_dict).name)
    localize_file(args.reference_dict, sequence_dict_local_filename)
    initialize_transformers(args.working_dir, args.reference_fasta, args.reference_fasta_index,
                            sequence_dict_local_filename, args.picard_path)

    update_workflow_json(
        working_dir=args.working_dir,
        input_filename=args.input_workflow_json,
        output_filename=output_workflow_json,
        output_filename_prefix=args.output_filename_prefix,
        regions=regions,
        bucket_name=args.bucket_name,
        blob_name=args.blob_name)

    logging.info("All process finished successfully.")


if __name__ == '__main__':
    main()
