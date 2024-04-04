import argparse
import downsamplers
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
    downsampler: Union[downsamplers.BaseDownsampler, Type[downsamplers.BaseDownsampler]]
    callback: Callable[[str, ...], dict]


SUBJECT_WORKFLOW_INPUTS = {
    "GatherSampleEvidence": {
        "bam_or_cram_file": Handler(
            downsamplers.CramDownsampler,
            lambda cram, index: {"bam_or_cram_file": cram, "bam_or_cram_index": index}
        ),
        "preprocessed_intervals": Handler(
            downsamplers.IntervalListDownsampler,
            lambda x: {"preprocessed_intervals": x}
        ),
        "sd_locs_vcf": Handler(
            downsamplers.VcfDownsampler,
            lambda x: {"sd_locs_vcf": x}
        ),
        "primary_contigs_list": Handler(
            downsamplers.PrimaryContigsDownsampler,
            lambda x: {"primary_contigs_list": x}
        ),
        "primary_contigs_fai": Handler(
            downsamplers.PrimaryContigsDownsampler,
            lambda x: {"primary_contigs_fai": x}
        )
    }
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
        return
    if input_filename.startswith("gs://"):
        logging.info(f"Localizing from GCP; blob: {input_filename} ...")
        download_blob_from_gs(input_filename, output_filename)
        logging.info(f"Finished localizing blob {input_filename}")
    else:
        raise NotImplementedError()


def initialize_downsamplers(working_dir: str):
    for _, inputs in SUBJECT_WORKFLOW_INPUTS.items():
        for _, handler in inputs.items():
            handler.downsampler = handler.downsampler(working_dir, handler.callback)


def update_workflow_json(
        working_dir: str, input_filename: str, output_filename: str, output_filename_prefix: str,
        regions: List[Region], bucket_name: str = None, blob_name: str = None):
    with open(input_filename, "r") as f:
        workflow_inputs = json.load(f)

    for k, v in dict(workflow_inputs).items():
        # Example of the following split:
        # k="a.b.c" --> workflow_name="a.b", input_var="c"
        workflow_name, input_var = k.rsplit(".", maxsplit=1)

        try:
            handler = SUBJECT_WORKFLOW_INPUTS[workflow_name][input_var]
        except KeyError:
            # This workflow input is not set to be downsampled.
            continue

        logging.info(f"Processing input {k}.")
        workflow_input_local_filename = Path(working_dir).joinpath(Path(v).name)
        localize_file(v, workflow_input_local_filename)
        updated_files = handler.downsampler.downsample(workflow_input_local_filename, output_filename_prefix, regions)
        if bucket_name is not None and blob_name is not None:
            for varname, filename in updated_files.items():
                logging.info(f"Uploading downsampled file {filename} to bucket {bucket_name}.")
                blob = upload_to_gs_blob(filename, bucket_name, blob_name)
                logging.info(f"Finished uploading {filename}.")
                workflow_inputs[f"{workflow_name}.{varname}"] = blob
    logging.info(f"Creating output JSON {output_filename}.")
    with open(output_filename, "w") as f:
        json.dump(workflow_inputs, f, indent=4)
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
                    "file with the updated downsampled inputs.",
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

    args = parser.parse_args()

    regions = parse_target_regions(args.target_regions)
    logging.info(f"Found {len(regions)} target regions for downsampling.")

    output_workflow_json = args.output_workflow_json
    if not output_workflow_json:
        output_workflow_json = os.path.join(
            args.working_dir,
            f"{args.output_filename_prefix}{Path(args.input_workflow_json).name}"
        )

    initialize_downsamplers(args.working_dir)

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
