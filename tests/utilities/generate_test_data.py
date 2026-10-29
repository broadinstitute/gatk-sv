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
    # Set to True for inputs whose filenames are used by the tool (e.g., the genome tracks
    # of GQ recalibration, whose filenames are the names of the features of a trained model),
    # hence, the downsampled file is written with the same filename as the input.
    keep_filename: bool = False


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
    },
    # Inputs of sv-shell's `src/sv_shell/cnmops.sh`; see `gather_sample_evidence` for details.
    # Note that `bincov_matrix` and `ped_file` are created by upstream tasks (`make_bincov_matrix`
    # and `gather_batch_evidence`); `ped_file` is not downsampled, as it is not region-specific,
    # but you may use `--output-subdir` to set where it is copied.
    "cnmops": {
        "bincov_matrix": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"bincov_matrix": bed_gz, "bincov_matrix_index": bed_gz_tbi}
        ),
        "chrom_file": Handler(
            transformers.CnmopsContigsFaiDownsampler,
            lambda x: {"chrom_file": x}
        ),
        "allo_file": Handler(
            transformers.FaiDownsampler,
            lambda x: {"allo_file": x}
        ),
        "exclude_list": Handler(
            transformers.BedDownsampler,
            lambda x: {"exclude_list": x}
        ),
        "ref_dict": Handler(
            transformers.SequenceDictDownsampler,
            lambda x: {"ref_dict": x}
        ),
        "cytobands": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"cytobands": bed_gz}
        )
    },
    # Inputs of sv-shell's `src/sv_shell/cnv_germline_case.sh` (i.e., gCNV in case mode); see
    # `gather_sample_evidence` for details. Note that `counts` are created by an upstream task
    # (`condense_read_counts`). The gCNV model shards with no intervals within the target regions are dropped.
    "cnv_germline_case": {
        "counts": Handler(
            transformers.GcnvReadCountsDownsampler,
            lambda tsv_gz, tsv_gz_tbi: {"counts": tsv_gz}
        ),
        "contig_ploidy_model_tar": Handler(
            transformers.GcnvModelDownsampler,
            lambda x: {"contig_ploidy_model_tar": x}
        ),
        "gcnv_model_tars": Handler(
            transformers.GcnvModelShardDownsampler,
            lambda x: {"gcnv_model_tars": x}
        )
    },
    # Inputs of sv-shell's `src/sv_shell/cluster_batch.sh`; see `gather_sample_evidence` for details.
    # Note that the VCF tarballs and the DEL/DUP beds are created by upstream tasks
    # (in `single_sample_pipeline`, combining the outputs of `gather_batch_evidence` with the reference panel).
    "cluster_batch": {
        "reference_fasta": Handler(
            transformers.ReferenceDownsampler,
            lambda fasta, fai, dict, alt, amb, ann, bwt, pac, sa: {
                "reference_fasta": fasta, "reference_fasta_fai": fai, "reference_dict": dict}
        ),
        "contig_list": Handler(
            transformers.PrimaryContigsDownsampler,
            lambda x, _: {"contig_list": x}
        ),
        "manta_vcf_tar": Handler(
            transformers.SvVcfTarDownsampler,
            lambda x: {"manta_vcf_tar": x}
        ),
        "wham_vcf_tar": Handler(
            transformers.SvVcfTarDownsampler,
            lambda x: {"wham_vcf_tar": x}
        ),
        "scramble_vcf_tar": Handler(
            transformers.SvVcfTarDownsampler,
            lambda x: {"scramble_vcf_tar": x}
        ),
        "del_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"del_bed": bed_gz}
        ),
        "dup_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"dup_bed": bed_gz}
        ),
        "depth_exclude_intervals": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"depth_exclude_intervals": bed_gz}
        ),
        "pesr_exclude_intervals": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"pesr_exclude_intervals": bed_gz}
        )
    },
    # Inputs of sv-shell's `src/sv_shell/make_cohort_vcf.sh`; see `gather_sample_evidence` for details.
    # Note that the genotyped VCFs, PE evidence, binned coverage, and median coverage are created by upstream
    # tasks (`single_sample_pipeline` and `gather_batch_evidence`); the inputs without coordinates
    # (e.g., the median coverage, clustering configs, and genotyping tables) are copied.
    "make_cohort_vcf": {
        "reference_fasta": Handler(
            transformers.ReferenceDownsampler,
            lambda fasta, fai, dict, alt, amb, ann, bwt, pac, sa: {
                "reference_fasta": fasta, "reference_fasta_fai": fai, "reference_dict": dict}
        ),
        "pesr_vcfs": Handler(
            transformers.SvVcfDownsampler,
            lambda vcf_gz, vcf_gz_tbi: {"pesr_vcfs": vcf_gz}
        ),
        "depth_vcfs": Handler(
            transformers.SvVcfDownsampler,
            lambda vcf_gz, vcf_gz_tbi: {"depth_vcfs": vcf_gz}
        ),
        "disc_files": Handler(
            transformers.PeEvidenceDownsampler,
            lambda txt_gz, txt_gz_tbi: {"disc_files": txt_gz}
        ),
        "bincov_files": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"bincov_files": bed_gz}
        ),
        "contig_list": Handler(
            transformers.FaiDownsampler,
            lambda x: {"contig_list": x}
        ),
        "allosome_fai": Handler(
            transformers.FaiDownsampler,
            lambda x: {"allosome_fai": x}
        ),
        "primary_contigs_list": Handler(
            transformers.PrimaryContigsDownsampler,
            lambda x, _: {"primary_contigs_list": x}
        ),
        "cytobands": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"cytobands": bed_gz}
        ),
        "bin_exclude": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"bin_exclude": bed_gz}
        ),
        "mei_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"mei_bed": bed_gz}
        ),
        "pe_exclude_list": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"pe_exclude_list": bed_gz}
        ),
        "HERVK_reference": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"HERVK_reference": bed_gz}
        ),
        "LINE1_reference": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"LINE1_reference": bed_gz}
        ),
        "track_bed_files": Handler(
            transformers.BedDownsampler,
            lambda x: {"track_bed_files": x}
        ),
        "intron_reference": Handler(
            transformers.IntronTableDownsampler,
            lambda x: {"intron_reference": x}
        )
    },
    # Inputs of sv-shell's `src/sv_shell/single_sample_pipeline.sh` (non-DRAGEN mode); see `gather_sample_evidence`
    # for details. The lists of reference panel files (e.g., `ref_pesr_disc_files_list`) are downsampled by
    # downsampling each of the listed files. The inputs without coordinates (e.g., the reference panel's sample
    # list, PED file, metrics, cutoffs, genotyping tables, and models) are copied.
    "single_sample_pipeline": {
        "reference_fasta": Handler(
            transformers.ReferenceDownsampler,
            lambda fasta, fai, dict, alt, amb, ann, bwt, pac, sa: {
                "reference_fasta": fasta, "reference_index": fai, "reference_fasta_fai": fai, "reference_dict": dict,
                "reference_bwa_alt": alt, "reference_bwa_amb": amb, "reference_bwa_ann": ann,
                "reference_bwa_bwt": bwt, "reference_bwa_pac": pac, "reference_bwa_sa": sa}
        ),
        "bam_or_cram_file": Handler(
            transformers.CramDownsampler,
            lambda cram, index: {"bam_or_cram_file": cram, "bam_or_cram_index": index}
        ),
        "sample_bincov_matrix": Handler(
            transformers.IntervalListDownsampler,
            lambda tsv_gz, tsv_gz_tbi: {"sample_bincov_matrix": tsv_gz, "sample_bincov_matrix_idx": tsv_gz_tbi}
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
        "autosome_file": Handler(
            # In single_sample_pipeline, the autosomes are used only as the contigs that cn.mops runs on.
            transformers.CnmopsContigsFaiDownsampler,
            lambda x: {"autosome_file": x}
        ),
        "allosome_file": Handler(
            transformers.FaiDownsampler,
            lambda x: {"allosome_file": x}
        ),
        "genome_file": Handler(
            transformers.FaiDownsampler,
            lambda x: {"genome_file": x}
        ),
        "manta_region_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"manta_region_bed": bed_gz, "manta_region_bed_index": bed_gz_tbi}
        ),
        "mei_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"mei_bed": bed_gz}
        ),
        "wham_include_list_bed_file": Handler(
            transformers.BedDownsampler,
            lambda x: {"wham_include_list_bed_file": x}
        ),
        "cytobands": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"cytobands": bed_gz}
        ),
        "cnmops_exclude_list": Handler(
            transformers.BedDownsampler,
            lambda x: {"cnmops_exclude_list": x}
        ),
        "depth_exclude_list": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"depth_exclude_list": bed_gz}
        ),
        "pesr_exclude_intervals": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"pesr_exclude_intervals": bed_gz}
        ),
        "bin_exclude": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"bin_exclude": bed_gz}
        ),
        "HERVK_reference": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"HERVK_reference": bed_gz}
        ),
        "LINE1_reference": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"LINE1_reference": bed_gz}
        ),
        "wgd_scoring_mask": Handler(
            transformers.BedDownsampler,
            lambda x: {"wgd_scoring_mask": x}
        ),
        "rmsk": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"rmsk": bed_gz}
        ),
        "segdups": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"segdups": bed_gz}
        ),
        "noncoding_bed": Handler(
            transformers.BedDownsampler,
            lambda x: {"noncoding_bed": x}
        ),
        "par_bed": Handler(
            transformers.BedDownsampler,
            lambda x: {"par_bed": x}
        ),
        "external_af_ref_bed": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"external_af_ref_bed": bed_gz}
        ),
        "ref_panel_bincov_matrix": Handler(
            transformers.BincovMatrixDownsampler,
            lambda bed_gz, bed_gz_tbi: {"ref_panel_bincov_matrix": bed_gz}
        ),
        "ref_panel_del_bed": Handler(
            transformers.RefPanelDepthCallsDownsampler,
            lambda bed_gz, bed_gz_tbi: {"ref_panel_del_bed": bed_gz}
        ),
        "ref_panel_dup_bed": Handler(
            transformers.RefPanelDepthCallsDownsampler,
            lambda bed_gz, bed_gz_tbi: {"ref_panel_dup_bed": bed_gz}
        ),
        "clustering_track_bed_files": Handler(
            transformers.BedDownsampler,
            lambda x: {"clustering_track_bed_files": x}
        ),
        "genome_tracks": Handler(
            transformers.BedDownsampler,
            lambda bed_gz, bed_gz_tbi: {"genome_tracks": bed_gz},
            keep_filename=True
        ),
        "intron_reference": Handler(
            transformers.IntronTableDownsampler,
            lambda x: {"intron_reference": x}
        ),
        "protein_coding_gtf": Handler(
            transformers.GtfDownsampler,
            lambda x: {"protein_coding_gtf": x}
        ),
        "contig_ploidy_model_tar": Handler(
            transformers.GcnvModelDownsampler,
            lambda x: {"contig_ploidy_model_tar": x}
        ),
        "gcnv_model_tars_list": Handler(
            transformers.GcnvModelListDownsampler,
            lambda x: {"gcnv_model_tars_list": x}
        ),
        "ref_std_manta_vcf_tar": Handler(
            transformers.SvVcfTarDownsampler,
            lambda x: {"ref_std_manta_vcf_tar": x}
        ),
        "ref_std_wham_vcf_tar": Handler(
            transformers.SvVcfTarDownsampler,
            lambda x: {"ref_std_wham_vcf_tar": x}
        ),
        "ref_std_scramble_vcf_tar": Handler(
            transformers.SvVcfTarDownsampler,
            lambda x: {"ref_std_scramble_vcf_tar": x}
        ),
        "ref_pesr_disc_files_list": Handler(
            transformers.PeEvidenceListDownsampler,
            lambda x: {"ref_pesr_disc_files_list": x}
        ),
        "ref_pesr_split_files_list": Handler(
            transformers.PositionEvidenceListDownsampler,
            lambda x: {"ref_pesr_split_files_list": x}
        ),
        "ref_pesr_sd_files_list": Handler(
            transformers.PositionEvidenceListDownsampler,
            lambda x: {"ref_pesr_sd_files_list": x}
        ),
        "ref_panel_vcf": Handler(
            transformers.SvVcfDownsampler,
            lambda vcf_gz, vcf_gz_tbi: {"ref_panel_vcf": vcf_gz}
        ),
        "ref_samples_list": Handler(
            transformers.SampleListDownsampler,
            lambda x: {"ref_samples_list": x}
        ),
        "ref_ped_file": Handler(
            transformers.PedDownsampler,
            lambda x: {"ref_ped_file": x}
        ),
        "ref_panel_median_cov": Handler(
            transformers.MedianCovDownsampler,
            lambda x: {"ref_panel_median_cov": x}
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


# The inputs set to a fixed value in the downsampled workflow JSON, for the steps that are not tested.
WORKFLOW_INPUTS_TO_SET = {
    "single_sample_pipeline": {
        # STRipy is slow and is not tested on the downsampled data.
        "run_stripy": False
    }
}


# The workflows whose inputs are downsampled to a mini reference (see `transformers.MiniReference`),
# and their inputs that are the reference or its index files. When the reference is not downsampled
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
    ],
    "cnmops": [
        "ref_dict"
    ],
    "cnv_germline_case": [],
    "cluster_batch": [
        "reference_fasta",
        "reference_fasta_fai",
        "reference_dict"
    ],
    "make_cohort_vcf": [
        "reference_fasta",
        "reference_fasta_fai",
        "reference_dict"
    ],
    "single_sample_pipeline": [
        "reference_fasta",
        "reference_index",
        "reference_fasta_fai",
        "reference_dict",
        "reference_bwa_alt",
        "reference_bwa_amb",
        "reference_bwa_ann",
        "reference_bwa_bwt",
        "reference_bwa_pac",
        "reference_bwa_sa"
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


def parse_ref_panel_samples(input_filename: str = None) -> List[str]:
    if input_filename is None:
        return None
    with open(input_filename, "r") as f:
        samples = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    logging.info(f"Subsetting the reference panel to {len(samples)} samples.")
    return samples


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
        raise FileNotFoundError(f"Input file {input_filename} is neither a local file nor a gs:// path.")


def initialize_transformers(
        working_dir: str,
        reference_fasta: str,
        reference_index: str,
        sequence_dict_filename: str,
        picard_path: str,
        input_root: str = None,
        mini_reference: transformers.MiniReference = None,
        output_subdirs: dict = None,
        cnmops_contigs: List[str] = None,
        gcnv_contigs: List[str] = None,
        ref_panel_samples: List[str] = None):
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
                output_subdir=output_subdirs.get(input_var),
                cnmops_contigs=cnmops_contigs,
                gcnv_contigs=gcnv_contigs,
                ref_panel_samples=ref_panel_samples
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


def load_workflow_inputs(input_filename: str) -> dict:
    """
    Loads the given workflow JSON; the relative paths of local files in it
    are resolved relative to the JSON file's directory.
    """
    with open(input_filename, "r") as f:
        workflow_inputs = json.load(f)
    base_dir = Path(input_filename).resolve().parent

    def resolve(v):
        if isinstance(v, list):
            return [resolve(x) for x in v]
        if isinstance(v, str) and not os.path.isabs(v) and base_dir.joinpath(v).is_file():
            return str(base_dir.joinpath(v))
        return v
    return {k: resolve(v) for k, v in workflow_inputs.items()}


def relocate_output_paths(output_filename: str, working_dir: str, output_root: str):
    """
    Replaces the working directory with the given root (e.g., `/test_inputs`, where the working directory is
    mounted in a container) in the paths of the given output JSON, and of the list files it references.
    """
    working_dir = str(Path(working_dir).resolve())
    relocate = lambda x: x.replace(working_dir, output_root.rstrip("/"))
    with open(output_filename, "r") as f:
        workflow_inputs = json.load(f)
    for v in workflow_inputs.values():
        for x in (v if isinstance(v, list) else [v]):
            if isinstance(x, str) and x.startswith(working_dir) and x.endswith(".list") and os.path.isfile(x):
                with open(x, "r") as f:
                    content = f.read()
                with open(x, "w") as f:
                    f.write(relocate(content))
    workflow_inputs = {
        k: [relocate(x) if isinstance(x, str) else x for x in v] if isinstance(v, list)
        else relocate(v) if isinstance(v, str) else v
        for k, v in workflow_inputs.items()}
    logging.info(f"Relocating the paths in {output_filename} from {working_dir} to {output_root}.")
    with open(output_filename, "w") as f:
        json.dump(workflow_inputs, f, indent=4)


def get_local_files(workflow_inputs: dict) -> List[str]:
    values = []
    for v in workflow_inputs.values():
        values.extend(v if isinstance(v, list) else [v])
    return [v for v in values if isinstance(v, str) and os.path.isfile(v)]


def uses_mini_reference(workflow_inputs: dict, workflow_name: str = None) -> bool:
    """
    Returns True if the inputs belong to a workflow listed in `REFERENCE_INPUTS`,
    hence, the reference should be downsampled to a mini reference.
    """
    for k in workflow_inputs.keys():
        key_workflow_name, _, _ = split_input_key(k, workflow_name)
        if key_workflow_name in REFERENCE_INPUTS:
            return True
    return False


def parse_output_subdir(value: str):
    input_var, subdir = value.split("=", maxsplit=1)
    return input_var, subdir


def copy_local_file(working_dir: str, input_root: str, filename: str, output_subdir: str = None) -> str:
    """
    Copies the given file to the given subdirectory of the working directory if provided, otherwise,
    following the same directory hierarchy as the file has relative to the `input_root` directory.
    """
    if output_subdir is not None:
        relative_path = Path(output_subdir).joinpath(Path(filename).name)
    elif Path(filename).resolve().is_relative_to(Path(input_root).resolve()):
        relative_path = Path(filename).resolve().relative_to(Path(input_root).resolve())
    else:
        relative_path = Path(filename).name
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
        workflow_name: str = None, input_root: str = None, skip_reference_downsampling: bool = False,
        output_subdirs: dict = None):
    output_subdirs = output_subdirs or {}
    workflow_inputs = load_workflow_inputs(input_filename)

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
            # Inputs downsampled with their filename unchanged are localized in a subdirectory,
            # so the downsampled file does not overwrite its input.
            localize_dir = Path(working_dir).joinpath("localized") if handler.keep_filename else Path(working_dir)
            localize_dir.mkdir(parents=True, exist_ok=True)
            workflow_input_local_filename = localize_file(x, localize_dir.joinpath(Path(x).name))
            updated_files.append(handler.transformer.transform(
                input_filename=workflow_input_local_filename,
                output_prefix="" if handler.keep_filename else output_filename_prefix,
                regions=regions
            ))

        # A transformer returns None if the file has no data within the target regions,
        # (e.g., a gCNV model shard), in which case, the file is dropped from the list.
        updated_files = [x for x in updated_files if x is not None]
        if len(updated_files) == 0:
            if not isinstance(v, list):
                raise ValueError(f"Input {k} has no data within the target regions.")
            logging.warning(f"None of the files of input {k} has data within the target regions.")
            updated_workflow_inputs[k] = []
            continue

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
            updated_workflow_inputs[k] = copy_local_file(working_dir, input_root, v, output_subdirs.get(input_var))
        elif input_root is not None and isinstance(v, list) and len(v) > 0 \
                and all(isinstance(x, str) and os.path.isfile(x) for x in v):
            updated_workflow_inputs[k] = [
                copy_local_file(working_dir, input_root, x, output_subdirs.get(input_var)) for x in v]
        else:
            updated_workflow_inputs[k] = v
            logging.info(f"Leaving {k} unchanged.")

    for k, v in WORKFLOW_INPUTS_TO_SET.get(workflow_name, {}).items():
        logging.info(f"Setting {k} to {v}.")
        updated_workflow_inputs[k] = v

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
        "--output-root-in-json",
        help="Sets the path that replaces the working directory in the paths of the output JSON and of the "
             "list files it references; e.g., `/test_inputs` when the working directory is mounted there "
             "in the container that runs the workflow. The default is the working directory's path."
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

    parser.add_argument(
        "--cnmops-contigs",
        nargs="+",
        help="Sets the contigs that cn.mops runs on (i.e., the contigs kept in the `chrom_file` of cnmops and "
             "the `autosome_file` of single_sample_pipeline), e.g., `chr16 chr20`. Since cn.mops runs on each "
             "contig separately, fewer contigs shorten its runtime. The default is all the contigs of the "
             "target regions in these files."
    )

    parser.add_argument(
        "--gcnv-contigs",
        nargs="+",
        help="Sets the contigs that gCNV runs on, e.g., `chr16 chr20 chrX chrY`; i.e., only the intervals "
             "on these contigs are kept in the gCNV model shards (`gcnv_model_tars` of cnv_germline_case "
             "and `gcnv_model_tars_list` of single_sample_pipeline), and the shards with no such intervals "
             "are dropped. Since gCNV runs on each model shard separately, fewer shards shorten its runtime. "
             "The contig-ploidy model is not affected. The default is all the contigs of the target regions."
    )

    parser.add_argument(
        "--ref-panel-samples",
        help="Sets a file with the IDs of the reference panel samples to keep, one per line (e.g., "
             "`sv_shell_ref_panel_samples.list`); i.e., the reference panel inputs of single_sample_pipeline "
             "(the samples list, PED file, evidence files, standardized VCFs, bin coverage matrix, median "
             "coverages, and depth calls) are subset to these samples. Since the batch steps (e.g., cn.mops "
             "and clustering) scale with the number of samples, fewer samples shorten their runtime. "
             "The default is all the samples."
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

    workflow_inputs = load_workflow_inputs(args.input_workflow_json)

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
                            dict(args.output_subdir), args.cnmops_contigs,
                            args.gcnv_contigs, parse_ref_panel_samples(args.ref_panel_samples))

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
        skip_reference_downsampling=args.skip_reference_downsampling,
        output_subdirs=dict(args.output_subdir))

    if args.output_root_in_json is not None:
        relocate_output_paths(output_workflow_json, args.working_dir, args.output_root_in_json)

    if mini_reference is not None:
        mini_reference.cleanup()

    logging.info("All process finished successfully.")


if __name__ == '__main__':
    main()
