#!/usr/bin/env python
import shutil
import sys
from pathlib import Path
import tempfile
import argparse
import warnings

import attr
from tqdm.auto import tqdm as tqdm
import numpy
import pandas.core.arrays
import pysam
import torch
import concurrent.futures

from sv_utils import common, genomics_io
from gq_recalibrator import (
    dask_utils, training_utils, tarred_properties_to_parquet, vcf_filter_data_loader,
    gq_recalibrator_net, train_sv_gq_recalibrator
)
from gq_recalibrator.tarred_properties_to_parquet import PropertiesScaling, PropertiesSummary
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
from collections.abc import Iterator
from typing import Optional, TypeVar

tqdm.monitor_interval = 0

DaskDataFrame = dask.dataframe.DataFrame
PandasDataFrame = pandas.DataFrame
DataFrame = DaskDataFrame | PandasDataFrame
ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
PandasArray = pandas.core.arrays.PandasArray
ArrayType = TypeVar("ArrayType", numpy.ndarray, PandasArray, torch.Tensor)


class Keys:
    id = tarred_properties_to_parquet.Keys.id
    row = tarred_properties_to_parquet.Keys.row
    best_net_kwargs = train_sv_gq_recalibrator.Keys.best_net_kwargs
    net_kwargs = train_sv_gq_recalibrator.Keys.net_kwargs
    properties_scaling = train_sv_gq_recalibrator.Keys.properties_scaling
    properties_summary = train_sv_gq_recalibrator.Keys.properties_summary
    gq = genomics_io.Keys.gq
    ogq = genomics_io.Keys.ogq
    sl = genomics_io.Keys.sl
    allele_count = genomics_io.Keys.allele_count


class VcfKeys:
    gq = genomics_io.VcfKeys.gq
    ogq = genomics_io.VcfKeys.ogq
    sl = genomics_io.VcfKeys.sl


class Default:
    temp_dir = Path(tempfile.gettempdir())
    variants_per_batch = 1000
    non_filtration_properties = vcf_filter_data_loader.Default.non_filtration_properties
    torch_device_kind = training_utils.TorchDeviceKind.cpu
    scale_bool_properties = vcf_filter_data_loader.Default.scale_bool_properties
    percent_likelihood_per_logit = 0.1
    num_processes = vcf_filter_data_loader.Default.max_look_ahead_batches + 1
    encoding = "utf-8"
    original_gq_field = Keys.ogq
    scaled_logits_field = Keys.sl
    keep_multiallelic = True
    keep_homref = False
    keep_homvar = False
    compression_algorithm = tarred_properties_to_parquet.Default.compression_algorithm


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict likelihood each genotype is good using trained neural network and "
                    "extracted VCF properties. Output predicted VCF annotations to parquet",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--properties", "-p", type=Path, required=True,
                        help="path to tarred parquet file with variant properties")
    parser.add_argument("--model", "-m", type=Path, required=True,
                        help="path to input file with trained model")
    parser.add_argument("--output-parquet", "-o", type=Path, required=True,
                        help="path to output parquet files with GQs/SLs. If the path ends with"
                             ".tar, it will be tarred, otherwise, a folder with parquet files")
    parser.add_argument("--keep-multiallelic", type=bool, default=Default.keep_multiallelic,
                        help="If True, do not recalibrate multiallelic variants; if False, do so")
    parser.add_argument("--keep-homref", type=bool, default=Default.keep_homref,
                        help="If True, do not recalibrate HOMREF genotypes; if False, do so")
    parser.add_argument("--keep-homvar", type=bool, default=Default.keep_homvar,
                        help="If True, do not recalibrate HOMVAR genotypes; if False, do so")
    parser.add_argument("--variants-per-batch", type=int, default=Default.variants_per_batch,
                        help="number of variants used in each filter batch")
    parser.add_argument("--temp-dir", type=Path, default=Default.temp_dir,
                        help="preferred path to folder for temporary files")
    parser.add_argument(
        "--excluded-properties", type=str, default=','.join(Default.non_filtration_properties),
        help="comma-separated list of additional properties to not use in training or filtering. "
             f"Default excluded properties cannot be used (will remain excluded), passing"
             f"--exclude-properties will add one or more properties to exclude."
    )
    parser.add_argument("--original-gq-field", type=str, default=Default.original_gq_field,
                        help="Name for FORMAT field with original GQ values")
    parser.add_argument("--scaled-logits-field", type=str, default=Default.scaled_logits_field,
                        help="Name for FORMAT field with scaled logit values")
    # noinspection PyTypeChecker
    parser.add_argument(
        "--torch-device", type=str, default=Default.torch_device_kind,
        choices=training_utils.TorchDeviceKind.choices, help="Device on which to perform training"
    )
    parser.add_argument(
        "--scale-bool-properties", type=common.argparse_bool,
        default=Default.scale_bool_properties,
        help="If true, z-score bool variables too. If false, leave as 0.0 or 1.0. Scaling can "
             "result in large values if the bool variable is mostly true or mostly false."
    )
    parser.add_argument(
        "--num-processes", type=int, default=Default.num_processes,
        help="Number of cpu processes to use. One for filtering, the rest for loading data."
    )
    parser.add_argument(
        "--compression-algorithm", type=str, default=Default.compression_algorithm,
        choices=tarred_properties_to_parquet.CompressionAlgorithms.list(),
        help="compression algorithm for parquet data"
    )
    parser.add_argument(
        "--percent-likelihood-per-logit", type=float, default=Default.percent_likelihood_per_logit,
        help="Logit scale: the change in percent likelihood when going from 0 logits to 1"
    )
    parser.add_argument(
        "--use-best-net", action="store_true", default=False,
        help="If true, use the model rom the round with the best cross-validated loss. If false, "
             "use the final trained model. Because the cross-validated loss is somewhat noisy, the"
             "final model is probably the better option."
    )

    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    recalibrate_gq(
        properties_path=args.properties,
        model_path=args.model,
        output_parquet=args.output_parquet,
        temp_dir=args.temp_dir,
        keep_multiallelic=args.keep_multiallelic,
        keep_homref=args.keep_homref,
        keep_homvar=args.keep_homvar,
        variants_per_batch=args.variants_per_batch,
        excluded_properties=set(args.excluded_properties.split(',')),
        torch_device_kind=training_utils.TorchDeviceKind(args.torch_device),
        scale_bool_properties=args.scale_bool_properties,
        percent_likelihood_per_logit=args.percent_likelihood_per_logit,
        original_gq_field=args.original_gq_field,
        scaled_logits_field=args.scaled_logits_field,
        num_processes=args.num_processes,
        compression_algorithm=args.compression_algorithm,
        use_best_net=args.use_best_net,
    )


@attr.frozen(slots=True, weakref_slot=False)
class QualityCalculator:
    """Converts between probability, scaled-logits, and phred (GQ). Ensures appropriate device is
    used.
    Should be produced by the QualityCalculator.build method, which sets auxiliary values to
    handle round-off errors and integer min/max values correctly.
    """
    logit_scale: torch.Tensor
    min_sl_value: torch.Tensor
    max_sl_value: torch.Tensor
    min_sl_prob: torch.Tensor
    max_sl_prob: torch.Tensor
    phred_coef: torch.Tensor
    max_gq_value: torch.Tensor
    max_gq_prob: torch.Tensor
    out_device: torch.device
    quality_dtype: torch.dtype

    @staticmethod
    def build(
            in_device: torch.device,
            out_device: torch.device,
            percent_likelihood_per_logit: float = Default.percent_likelihood_per_logit,
            quality_dtype: torch.dtype = torch.short
    ) -> "QualityCalculator":
        """Produce a new QualityCalculator using converting between integer quality values and
        float32 torch tensors.
        """
        # compute logit scale
        logit_scale = QualityCalculator.calc_logit_scale(
            torch.tensor(percent_likelihood_per_logit, dtype=torch.float32, device=in_device)
        )
        # get boundaries to prevent overflow errors (torch will e.g. cast over-large float to
        # negative integers)
        iinfo = torch.iinfo(quality_dtype)
        i_min = torch.tensor(iinfo.min, dtype=quality_dtype, device=in_device)
        i_max = torch.tensor(iinfo.max, dtype=quality_dtype, device=in_device)
        finfo = torch.finfo(torch.float32)
        eps = torch.tensor(finfo.eps, dtype=torch.float32, device=in_device)
        min_sl_value = QualityCalculator._p_to_scaled_logits(
            predicted_probabilities=eps, logit_scale=logit_scale
        )
        if min_sl_value > i_min:
            min_sl_prob = eps
            min_sl_value = min_sl_value.to(dtype=quality_dtype)
        else:
            min_sl_value = i_min
            min_sl_prob = QualityCalculator._scaled_logits_to_p(
                scaled_logits=min_sl_value, logit_scale=logit_scale
            )
        max_sl_value = QualityCalculator._p_to_scaled_logits(
            predicted_probabilities=1.0 - eps, logit_scale=logit_scale
        )
        if max_sl_value < i_max:
            max_sl_prob = 1.0 - eps
            max_sl_value = max_sl_value.to(dtype=quality_dtype)
        else:
            max_sl_value = i_max
            max_sl_prob = QualityCalculator._scaled_logits_to_p(
                scaled_logits=max_sl_value, logit_scale=logit_scale
            )

        # compute phred_coef (only needed for phred_to_p)
        phred_coef = torch.tensor(numpy.log(10.0) / 10.0, dtype=torch.float32, device=in_device)
        # get boundaries to prevent overflow error
        max_gq_value = QualityCalculator._p_to_phred(predicted_probabilities=1.0 - eps)
        if max_gq_value < i_max:
            max_gq_prob = 1 - eps
            max_gq_value = max_gq_value.to(dtype=quality_dtype)
        else:
            max_gq_value = i_max
            max_gq_prob = QualityCalculator._phred_to_p(phred_coef=phred_coef, phred=max_sl_value)

        return QualityCalculator(
            logit_scale=logit_scale, min_sl_value=min_sl_value, max_sl_value=max_sl_value,
            min_sl_prob=min_sl_prob, max_sl_prob=max_sl_prob,
            phred_coef=phred_coef, max_gq_value=max_gq_value, max_gq_prob=max_gq_prob,
            out_device=out_device, quality_dtype=quality_dtype
        )

    @staticmethod
    def calc_logit_scale(percent_likelihood_per_logit: torch.Tensor) -> torch.Tensor:
        """Compute logit scale given the requested % change in likelihood / logit at p=0.5"""
        return 1.0 / torch.log(
            (0.5 + percent_likelihood_per_logit / 100) / (0.5 - percent_likelihood_per_logit / 100)
        )

    @staticmethod
    def _scaled_logits_to_p(
            scaled_logits: torch.Tensor,
            logit_scale: torch.Tensor
    ) -> torch.Tensor:
        """Static private method to convert scaled logits to probability"""
        return 1.0 / (1.0 + torch.exp(-scaled_logits / logit_scale))

    @staticmethod
    def _p_to_scaled_logits(
        predicted_probabilities: torch.Tensor,
        logit_scale: torch.Tensor
    ) -> torch.Tensor:
        """Static private method to convert probability to scaled logits"""
        return 1 + torch.floor(
            logit_scale * torch.log(
                predicted_probabilities / (1.0 - predicted_probabilities)
            )
        )

    def p_to_scaled_logits(
        self,
        probabilities: torch.Tensor,
        to_out_device: bool = True
    ) -> torch.Tensor:
        """Convert probability tensor to scaled logits, handling extreme values correctly

        Args:
            probabilities: tensor of probabilities
            to_out_device: if True, move tensor to output torch Device. If False, leave on input
                           torch Device.
        """
        scaled_logits = torch.where(
            probabilities <= self.min_sl_prob,
            self.min_sl_value,
            torch.where(
                probabilities >= self.max_sl_prob,
                self.max_sl_value,
                QualityCalculator._p_to_scaled_logits(
                    predicted_probabilities=probabilities, logit_scale=self.logit_scale,
                ).to(self.quality_dtype)
            )
        )
        return scaled_logits.to(device=self.out_device) if to_out_device else scaled_logits

    def get_output_scaled_logits(
        self,
        predicted_probabilities: torch.Tensor,
        is_filterable: torch.Tensor,
        original_gqs: torch.Tensor
    ) -> torch.Tensor:
        return self.p_to_scaled_logits(
            probabilities=torch.where(
                is_filterable,
                predicted_probabilities,
                self.phred_to_p(phred=original_gqs, to_out_device=False)
            )
        )

    @staticmethod
    def _phred_to_p(phred_coef: torch.Tensor, phred: torch.Tensor) -> torch.Tensor:
        return 1.0 - 0.5 * torch.exp(phred_coef * (1.0 - phred))

    def phred_to_p(self, phred: torch.Tensor, to_out_device: bool = True) -> torch.Tensor:
        p = QualityCalculator._phred_to_p(phred_coef=self.phred_coef, phred=phred)
        return p.to(self.out_device) if to_out_device else p

    @staticmethod
    def _p_to_phred(predicted_probabilities: torch.Tensor) -> torch.Tensor:
        return torch.floor(
            1.0 - 10 * torch.log10(2.0 * (1.0 - predicted_probabilities))
        )

    def p_to_phred(
        self,
        predicted_probabilities: torch.Tensor,
        to_out_device: bool = True
    ) -> torch.Tensor:
        phred = torch.where(
            predicted_probabilities > self.max_gq_prob,
            self.max_gq_value,
            QualityCalculator._p_to_phred(
                predicted_probabilities=predicted_probabilities
            ).to(dtype=self.quality_dtype)
        )
        return phred.to(self.out_device) if to_out_device else phred

    def get_output_gqs(
        self,
        predicted_probabilities: torch.Tensor,
        is_filterable: torch.Tensor,
        original_gqs: torch.Tensor
    ) -> torch.Tensor:
        return torch.where(
            is_filterable,
            self.p_to_phred(predicted_probabilities, to_out_device=False),
            original_gqs
        ).to(self.out_device)


def recalibrate_gq(
        properties_path: Path,
        model_path: Path,
        output_parquet: Path,
        temp_dir: Path = Default.temp_dir,
        keep_multiallelic: bool = Default.keep_multiallelic,
        keep_homref: bool = Default.keep_homref,
        keep_homvar: bool = Default.keep_homvar,
        variants_per_batch: int = Default.variants_per_batch,
        excluded_properties: set[str] = frozenset({}),
        torch_device_kind: training_utils.TorchDeviceKind = Default.torch_device_kind,
        scale_bool_properties: bool = Default.scale_bool_properties,
        percent_likelihood_per_logit: float = Default.percent_likelihood_per_logit,
        original_gq_field: str = Default.original_gq_field,
        scaled_logits_field: str = Default.scaled_logits_field,
        num_processes: int = Default.num_processes,
        compression_algorithm: str = Default.compression_algorithm,
        use_best_net: bool = False,
):
    non_filtration_properties = Default.non_filtration_properties.union(excluded_properties)
    num_batches = 0
    num_variants = 0
    filter_logger = training_utils.ProgressLogger()
    cpu_torch_device = torch_device_kind.cpu.get_device(progress_logger=None)

    if output_parquet.suffix == ".tar":
        output_tar = output_parquet
        parquet_folder = output_parquet.with_suffix("")
    else:
        output_tar = None
        parquet_folder = output_parquet
    shutil.rmtree(parquet_folder, ignore_errors=True)
    parquet_folder.mkdir(parents=True, exist_ok=False)

    with (
        torch.no_grad(),
        concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as process_executor
    ):
        recalibrator_net, properties_scaling, properties_summary \
            = load_recalibrator_and_data_transformer(
                model_path=model_path,
                use_best_net=use_best_net,
                torch_device=torch_device_kind.get_device(progress_logger=filter_logger)
            )
        vcf_tensor_data_loader = vcf_filter_data_loader.VcfFilterTensorDataLoader(
            parquet_path=properties_path,
            properties_scaling=properties_scaling,
            properties_summary=properties_summary,
            variants_per_batch=variants_per_batch,
            keep_multiallelic=keep_multiallelic,
            keep_homref=keep_homref,
            keep_homvar=keep_homvar,
            torch_device_kind=torch_device_kind,
            progress_logger=filter_logger,
            temp_dir=temp_dir,
            scale_bool_values=scale_bool_properties,
            non_filtration_properties=non_filtration_properties,
            process_executor=process_executor
        )

        quality_calculator = QualityCalculator.build(
            in_device=vcf_tensor_data_loader.torch_device,
            out_device=cpu_torch_device,
            percent_likelihood_per_logit=percent_likelihood_per_logit,
            quality_dtype=torch.int16
        )

        for batch_tensor, variant_ids, original_gqs, is_filterable in tqdm(
                vcf_tensor_data_loader, desc="batch", mininterval=0.5, maxinterval=float('inf'),
                smoothing=0
        ):
            predicted_probabilities = recalibrator_net(batch_tensor)
            output_batch(
                parquet_file=parquet_folder / f"part.{num_batches}.parquet",
                start_row=num_variants,
                sample_ids=vcf_tensor_data_loader.sample_ids,
                predicted_probabilities=predicted_probabilities,
                variant_ids=variant_ids,
                original_gqs=original_gqs,
                is_filterable=is_filterable,
                quality_calculator=quality_calculator,
                original_gq_field=original_gq_field,
                scaled_logits_field=scaled_logits_field,
                compression_algorithm=compression_algorithm
            )
            num_batches += 1
            num_variants += len(predicted_probabilities)

    if output_tar is not None:
        tarred_properties_to_parquet.archive_to_tar(
            folder_to_archive=parquet_folder, tar_file=output_tar, remove_files=True
        )

    print(f"Filtered {num_variants} in {num_batches} batches in {filter_logger.elapsed_time}")


def _get_output_header(
        input_header: pysam.VariantHeader,
        original_gq_field: str,
        scaled_logits_field: str
) -> pysam.VariantHeader:
    output_header = input_header
    output_header.formats.add(id=original_gq_field, number=1, type="Integer",
                              description="GQ values from original unrecalibrated VCF")
    output_header.formats.add(id=scaled_logits_field, number=1, type="Integer",
                              description="Scaled logits based on recalibrated probabilities")
    return output_header


def load_recalibrator_and_data_transformer(
        model_path: Path,
        use_best_net: bool,
        torch_device: torch.device
) -> tuple[gq_recalibrator_net.GqRecalibratorNet, PropertiesScaling, PropertiesSummary]:
    with open(model_path, "rb") as f_in:
        pickle_dict = torch.load(f_in)
    net_kwargs = Keys.best_net_kwargs if use_best_net else Keys.net_kwargs
    return (
        gq_recalibrator_net.GqRecalibratorNet(**pickle_dict[net_kwargs]).to(torch_device).eval(),
        pickle_dict[Keys.properties_scaling],
        pickle_dict[Keys.properties_summary],
    )


def output_batch(
        parquet_file: Path,
        start_row: int,
        sample_ids: list[str],
        predicted_probabilities: torch.Tensor,
        variant_ids: ArrowStringArray,
        original_gqs: torch.Tensor,
        is_filterable: torch.Tensor,
        quality_calculator: QualityCalculator,
        original_gq_field: str,
        scaled_logits_field: str,
        compression_algorithm: str = Default.compression_algorithm
):
    """Write parquet file for this batch with variant ID, allele_counts, GQ, original GQ,
    and scaled logits for each sample

    Args:
        parquet_file: Path to output parquet file for this batch
        start_row: first row of index for this dataframe
        sample_ids: list of sample IDs in order
        predicted_probabilities: tensor of predicted likelihood each genotype is good,
        variant_ids: array of variant IDs (one per row)
        original_gqs: tensor of original genotype qualities
        is_filterable: boolean tensor, True where genotype is filterable, False otherwise
        quality_calculator: object for translating between probabilities and quality scores
        original_gq_field: string with name of output field for original GQs
        scaled_logits_field: string with name of output field for scaled logits
        compression_algorithm: compression algorithm to use
    """
    gqs = quality_calculator.get_output_gqs(
        predicted_probabilities=predicted_probabilities,
        is_filterable=is_filterable,
        original_gqs=original_gqs
    )
    scaled_logits = quality_calculator.get_output_scaled_logits(
        predicted_probabilities=predicted_probabilities,
        is_filterable=is_filterable,
        original_gqs=original_gqs
    )
    out_device = quality_calculator.out_device
    original_gqs = original_gqs.detach().to(out_device)

    def _dict_items_iter() -> Iterator[tuple[Optional[str], str], ArrayType]:
        """Yield column_name, column pairs for constructing pandas DataFrame"""
        yield (None, Keys.id), variant_ids
        for column_ind, sample_id in enumerate(sample_ids):
            yield (sample_id, original_gq_field), original_gqs[:, column_ind]
            yield (sample_id, Keys.gq), gqs[:, column_ind]
            yield (sample_id, scaled_logits_field), scaled_logits[:, column_ind]

    stop_row = start_row + gqs.shape[0]
    df = pandas.DataFrame(
        {column_name: column_values for column_name, column_values in _dict_items_iter()},
        index=pandas.Index(numpy.arange(start_row, stop_row, dtype=numpy.int64), name=Keys.row)
    )
    dask_utils.flatten_columns(df)
    df.to_parquet(
        f"{parquet_file}", engine="pyarrow", compression=compression_algorithm, index=True
    )


if __name__ == "__main__":
    main()
