#!/usr/bin/env python
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
    tarred_properties_to_parquet, vcf_filter_data_loader, training_utils, gq_recalibrator_net
)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
from typing import Optional, Union

tqdm.monitor_interval = 0

DaskDataFrame = dask.dataframe.DataFrame
PandasDataFrame = pandas.DataFrame
DataFrame = Union[DaskDataFrame, PandasDataFrame]
ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray


class Keys:
    id = tarred_properties_to_parquet.Keys.id
    variant_weights = tarred_properties_to_parquet.Keys.variant_weights
    positive_value = tarred_properties_to_parquet.Keys.positive_value
    negative_value = tarred_properties_to_parquet.Keys.negative_value
    baseline = tarred_properties_to_parquet.Keys.baseline
    scale = tarred_properties_to_parquet.Keys.scale
    logger_kwargs = "logger_kwargs"
    data_loader_state = "data_loader_state"
    net_kwargs = "neural_net_kwargs"
    optimizer_kwargs = "optimizer_kwargs"
    optimizer_type = "optimizer_type"
    training_losses = "training_losses"
    validation_losses = "validation_losses"


class VcfKeys:
    gq = genomics_io.VcfKeys.gq


class Default:
    temp_dir = Path(tempfile.gettempdir())
    variants_per_batch = 100
    batches_per_mini_epoch = 50
    non_filtration_properties = vcf_filter_data_loader.Default.non_filtration_properties
    torch_device_kind = training_utils.TorchDeviceKind.cpu
    scale_bool_properties = vcf_filter_data_loader.Default.scale_bool_properties
    percent_likelihood_per_logit = 0.1
    num_processes = vcf_filter_data_loader.Default.max_look_ahead_batches + 1
    encoding = "utf-8"
    original_gq_field = "OGQ"
    scaled_logits_field = "SL"
    keep_multiallelic = True
    keep_homref = False
    keep_homvar = False


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict likelihood each genotype is good using trained neural network and "
                    "extracted VCF properties and update VCF with recalibrated qualities",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--properties", "-p", type=Path, required=True,
                        help="path to tarred parquet file with variant properties")
    parser.add_argument("--properties-scaling-json", "-j", type=Path, required=True,
                        help="path to JSON with baseline and scales needed for training / "
                             "filtering properties")
    parser.add_argument("--model", "-m", type=Path, required=True,
                        help="path to input file with trained model")
    parser.add_argument("--input-vcf", "-i", type=Path, required=True,
                        help="path to input VCF to be recalibrated")
    parser.add_argument("--output-vcf", "-o", type=Path, required=True,
                        help="path to output VCF with recalibrated GQs/SLs")
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
        help="Number of cpu processes to use. One for training, the rest for loading data."
    )
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)

    apply_sv_gq_recalibrator(
        properties_path=args.properties,
        properties_scaling_json=args.properties_scaling_json,
        model_path=args.model,
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        variants_per_batch=args.variants_per_batch,
        excluded_properties=set(args.excluded_properties.split(',')),
        torch_device_kind=training_utils.TorchDeviceKind(args.torch_device),
        num_processes=args.num_processes
    )


@attr.frozen
class QualityCalculator:
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
    def make(
            in_device: torch.device,
            out_device: torch.device,
            percent_likelihood_per_logit: float = Default.percent_likelihood_per_logit,
            quality_dtype: torch.dtype = torch.short
    ) -> "QualityCalculator":
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
        if min_sl_value >= i_min:
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
        if max_sl_value <= i_max:
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
        return 1.0 / torch.log(
            (0.5 + percent_likelihood_per_logit / 100) / (0.5 - percent_likelihood_per_logit / 100)
        )

    @staticmethod
    def _scaled_logits_to_p(
            scaled_logits: torch.Tensor,
            logit_scale: torch.Tensor
    ) -> torch.Tensor:
        return 1.0 / (1.0 + torch.exp(-scaled_logits / logit_scale))

    @staticmethod
    def _p_to_scaled_logits(
        predicted_probabilities: torch.Tensor,
        logit_scale: torch.Tensor
    ) -> torch.Tensor:
        return 1 + torch.floor(
            logit_scale * torch.log(
                predicted_probabilities / (1.0 - predicted_probabilities)
            )
        )

    def p_to_scaled_logits(
        self,
        predicted_probabilities: torch.Tensor,
        to_out_device: bool = True
    ) -> torch.Tensor:
        scaled_logits = torch.where(
            predicted_probabilities <= self.min_sl_prob,
            self.min_sl_value,
            torch.where(
                predicted_probabilities >= self.max_sl_prob,
                self.max_sl_value,
                QualityCalculator._p_to_scaled_logits(
                    predicted_probabilities=predicted_probabilities, logit_scale=self.logit_scale,
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
            predicted_probabilities=torch.where(
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


def apply_sv_gq_recalibrator(
        properties_path: Path,
        properties_scaling_json: Path,
        model_path: Path,
        input_vcf: Path,
        output_vcf: Path,
        keep_multiallelic: bool = Default.keep_multiallelic,
        keep_homref: bool = Default.keep_homref,
        keep_homvar: bool = Default.keep_homvar,
        variants_per_batch: int = Default.variants_per_batch,
        excluded_properties: set[str] = frozenset({}),
        torch_device_kind: training_utils.TorchDeviceKind = Default.torch_device_kind,
        scale_bool_values: bool = Default.scale_bool_properties,
        percent_likelihood_per_logit: float = Default.percent_likelihood_per_logit,
        original_gq_field: str = Default.original_gq_field,
        scaled_logits_field: str = Default.scaled_logits_field,
        num_processes: int = Default.num_processes
):
    non_filtration_properties = Default.non_filtration_properties.union(excluded_properties)
    num_batches = 0
    num_variants = 0
    filter_logger = training_utils.ProgressLogger()
    cpu_torch_device = torch_device_kind.cpu.get_device(progress_logger=None)

    with (
        torch.no_grad(),
        concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as process_executor,
        pysam.VariantFile(f"{input_vcf}", mode='r', threads=2) as vcf_in,
        pysam.VariantFile(
            f"{output_vcf}",
            mode='w',
            threads=2,
            header=_get_output_header(
                input_header=vcf_in.header,
                original_gq_field=original_gq_field,
                scaled_logits_field=scaled_logits_field
            )
        ) as vcf_out
    ):
        vcf_tensor_data_loader = vcf_filter_data_loader.VcfFilterTensorDataLoader(
            parquet_path=properties_path,
            properties_scaling_json=properties_scaling_json,
            variants_per_batch=variants_per_batch,
            keep_multiallelic=keep_multiallelic,
            keep_homref=keep_homref,
            keep_homvar=keep_homvar,
            torch_device_kind=torch_device_kind,
            progress_logger=filter_logger,
            scale_bool_values=scale_bool_values,
            non_filtration_properties=non_filtration_properties,
            process_executor=process_executor
        )
        gq_recalibrator = load_gq_recalibrator(model_path=model_path,
                                               torch_device=vcf_tensor_data_loader.torch_device)
        quality_calculator = QualityCalculator.make(
            in_device=vcf_tensor_data_loader.torch_device,
            out_device=cpu_torch_device,
            percent_likelihood_per_logit=percent_likelihood_per_logit,
            quality_dtype=torch.short
        )

        for batch_tensor, original_gqs, is_filterable in tqdm(
                vcf_tensor_data_loader, desc="batch", mininterval=0.5, maxinterval=float('inf'),
                smoothing=0
        ):
            predicted_probabilities = gq_recalibrator(batch_tensor)
            output_batch(
                vcf_in=vcf_in,
                vcf_out=vcf_out,
                predicted_probabilities=predicted_probabilities,
                original_gqs=original_gqs,
                is_filterable=is_filterable,
                quality_calculator=quality_calculator,
                original_gq_field=original_gq_field,
                scaled_logits_field=scaled_logits_field
            )
            num_batches += 1
            num_variants += len(predicted_probabilities)
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


def load_gq_recalibrator(
        model_path: Path,
        torch_device: torch.device
) -> gq_recalibrator_net.GqRecalibratorNet:
    with open(model_path, "rb") as f_in:
        pickle_dict = torch.load(f_in)
    return gq_recalibrator_net.GqRecalibratorNet(
        **pickle_dict[Keys.net_kwargs]
    ).to(torch_device).eval()


def output_batch(
        vcf_in: pysam.VariantFile,
        vcf_out: pysam.VariantFile,
        predicted_probabilities: torch.Tensor,
        original_gqs: torch.Tensor,
        is_filterable: torch.Tensor,
        quality_calculator: QualityCalculator,
        original_gq_field: str,
        scaled_logits_field: str
):
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
    num_rows = gqs.shape[0]
    for batch_row in range(num_rows):
        in_record = next(vcf_in)
        out_record = in_record.copy()
        for column_ind, out_gt in enumerate(out_record.samples.itervalues()):
            out_gt.update({
                VcfKeys.gq: int(gqs[batch_row, column_ind]),
                original_gq_field: int(original_gqs[batch_row, column_ind]),
                scaled_logits_field: int(scaled_logits[batch_row, column_ind])
            })
        vcf_out.write(out_record)


if __name__ == "__main__":
    main()
