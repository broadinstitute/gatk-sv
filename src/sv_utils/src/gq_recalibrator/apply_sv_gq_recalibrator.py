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

from sv_utils import common
from gq_recalibrator import (
    tarred_properties_to_parquet, vcf_filter_data_loader, training_utils, gq_recalibrator_net
)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe

from types import MappingProxyType
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


class Default:
    temp_dir = Path(tempfile.gettempdir())
    variants_per_batch = 100
    batches_per_mini_epoch = 50
    shuffle = vcf_filter_data_loader.Default.shuffle
    random_state = vcf_filter_data_loader.Default.random_seed
    excluded_properties = vcf_filter_data_loader.Default.excluded_properties
    torch_device_kind = training_utils.TorchDeviceKind.cpu
    scale_bool_properties = vcf_filter_data_loader.Default.scale_bool_properties
    percent_likelihood_per_logit = 0.1
    validation_proportion = vcf_filter_data_loader.Default.validation_proportion
    num_hidden_layers = gq_recalibrator_net.Default.num_hidden_layers
    layer_expansion_factor = gq_recalibrator_net.Default.layer_expansion_factor
    bias = gq_recalibrator_net.Default.bias
    optim_type = torch.optim.AdamW
    optim_kwargs = MappingProxyType({"weight_decay": 0.1, "lr": 1.0e-3})
    hidden_nonlinearity = gq_recalibrator_net.Default.hidden_nonlinearity
    output_nonlinearity = gq_recalibrator_net.Default.output_nonlinearity
    max_gradient_value = 0.0  # 0.01
    max_gradient_norm = 0.0  # 1.0e-3
    correlation_loss_weight = 0.1
    num_processes = vcf_filter_data_loader.Default.max_look_ahead_batches + 1
    encoding = "utf-8"


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict likelihood each genotype is good using trained neural network and "
                    "extracted VCF properties",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--properties", "-p", type=Path, required=True,
                        help="path to tarred parquet file with variant properties")
    parser.add_argument("--properties-scaling-json", "-j", type=Path, required=True,
                        help="path to JSON with baseline and scales needed for training / "
                             "filtering properties")
    parser.add_argument("--model", "-m", type=Path, required=False,
                        help="path to input file with trained model")
    parser.add_argument("--output", "-o", type=Path, required=True,
                        help="path to output predicted probabilities TSV")
    parser.add_argument("--variants-per-batch", type=int, default=Default.variants_per_batch,
                        help="number of variants used in each filter batch")
    parser.add_argument("--temp-dir", type=Path, default=Default.temp_dir,
                        help="preferred path to folder for temporary files")
    parser.add_argument(
        "--excluded-properties", type=str, default=','.join(Default.excluded_properties),
        help="comma-separated list of properties to not use in tensor"
    )
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
        output_path=args.output,
        variants_per_batch=args.variants_per_batch,
        excluded_properties=set(args.excluded_properties.split(',')),
        torch_device_kind=training_utils.TorchDeviceKind(args.torch_device),
        num_processes=args.num_processes
    )


@attr.frozen
class QualityCalculator:
    min_value: torch.Tensor
    max_value: torch.Tensor
    logit_scale: torch.Tensor
    min_sl_prob: torch.Tensor
    max_sl_prob: torch.Tensor
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
        iinfo = torch.iinfo(quality_dtype)
        min_value = torch.tensor(iinfo.min, dtype=quality_dtype, device=in_device)
        max_value = torch.tensor(iinfo.max, dtype=quality_dtype, device=in_device)
        logit_scale = QualityCalculator.calc_logit_scale(
            torch.tensor(percent_likelihood_per_logit, dtype=torch.float32, device=in_device)
        )
        min_sl_prob = QualityCalculator._scaled_logits_to_p(
            scaled_logits=min_value, logit_scale=logit_scale
        )
        max_sl_prob = QualityCalculator._scaled_logits_to_p(
            scaled_logits=max_value, logit_scale=logit_scale
        )
        max_gq_prob = QualityCalculator._phred_to_p(phred=max_value)

        return QualityCalculator(
            min_value=min_value, max_value=max_value, logit_scale=logit_scale,
            min_sl_prob=min_sl_prob, max_sl_prob=max_sl_prob, max_gq_prob=max_gq_prob,
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

    def p_to_scaled_logits(self, predicted_probabilities: torch.Tensor) -> torch.Tensor:
        return torch.where(
            predicted_probabilities < self.min_sl_prob,
            self.min_value,
            torch.where(
                predicted_probabilities > self.max_sl_prob,
                self.max_value,
                torch.floor(
                    self.logit_scale * torch.log(
                        predicted_probabilities / (1.0 - predicted_probabilities)
                    )
                ).to(dtype=self.quality_dtype)
            )
        ).to(device=self.out_device)

    @staticmethod
    def _phred_to_p(phred: torch.Tensor) -> torch.Tensor:
        phred_coef = torch.tensor(numpy.log(10.0) / 10.0, dtype=torch.float32, device=phred.device)
        return 1.0 - 0.5 * torch.exp(phred_coef * (1.0 - phred))

    def p_to_phred(self, predicted_probabilities: torch.Tensor) -> torch.Tensor:
        return torch.where(
            predicted_probabilities > self.max_gq_prob,
            self.max_value,
            torch.floor(
                self.logit_scale * torch.log(
                    predicted_probabilities / (1.0 - predicted_probabilities)
                )
            ).to(dtype=self.quality_dtype)
        ).to(device=self.out_device)


def apply_sv_gq_recalibrator(
        properties_path: Path,
        properties_scaling_json: Path,
        model_path: Path,
        output_path: Path,
        variants_per_batch: int = Default.variants_per_batch,
        excluded_properties: set[str] = Default.excluded_properties,
        torch_device_kind: training_utils.TorchDeviceKind = Default.torch_device_kind,
        scale_bool_values: bool = Default.scale_bool_properties,
        percent_likelihood_per_logit: float = Default.percent_likelihood_per_logit,
        num_processes: int = Default.num_processes,
        encoding: str = Default.encoding
):
    num_batches = 0
    num_variants = 0
    filter_logger = training_utils.ProgressLogger()
    cpu_torch_device = torch_device_kind.cpu.get_device(progress_logger=None)

    with (
        concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as process_executor,
        pysam.BGZFile(f"{output_path}", mode="wb", index=None) as f_out,
        torch.no_grad()
    ):
        vcf_tensor_data_loader = vcf_filter_data_loader.VcfFilterTensorDataLoader(
            parquet_path=properties_path, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values, variants_per_batch=variants_per_batch,
            excluded_properties=excluded_properties, torch_device_kind=torch_device_kind,
            progress_logger=filter_logger, process_executor=process_executor
        )
        gq_recalibrator = load_gq_recalibrator(model_path=model_path,
                                               torch_device=vcf_tensor_data_loader.torch_device)
        quality_calculator = QualityCalculator.make(
            in_device=vcf_tensor_data_loader.torch_device,
            out_device=cpu_torch_device,
            percent_likelihood_per_logit=percent_likelihood_per_logit,
            quality_dtype=torch.short
        )

        for batch_tensor in tqdm(
                vcf_tensor_data_loader, desc="batch", mininterval=0.5, maxinterval=float('inf'),
                smoothing=0
        ):
            predicted_probabilities = gq_recalibrator(batch_tensor)
            output_batch(
                f_out=f_out,
                predicted_probabilities=predicted_probabilities,
                vcf_tensor_data_loader=vcf_tensor_data_loader,
                write_header=(num_batches == 0),
                quality_calculator=quality_calculator,
                encoding=encoding
            )
            num_batches += 1
            num_variants += len(predicted_probabilities)
    print(f"Filtered {num_variants} in {num_batches} batches in {filter_logger.elapsed_time}")


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
        f_out: pysam.BGZFile,
        predicted_probabilities: torch.tensor,
        vcf_tensor_data_loader: vcf_filter_data_loader.VcfFilterTensorDataLoader,
        write_header: bool,
        quality_calculator: QualityCalculator,
        encoding: str = Default.encoding,
):
    last_sample_ind = vcf_tensor_data_loader.num_samples - 1
    if write_header:
        for ind, sample_id in enumerate(vcf_tensor_data_loader.sample_ids):
            last_space = '\n' if ind == last_sample_ind else '\t'
            f_out.write(
                f"{sample_id},probability\t{sample_id},GQ\t{sample_id},SL{last_space}"
                .encode(encoding)
            )
        f_out.write("\n".encode(encoding))
    gqs = quality_calculator.p_to_phred(predicted_probabilities)
    sls = quality_calculator.p_to_scaled_logits(predicted_probabilities)
    out_device = quality_calculator.out_device
    for ind, (probability, gq, sl) in enumerate(
            zip(predicted_probabilities.to(out_device), gqs, sls)
    ):
        last_space = '\n' if ind == last_sample_ind else '\t'
        f_out.write(f"{probability}\t{gq}\t{sl}{last_space}".encode(encoding))


if __name__ == "__main__":
    main()
