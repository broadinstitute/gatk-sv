#!/usr/bin/env python

import sys
import os
from pathlib import Path
import tempfile
import argparse
import warnings
import time
from tqdm.auto import tqdm as tqdm
import numpy.random
import pandas.core.arrays
import torch
import concurrent.futures

from sv_utils import common
from gq_recalibrator import tarred_properties_to_parquet, vcf_tensor_data_loaders, training_utils, gq_recalibrator_net
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
from types import MappingProxyType
from typing import Optional, Union, Any
from collections.abc import Mapping

tqdm.monitor_interval = 0

DaskDataFrame = dask.dataframe.DataFrame
PandasDataFrame = pandas.DataFrame
DataFrame = Union[DaskDataFrame, PandasDataFrame]
ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray


class Keys:
    cuda = vcf_tensor_data_loaders.Keys.cuda
    cpu = vcf_tensor_data_loaders.Keys.cpu
    mps = vcf_tensor_data_loaders.Keys.mps
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
    shuffle = vcf_tensor_data_loaders.Default.shuffle
    random_state = vcf_tensor_data_loaders.Default.random_seed
    excluded_properties = vcf_tensor_data_loaders.Default.excluded_properties
    torch_device = Keys.cpu
    scale_bool_properties = vcf_tensor_data_loaders.Default.scale_bool_properties
    validation_proportion = vcf_tensor_data_loaders.Default.validation_proportion
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
    num_processes = vcf_tensor_data_loaders.Default.max_look_ahead_batches + 1


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train filter to recalibrate genotype quality from extracted VCF properties",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--properties", "-p", type=Path, required=True,
                        help="full path to tarred parquet file with variant properties")
    parser.add_argument("--properties-scaling-json", "-j", type=Path, required=True,
                        help="full path to JSON with baseline and scales needed for training / "
                             "filtering properties")
    parser.add_argument(
        "--truth-json", type=Path, required=True,
        help="full path to JSON with info about which GTs are good for each variant x sample"
    )
    parser.add_argument("--output-model", "-o", type=Path, required=True,
                        help="path to output file with trained model")
    parser.add_argument("--input-model", "-i", type=Path, required=False,
                        help="path to input file with partially trained model")
    parser.add_argument("--variants-per-batch", type=int, default=Default.variants_per_batch,
                        help="number of variants used in each training batch")
    parser.add_argument("--batches-per-mini-epoch", type=int,
                        default=Default.batches_per_mini_epoch,
                        help="number of batches between progress updates and checkpoints")
    parser.add_argument("--max-train-variants", type=int, default=None,
                        help="maximum number of variants to train before ending training")
    parser.add_argument("--temp-dir", type=Path, default=Default.temp_dir,
                        help="preferred path to folder for temporary files")
    parser.add_argument("--shuffle", type=common.argparse_bool, default=Default.shuffle,
                        help="if True, shuffle order of variants while training")
    parser.add_argument("--random-seed", type=int, default=Default.random_state,
                        help="initial random seed")
    parser.add_argument(
        "--excluded-properties", type=str, default=','.join(Default.excluded_properties),
        help="comma-separated list of properties to not use in tensor"
    )
    parser.add_argument(
        "--torch-device", type=str, default=Default.torch_device,
        choices=[Keys.cpu, Keys.cuda, Keys.mps], help="Device on which to perform training"
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
    train_sv_gq_recalibrator(
        properties_path=args.properties,
        properties_scaling_json=args.properties_scaling_json,
        truth_json=args.truth_json,
        output_model_path=args.output_model,
        input_model_path=args.input_model,
        variants_per_batch=args.variants_per_batch,
        batches_per_mini_epoch=args.batches_per_mini_epoch,
        max_train_variants=args.max_train_variants,
        shuffle=args.shuffle,
        random_state=args.random_seed,
        excluded_properties=set(args.excluded_properties.split(',')),
        torch_device=args.torch_device,
        num_cpu_processes=args.num_processes
    )
    # filter_sv_gq_recalibrator(
    #     properties_path=args.properties,
    #     properties_scaling_json=args.properties_scaling_json,
    #     output_model_path=args.output_model,
    #     variants_per_batch=args.variants_per_batch,
    #     shuffle=args.shuffle,
    #     random_state=args.random_seed,
    #     excluded_properties=set(args.excluded_properties.split(',')),
    #     torch_device=args.torch_device
    # )


def train_sv_gq_recalibrator(
        properties_path: Path,
        properties_scaling_json: Path,
        truth_json: Path,
        output_model_path: Path,
        input_model_path: Optional[Path] = None,
        variants_per_batch: int = Default.variants_per_batch,
        batches_per_mini_epoch: int = Default.batches_per_mini_epoch,
        max_train_variants: Optional[int] = None,
        shuffle: bool = Default.shuffle,
        random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
        excluded_properties: set[str] = Default.excluded_properties,
        torch_device: str = Default.torch_device,
        scale_bool_values: bool = Default.scale_bool_properties,
        num_hidden_layers: int = Default.num_hidden_layers,
        layer_expansion_factor: float = Default.layer_expansion_factor,
        bias: bool = Default.bias,
        optim_type: type = Default.optim_type,
        optim_kwargs: Mapping[str, Any] = Default.optim_kwargs,
        hidden_nonlinearity: torch.nn.Module = Default.hidden_nonlinearity,
        output_nonlinearity: torch.nn.Module = Default.output_nonlinearity,
        model_state_dict: Optional[dict[str, Any]] = None,
        optimizer_state_dict: Optional[dict[str, Any]] = None,
        max_gradient_value: float = Default.max_gradient_value,
        max_gradient_norm: float = Default.max_gradient_norm,
        correlation_loss_weight: float = Default.correlation_loss_weight,
        num_cpu_processes: int = Default.num_processes
):
    variants_trained = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpu_processes) as process_executor, \
            concurrent.futures.ThreadPoolExecutor(max_workers=2) as thread_executor:
        # create the tensor data loader
        if input_model_path is None:
            training_logger = training_utils.ProgressLogger()
            vcf_tensor_data_loader = vcf_tensor_data_loaders.VcfTrainingTensorDataLoader(
                parquet_path=properties_path, properties_scaling_json=properties_scaling_json, truth_json=truth_json,
                scale_bool_values=scale_bool_values, variants_per_batch=variants_per_batch, shuffle=shuffle,
                random_generator=random_state, excluded_properties=excluded_properties, torch_device=torch_device,
                progress_logger=training_logger, process_executor=process_executor, thread_executor=thread_executor,
                max_look_ahead_batches=max(2, num_cpu_processes)
            )
            # create the neural net, move to requested device, and set to train
            gq_recalibrator = gq_recalibrator_net.GqRecalibratorNet(
                num_input_properties=vcf_tensor_data_loader.num_properties, num_hidden_layers=num_hidden_layers,
                layer_expansion_factor=layer_expansion_factor, bias=bias, hidden_nonlinearity=hidden_nonlinearity,
                output_nonlinearity=output_nonlinearity, model_state_dict=model_state_dict
            ).to(vcf_tensor_data_loader.torch_device).train()
            # create optimizer
            optimizer = optim_type(gq_recalibrator.parameters(), **optim_kwargs)
            if optimizer_state_dict is not None:
                optimizer.load_state_dict(optimizer_state_dict)
            training_losses = training_utils.BatchLosses()
            validation_losses = training_utils.BatchLosses()
        else:
            training_logger, vcf_tensor_data_loader, gq_recalibrator, optimizer, training_losses, validation_losses = \
                load(
                    input_model_path, parquet_path=properties_path, properties_scaling_json=properties_scaling_json,
                    truth_json=truth_json, scale_bool_values=scale_bool_values, variants_per_batch=variants_per_batch,
                    shuffle=shuffle, excluded_properties=excluded_properties, torch_device=torch_device,
                    process_executor=process_executor, thread_executor=thread_executor,
                    max_look_ahead_batches=max(2, num_cpu_processes)
                )

        gq_baseline, gq_scale = vcf_tensor_data_loader.gq_baseline, vcf_tensor_data_loader.gq_scale
        for (is_training_batch, batch_tensor, batch_weights, gq, is_good_gt, is_bad_gt) in vcf_tensor_data_loader:
            if is_training_batch:
                optimizer.zero_grad()
                batch_outputs = gq_recalibrator(batch_tensor)
                loss = loss_function(
                    predicted_probabilities=batch_outputs, variant_weights=batch_weights, gq_values=gq,
                    genotype_is_good=is_good_gt, genotype_is_bad=is_bad_gt,
                    correlation_loss_scale=correlation_loss_weight, gq_baseline=gq_baseline,
                    gq_scale=gq_scale
                )
                loss.backward()
                if max_gradient_value > 0:
                    torch.nn.utils.clip_grad_value_(gq_recalibrator.parameters(), max_gradient_value)
                if max_gradient_norm > 0:
                    torch.nn.utils.clip_grad_norm_(gq_recalibrator.parameters(), max_gradient_norm)
                optimizer.step()

                training_losses.append(loss.item())
                variants_trained += variants_per_batch
                if training_losses.num_mini_epoch_batches >= batches_per_mini_epoch:
                    training_logger.summarize_training_progress(training_losses, validation_losses)
                    training_losses.next_mini_epoch()
                    validation_losses.next_mini_epoch()
                    save(save_path=output_model_path, training_logger=training_logger,
                         vcf_tensor_data_loader=vcf_tensor_data_loader, neural_net=gq_recalibrator,
                         optimizer=optimizer, training_losses=training_losses, validation_losses=validation_losses)

                if max_train_variants is not None and variants_trained >= max_train_variants:
                    break
            else:
                with torch.no_grad():
                    batch_outputs = gq_recalibrator(batch_tensor)

                    validation_losses.append(
                        loss_function(
                            predicted_probabilities=batch_outputs, variant_weights=batch_weights, gq_values=gq,
                            genotype_is_good=is_good_gt, genotype_is_bad=is_bad_gt,
                            correlation_loss_scale=correlation_loss_weight, gq_baseline=gq_baseline,
                            gq_scale=gq_scale
                        ).item()
                    )

    training_logger.summarize_training_progress(training_losses, validation_losses)
    save(save_path=output_model_path, training_logger=training_logger, vcf_tensor_data_loader=vcf_tensor_data_loader,
         neural_net=gq_recalibrator, optimizer=optimizer, training_losses=training_losses,
         validation_losses=validation_losses)

    if torch_device == Keys.cuda:
        torch.cuda.empty_cache()

    training_logger.log(
        f"Got {len(training_losses)} training batches and {len(validation_losses)} validation"
        f" batches in {training_logger.elapsed_time}"
    )


def correlation_coefficient(
        target_values: torch.Tensor,
        predicted_values: torch.Tensor,
        variant_weights: torch.Tensor,
        gq_baseline: torch.Tensor,
        gq_scale: torch.Tensor
) -> torch.Tensor:
    """ Differentiable correlation coefficient between predicted and target values """
    # np.corrcoef in torch from @mdo
    # https://forum.numer.ai/t/custom-loss-functions-for-xgboost-using-pytorch/960
    # subtract mean and divide by norm
    predicted_values = 2.0 * (predicted_values - 0.5)
    target_values = (target_values - gq_baseline) * gq_scale
    #predicted_values = predicted_values / (1e-6 + predicted_values.norm(dim=1, keepdim=True))
    #target_values = target_values / (1e-6 + target_values.norm(dim=1, keepdim=True))
    # compute correlation across samples, then compute weighted sum across variants
    variant_correlation = (predicted_values * target_values).mean(dim=1)
    return (variant_correlation * variant_weights).mean()


def truth_agreement_loss(
        genotype_is_good: torch.Tensor,
        genotype_is_bad: torch.Tensor,
        predicted_probabilities: torch.Tensor,
        variant_weights: torch.Tensor
) -> torch.Tensor:
    """
    Differentiable strength of agreement between probabilities and truth
    Args:
        genotype_is_good:
        genotype_is_bad:
        predicted_probabilities:
        variant_weights:

    Returns:
    """
    # loss = genotype_is_good * (1.0 - predicted_probabilities) + genotype_is_bad * predicted_probabilities
    loss = predicted_probabilities * (genotype_is_bad - genotype_is_good) + genotype_is_good
    variant_loss = loss.sum(dim=1) / (genotype_is_good.sum(dim=1) + genotype_is_bad.sum(dim=1))
    return (variant_loss * variant_weights).mean()


def loss_function(
        predicted_probabilities: torch.Tensor,
        variant_weights: torch.Tensor,
        gq_values: torch.Tensor,
        genotype_is_good: torch.Tensor,
        genotype_is_bad: torch.Tensor,
        correlation_loss_scale: float,
        gq_baseline: torch.Tensor,
        gq_scale: torch.tensor
) -> torch.Tensor:
    tal = truth_agreement_loss(
        genotype_is_good=genotype_is_good, genotype_is_bad=genotype_is_bad,
        predicted_probabilities=predicted_probabilities, variant_weights=variant_weights
    )
    cc = correlation_coefficient(
        target_values=gq_values, predicted_values=predicted_probabilities,
        variant_weights=variant_weights, gq_baseline=gq_baseline, gq_scale=gq_scale
    )
    print(f"tal: {tal}, cc: {cc}")
    return tal + correlation_loss_scale * (1.0 - cc)


def save(
        save_path: Path,
        training_logger: training_utils.ProgressLogger,
        vcf_tensor_data_loader: vcf_tensor_data_loaders.VcfTrainingTensorDataLoader,
        neural_net: gq_recalibrator_net.GqRecalibratorNet,
        optimizer: torch.optim.Optimizer,
        training_losses: training_utils.BatchLosses,
        validation_losses: training_utils.BatchLosses
):
    """ Save neural net and state info for resuming training (or for filtering) """
    save_path.parent.mkdir(parents=True, exist_ok=True)
    pickle_dict = {
        Keys.logger_kwargs: training_logger.save_dict,
        Keys.data_loader_state: vcf_tensor_data_loader.state_dict,
        Keys.net_kwargs: neural_net.to(Keys.cpu).save_dict,
        Keys.optimizer_type: type(optimizer),
        Keys.optimizer_kwargs: state_dict_to(
            state_dict=optimizer.state_dict(),
            device=vcf_tensor_data_loaders.get_torch_device(Keys.cpu,  progress_logger=None)
        ),
        Keys.training_losses: training_losses.save_dict,
        Keys.validation_losses: validation_losses.save_dict
    }
    with open(save_path, "wb") as f_out:
        torch.save(pickle_dict, f_out)


def state_dict_to(state_dict: dict[str, Any], device: torch.device) -> dict[str, Any]:
    """
    move state_dict to desired device, necessary for optimizer
    idea taken from: https://discuss.pytorch.org/t/moving-optimizer-from-cpu-to-gpu/96068/3
    """
    moved_state_dict = {}
    for key, param in state_dict.items():
        if isinstance(param, torch.Tensor):
            moved_param = param.detach().to(device)
            if param._grad is not None:
                moved_param._grad.data = param.detach()._grad.data.to(device)
            moved_state_dict[key] = moved_param
        elif isinstance(param, dict):
            moved_state_dict[key] = state_dict_to(param, device)
        else:
            moved_state_dict[key] = param
    return moved_state_dict


def load(
        save_path: Path,
        parquet_path: Path,
        properties_scaling_json: Path,
        truth_json: Path,
        scale_bool_values: bool,
        variants_per_batch: int,
        shuffle: bool,
        excluded_properties: set[str],
        torch_device: str,
        process_executor: concurrent.futures.ProcessPoolExecutor,
        thread_executor: concurrent.futures.ThreadPoolExecutor,
        max_look_ahead_batches: int
) -> (training_utils.ProgressLogger, vcf_tensor_data_loaders.VcfTrainingTensorDataLoader,
      gq_recalibrator_net.GqRecalibratorNet, torch.optim.Optimizer,
      training_utils.BatchLosses, training_utils.BatchLosses):
    with open(save_path, "rb") as f_in:
        pickle_dict = torch.load(f_in)

    training_logger = training_utils.ProgressLogger(**pickle_dict[Keys.logger_kwargs])
    training_logger.replay()
    vcf_tensor_data_loader = vcf_tensor_data_loaders.VcfTrainingTensorDataLoader(
        parquet_path=parquet_path, properties_scaling_json=properties_scaling_json,
        truth_json=truth_json, scale_bool_values=scale_bool_values,
        variants_per_batch=variants_per_batch, shuffle=shuffle, random_generator=0,
        excluded_properties=excluded_properties, torch_device=torch_device,
        progress_logger=training_logger, process_executor=process_executor,
        thread_executor=thread_executor, max_look_ahead_batches=max_look_ahead_batches
    )
    vcf_tensor_data_loader.load_state_dict(pickle_dict[Keys.data_loader_state])
    neural_net = gq_recalibrator_net.GqRecalibratorNet(**pickle_dict[Keys.net_kwargs])\
        .to(vcf_tensor_data_loader.torch_device).train()
    optimizer = pickle_dict[Keys.optimizer_type](neural_net.parameters())
    optimizer.load_state_dict(
        state_dict_to(
            state_dict=pickle_dict[Keys.optimizer_kwargs],
            device=vcf_tensor_data_loader.torch_device
        )
    )
    training_losses = training_utils.BatchLosses(**pickle_dict[Keys.training_losses])
    validation_losses = training_utils.BatchLosses(**pickle_dict[Keys.validation_losses])
    return (
        training_logger, vcf_tensor_data_loader, neural_net, optimizer, training_losses,
        validation_losses
    )


def filter_sv_gq_recalibrator(
        properties_path: Path,
        properties_scaling_json: Path,
        input_model_path: Path,
        variants_per_batch: int = Default.variants_per_batch,
        shuffle: bool = Default.shuffle,
        random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
        excluded_properties: set[str] = Default.excluded_properties,
        torch_device: str = Default.torch_device,
        scale_bool_values: bool = Default.scale_bool_properties
):
    t0 = time.time()
    num_batches = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as process_executor, \
            concurrent.futures.ThreadPoolExecutor(max_workers=2) as thread_executor:
        vcf_tensor_data_loader = vcf_tensor_data_loaders.VcfFilterTensorDataLoader(
            parquet_path=properties_path, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values, variants_per_batch=variants_per_batch, shuffle=shuffle,
            random_generator=random_state, excluded_properties=excluded_properties, torch_device=torch_device,
            process_executor=process_executor, thread_executor=thread_executor
        )
        for batch_tensor, batch_ids in tqdm(
                vcf_tensor_data_loader, desc="batch", mininterval=0.5, maxinterval=float('inf'), smoothing=0
        ):
            if not isinstance(batch_tensor, torch.Tensor):
                raise ValueError(f"batch {num_batches} is a {type(batch_tensor)}")
            num_batches += 1
    t1 = time.time()
    print(f"Got {num_batches} batches in {t1 - t0} s")
    print(f"batch_ids:\n{batch_ids}")
    print(f"batch_tensor_shape:\n{batch_tensor.shape}")
    print(f"batch_tensor[0, 0, :] = {batch_tensor[0, 0, :]}")


if __name__ == "__main__":
    main()
