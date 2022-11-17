#!/usr/bin/env python

import sys
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
from gq_recalibrator import tarred_properties_to_parquet, vcf_tensor_data_loaders, gq_recalibrator_net
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
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
    id = tarred_properties_to_parquet.Keys.id
    variant_weights = tarred_properties_to_parquet.Keys.variant_weights
    positive_value = tarred_properties_to_parquet.Keys.positive_value
    negative_value = tarred_properties_to_parquet.Keys.negative_value
    baseline = tarred_properties_to_parquet.Keys.baseline
    scale = tarred_properties_to_parquet.Keys.scale


class Default:
    temp_dir = tempfile.gettempdir()
    variants_per_batch = 100
    batches_per_mini_epoch = 100
    shuffle = vcf_tensor_data_loaders.Default.shuffle
    random_state = vcf_tensor_data_loaders.Default.random_state
    excluded_properties = vcf_tensor_data_loaders.Default.excluded_properties
    torch_device = Keys.cpu
    scale_bool_properties = vcf_tensor_data_loaders.Default.scale_bool_properties
    validation_proportion = vcf_tensor_data_loaders.Default.validation_proportion
    num_hidden_layers = gq_recalibrator_net.Default.num_hidden_layers
    layer_expansion_factor = gq_recalibrator_net.Default.layer_expansion_factor
    bias = gq_recalibrator_net.Default.bias
    optim = gq_recalibrator_net.Default.optim
    optim_type = gq_recalibrator_net.Default.optim_type
    optim_kwargs = gq_recalibrator_net.Default.optim_kwargs
    hidden_nonlinearity = gq_recalibrator_net.Default.hidden_nonlinearity
    output_nonlinearity = gq_recalibrator_net.Default.output_nonlinearity
    max_gradient_value = 0.0 # 0.01
    max_gradient_norm = 0.0 # 1.0e-3


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train filter to recalibrate genotype quality from extracted VCF properties",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--properties", "-p", type=str, required=True,
                        help="full path to tarred parquet file with variant properties")
    parser.add_argument("--properties-scaling-json", "-j", type=str, required=True,
                        help="full path to JSON with baseline and scales needed for training / filtering properties")
    parser.add_argument("--output-model", "-o", type=str, required=True,
                        help="path to output file with trained model")
    parser.add_argument("--temp-dir", "-t", type=str, default=Default.temp_dir,
                        help="full path to preferred temporary directory")
    parser.add_argument("--variants-per-batch", type=int, default=Default.variants_per_batch,
                        help="number of variants used in each training batch")
    parser.add_argument("--shuffle", type=common.argparse_bool, default=Default.shuffle,
                        help="if True, shuffle order of variants while training")
    parser.add_argument("--random-seed", type=int, default=Default.random_state,
                        help="initial random seed")
    parser.add_argument("--excluded-properties", type=str, default=','.join(Default.excluded_properties),
                        help="comma-separated list of properties to not use in tensor")
    parser.add_argument("--torch-device", type=str, default=Default.torch_device, choices=[Keys.cpu, Keys.cuda],
                        help="Device on which to perform training")
    parser.add_argument("--scale-bool-properties", type=common.argparse_bool, default=Default.scale_bool_properties,
                        help="If true, z-score bool variables too. If false, leave as 0.0 or 1.0. Scaling can result in"
                             "large values if the bool variable is mostly true or mostly false.")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    train_sv_gq_recalibrator(
        properties_path=args.properties,
        properties_scaling_json=args.properties_scaling_json,
        output_model_path=args.output_model,
        temp_dir=args.temp_dir,
        variants_per_batch=args.variants_per_batch,
        shuffle=args.shuffle,
        random_state=args.random_seed,
        excluded_properties=set(args.excluded_properties.split(',')),
        torch_device=args.torch_device
    )
    # filter_sv_gq_recalibrator(
    #     properties_path=args.properties,
    #     properties_scaling_json=args.properties_scaling_json,
    #     output_model_path=args.output_model,
    #     temp_dir=args.temp_dir,
    #     variants_per_batch=args.variants_per_batch,
    #     shuffle=args.shuffle,
    #     random_state=args.random_seed,
    #     excluded_properties=set(args.excluded_properties.split(',')),
    #     torch_device=args.torch_device
    # )


def train_sv_gq_recalibrator(
        properties_path: str,
        properties_scaling_json: str,
        truth_json: str,
        output_model_path: str,
        temp_dir: str = Default.temp_dir,
        variants_per_batch: int = Default.variants_per_batch,
        batches_per_mini_epoch: int = Default.batches_per_mini_epoch,
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
        max_gradient_norm: float = Default.max_gradient_norm
):
    # some kind of checkpoint load here
    t0 = time.time()
    training_losses = BatchLosses()
    validation_losses = BatchLosses()

    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        # create the tensor data loader
        vcf_tensor_data_loader = vcf_tensor_data_loaders.VcfTrainingTensorDataLoader(
            parquet_path=properties_path, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values, variants_per_batch=variants_per_batch, shuffle=shuffle,
            random_state=random_state, excluded_properties=excluded_properties, torch_device=torch_device,
            executor=executor
        )
        # create the neural net, move to requested device, and set to train
        gq_recalibrator = gq_recalibrator_net.GqRecalibratorNet(
            num_input_properties=vcf_tensor_data_loader.num_properties, num_hidden_layers=num_hidden_layers,
            layer_expansion_factor=layer_expansion_factor, bias=bias, optim_type=optim_type, optim_kwargs=optim_kwargs,
            hidden_nonlinearity=hidden_nonlinearity, output_nonlinearity=output_nonlinearity,
            model_state_dict=model_state_dict
        ).to(vcf_tensor_data_loader.torch_device).train()
        # create optimizer
        optimizer = optim_type(gq_recalibrator.parameters(), **optim_kwargs)
        if optimizer_state_dict is not None:
            optimizer.load_state_dict(optimizer_state_dict)

        for (
                batch_tensor, batch_weights, batch_ids,
                validation_tensor, validation_weights, validation_ids
        ) in vcf_tensor_data_loader:
            optimizer.zero_grad()
            batch_outputs = gq_recalibrator(batch_tensor)
            loss = objective_function(batch_outputs, batch_targets, batch_weights)
            loss.backward()
            training_losses.append(loss.item())
            if max_gradient_value > 0:
                torch.nn.utils.clip_grad_value_(gq_recalibrator.parameters(), max_gradient_value)
            if max_gradient_norm > 0:
                torch.nn.utils.clip_grad_norm_(gq_recalibrator.parameters(), max_gradient_norm)
            optimizer.step()

            if validation_tensor is not None:
                with torch.no_grad():
                    validation_outputs = gq_recalibrator(validation_tensor)
                    validation_losses.append(
                        objective_function(validation_outputs, batch_targets, validation_weights).item()
                    )

            if training_losses.num_mini_epoch_batches >= batches_per_mini_epoch:
                _summarize_training_progress(training_losses, validation_losses, t0)
                # somehow save a checkpoint here
                training_losses.next_mini_epoch()
                validation_losses.next_mini_epoch()

    _summarize_training_progress(training_losses, validation_losses, t0)

    print(f"Got {len(training_losses)} training batches and {len(validation_losses)} validation batches in "
          f"{common.elapsed_time(time.time() - t0)}")
    print(f"batch_ids:\n{batch_ids}")
    print(f"batch_weights:\n{batch_weights}")
    print(f"batch_tensor_shape:\n{batch_tensor.shape}")
    print(f"batch_tensor[0, 0, :] = {batch_tensor[0, 0, :]}")


class BatchLosses:
    __slots__ = ("loss_history", "num_mini_epoch_batches", "num_mini_epochs")

    def __init__(self, loss_history: list[float] = (), num_mini_epoch_batches: int = 0, num_mini_epochs: int = 0):
        self.loss_history = list(loss_history)
        self.num_mini_epoch_batches = num_mini_epoch_batches
        self.num_mini_epochs = num_mini_epochs

    def append(self, loss: float):
        self.loss_history.append(loss)
        self.num_mini_epoch_batches += 1

    def __len__(self):
        return len(self.loss_history)

    def next_mini_epoch(self):
        self.num_mini_epoch_batches = 0
        self.num_mini_epochs += 1

    @property
    def mini_epoch_loss(self) -> float:
        return sum(self.loss_history[-self.num_mini_epoch_batches:]) / self.num_mini_epoch_batches


def _summarize_training_progress(training_losses: BatchLosses, validation_losses: BatchLosses, t0: float):
    elapsed_time = common.elapsed_time(time.time() - t0)
    if training_losses.num_mini_epochs == 0:
        # this is the first mini-epoch
        print(f"epoch: {'elapsed-time':11s} {'train-loss':10s} {'valid-loss':10s}")
    print(f"{training_losses.num_mini_epochs:5d}: {elapsed_time:11s} {training_losses.mini_epoch_loss:<10.3f}"
          f" {validation_losses:<10.3f}")


def correlation_coefficient(
        target_values: torch.Tensor,
        predicted_values: torch.Tensor,
        variant_weights: torch.Tensor
) -> torch.Tensor:
    """ Differentiable correlation coefficient between predicted and target values """
    # np.corrcoef in torch from @mdo
    # https://forum.numer.ai/t/custom-loss-functions-for-xgboost-using-pytorch/960
    # subtract mean and divide by norm
    predicted_values = predicted_values - predicted_values.mean(dim=1)
    target_values = target_values - target_values.mean(dim=1)
    predicted_values = predicted_values / predicted_values.norm(dim=1)
    target_values = target_values / target_values.norm(dim=1)
    # compute correlation across samples, then compute weighted sum across variants
    return ((predicted_values * target_values).mean(dim=1) * variant_weights).sum()


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
    loss = predicted_probabilities * (genotype_is_bad - genotype_is_good) - genotype_is_bad
    return (loss.sum(dim=1) / (genotype_is_good.sum(dim=1) + genotype_is_bad.sum(dim=1)) * variant_weights).sum()


def loss_function(
        gq_values: torch.Tensor,
        predicted_probabilities: torch.Tensor,
        genotype_is_good: torch.Tensor,
        genotype_is_bad: torch.Tensor,
        variant_weights: torch.Tensor,
        correlation_loss_scale: float
) -> torch.Tensor:
    return truth_agreement_loss(
        genotype_is_good=genotype_is_good, genotype_is_bad=genotype_is_bad,
        predicted_probabilities=predicted_probabilities, variant_weights=variant_weights
    ) + correlation_loss_scale * (
        1.0 - correlation_coefficient(target_values=gq_values, predicted_values=predicted_probabilities,
                                      variant_weights=variant_weights)
    )


def filter_sv_gq_recalibrator(
        properties_path: str,
        properties_scaling_json: str,
        output_model_path: str,
        temp_dir: str = Default.temp_dir,
        variants_per_batch: int = Default.variants_per_batch,
        shuffle: bool = Default.shuffle,
        random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
        excluded_properties: set[str] = Default.excluded_properties,
        torch_device: str = Default.torch_device,
        scale_bool_values: bool = Default.scale_bool_properties
):
    t0 = time.time()
    num_batches = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        vcf_tensor_data_loader = vcf_tensor_data_loaders.VcfFilterTensorDataLoader(
            parquet_path=properties_path, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values, variants_per_batch=variants_per_batch, shuffle=shuffle,
            random_state=random_state, excluded_properties=excluded_properties, torch_device=torch_device,
            executor=executor
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
