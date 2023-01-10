from collections.abc import Mapping
from typing import Any, Optional
import torch
from gq_recalibrator import vcf_tensor_data_loaders


class Keys:
    cpu: vcf_tensor_data_loaders.Keys.cpu


class Default:
    num_hidden_layers = 2
    layer_expansion_factor = 2.0
    bias = True
    leaky_slope = 0.1
    hidden_nonlinearity = torch.nn.LeakyReLU(negative_slope=leaky_slope) if leaky_slope >= 0 else \
        torch.nn.ReLU()
    output_nonlinearity = torch.nn.Sigmoid()


class GqRecalibratorNet(torch.nn.Module):
    __save_values__ = (
        "num_input_properties", "num_hidden_layers", "layer_expansion_factor", "bias",
        "hidden_nonlinearity", "output_nonlinearity", "model_state_dict"
    )

    def __init__(
            self,
            num_input_properties: int,
            num_hidden_layers: int = Default.num_hidden_layers,
            layer_expansion_factor: float = Default.layer_expansion_factor,
            bias: bool = Default.bias,
            hidden_nonlinearity: torch.nn.Module = Default.hidden_nonlinearity,
            output_nonlinearity: torch.nn.Module = Default.output_nonlinearity,
            model_state_dict: Optional[dict[str, Any]] = None

    ):
        super().__init__()
        self.num_input_properties = num_input_properties
        self.num_hidden_layers = num_hidden_layers
        self.layer_expansion_factor = layer_expansion_factor
        self.bias = bias
        self.hidden_nonlinearity = hidden_nonlinearity
        self.output_nonlinearity = output_nonlinearity

        self.layers = self.build()

        if model_state_dict is None:
            self.apply(self.__class__.weights_init)
        else:
            # noinspection PyTypeChecker
            self.load_state_dict(model_state_dict)

    def build(self) -> torch.nn.ModuleList:
        layers_size = self.num_input_properties
        layers = []
        for _ in range(self.num_hidden_layers):
            previous_size, layers_size = layers_size, round(self.layer_expansion_factor * layers_size)
            layers.append(torch.nn.Linear(in_features=previous_size, out_features=layers_size, bias=self.bias))
            layers.append(self.hidden_nonlinearity)
        previous_size, layers_size = layers_size, 1
        layers.append(torch.nn.Linear(in_features=previous_size, out_features=layers_size, bias=self.bias))
        layers.append(self.output_nonlinearity)
        return torch.nn.ModuleList(layers)

    @staticmethod
    def weights_init(m: torch.nn.Module, leaky_slope: float = Default.leaky_slope):
        f"""
        Recursively initialize weights in torch module
        Args:
            m: torch.nn.Module
                Torch module with weights that need to be initialized (possibly in sub-modules)
            leaky_slope: float (default={Default.leaky_slope})
                Slope for Leaky ReLU units. If 0, then basic ReLU is used  
        """
        if isinstance(m, torch.nn.Linear):
            torch.nn.init.kaiming_normal_(m.weight.data, a=leaky_slope, mode='fan_in',
                                          nonlinearity="leaky_relu" if leaky_slope > 0 else "relu")
            if m.bias is not None:
                torch.nn.init.zeros_(m.bias.data)
        elif isinstance(m, torch.nn.BatchNorm1d):
            torch.nn.init.normal_(m.weight.data, mean=0.0, std=1.0)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        for m in self.layers:
            x = m(x)
        return x.squeeze(2)  # remove trivial final dimension

    @property
    def model_state_dict(self) -> dict[str, Any]:
        return self.state_dict()

    @property
    def save_dict(self) -> dict[str, Any]:
        return {key: getattr(self, key) for key in GqRecalibratorNet.__save_values__}

