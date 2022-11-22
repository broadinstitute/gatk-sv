import time
import numpy
from sv_utils import common
from typing import Optional, Any


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
        # noinspection PyTypeChecker
        return numpy.mean(self.loss_history[-self.num_mini_epoch_batches:])

    @property
    def save_dict(self) -> dict[str, Any]:
        return {_slot: getattr(self, _slot) for _slot in self.__slots__}


class ProgressLogger:
    __slots__ = ("logs", "t_start", "t_last", "_has_header")

    def __init__(
            self,
            logs: list[str] = (),
            t_start: Optional[float] = None,
            t_last: Optional[float] = None,
            _has_header: bool = False
    ):
        self.logs = list(logs)
        if t_start is None:
            # start now
            self.t_start = time.time()
            self.t_last = self.t_start
        else:
            # construct pseudo timeline where we started delta_t ago
            delta_t = t_last - t_start
            self.t_last = time.time()
            self.t_start = self.t_last - delta_t
        self._has_header = _has_header

    @property
    def delta_t(self) -> float:
        return time.time() - self.t_start

    @property
    def elapsed_time(self) -> str:
        return common.elapsed_time(self.delta_t, seconds_precision=0)

    def log(self, log_str, add_timestamp: bool = True):
        self.t_last = time.time()
        if add_timestamp:
            log_str = time.strftime("%Y-%m-%d %H:%M:%S %Z", time.localtime(self.t_last)) + " " + log_str
        self.logs.append(log_str)
        print(log_str)

    def replay(self):
        for log in self.logs:
            print(log)
        self.log("RESTART")

    def __len__(self) -> int:
        return len(self.logs)

    @property
    def save_dict(self) -> dict[str, Any]:
        return {_slot: getattr(self, _slot) for _slot in self.__slots__}

    def summarize_training_progress(self, training_losses: BatchLosses, validation_losses: BatchLosses):
        if not self._has_header:
            # this is the first log, print a header
            self.log(f"{'timestamp':23s} epoch {'elapsed-time':12s} {'train-loss':10s} {'valid-loss':10s}",
                     add_timestamp=False)
            self._has_header = True
        self.log(f"{training_losses.num_mini_epochs:5d} {self.elapsed_time:12s}"
                 f" {training_losses.mini_epoch_loss:<10.3f} {validation_losses.mini_epoch_loss:<10.3f}")
