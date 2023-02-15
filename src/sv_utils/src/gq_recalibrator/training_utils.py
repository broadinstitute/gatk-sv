import sys
import os
import traceback
import time
import enum
import signal
import psutil
import numpy
import torch
import concurrent.futures
import functools
from sv_utils import common
from typing import Optional, Any


def kill_child_processes(parent_pid: int, sig: signal.Signals = signal.SIGTERM):
    """
    If a concurrent job is stuck and needs to be canceled, you may need to kill all child processes
    from:
    https://stackoverflow.com/questions/42782953/python-concurrent-futures-how-to-make-it-cancelable
    """
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        return
    children = parent.children(recursive=True)
    for process in children:
        process.send_signal(sig)


def reraise_with_stack(func):
    """decorator to allow preservation of stack info for Exceptions raised in futures"""

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            raise sys.exc_info()[0](traceback.format_exc())

    return wrapped


class _ExceptionStrWithNewlines(str):
    def __repr__(self): return f"{self}"


def get_result_or_raise(
        future: concurrent.futures.Future,
        executor: Optional[concurrent.futures.Executor] = None,
        context: Optional[str] = None
) -> Any:
    """Return result from future or raise the exception that terminated it.
    On exception:
     -Shut down execution of all remaining jobs for a clean exit.
     -When re-raising the exception, add the provided context description.

    Args:
        future: the future object encapsulating the concurrent job
        executor: the executor executing the concurrent job. If not None, shut down all remaining
                  jobs on this executor for a clean exit.
        context: If not None, add context to exception stack
    """
    # because timeout is None, this check for exception will wait until execution is complete:
    if future.exception(timeout=None) is not None:
        # need to shut down other processes before we can propagate exception
        if executor is not None:
            executor.shutdown(wait=False, cancel_futures=True)
            kill_child_processes(parent_pid=os.getpid())
        if context is None:
            context = "in child process"
        raise RuntimeError(
            _ExceptionStrWithNewlines(f"{context}: {future.exception()}")
        )
    else:
        return future.result()


class BatchLosses:
    __slots__ = ("loss_history", "num_mini_epoch_batches", "num_mini_epochs")

    def __init__(
            self,
            loss_history: list[float] = (),
            num_mini_epoch_batches: int = 0,
            num_mini_epochs: int = 0
    ):
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
            log_str = time.strftime(
                "%Y-%m-%d %H:%M:%S %Z", time.localtime(self.t_last)
            ) + " " + log_str
        self.logs.append(log_str)
        print(log_str)

    def __call__(self, *args, **kwargs):
        self.log(*args, **kwargs)

    def replay(self):
        for log in self.logs:
            print(log)
        self.log("RESTART")

    def __len__(self) -> int:
        return len(self.logs)

    @property
    def save_dict(self) -> dict[str, Any]:
        return {_slot: getattr(self, _slot) for _slot in self.__slots__}

    def summarize_training_progress(
            self,
            training_losses: BatchLosses,
            validation_losses: BatchLosses
    ):
        if not self._has_header:
            # this is the first log, print a header
            self.log(
                f"{'timestamp':23s} epoch {'elapsed-time':12s} {'train-loss':10s} "
                f"{'valid-loss':10s}",
                add_timestamp=False
            )
            self._has_header = True
        self.log(
            f"{training_losses.num_mini_epochs:5d} {self.elapsed_time:12s} "
            f"{training_losses.mini_epoch_loss:<10.3f} "
            f"{validation_losses.mini_epoch_loss:<10.3f}"
        )


class TorchDeviceKind(enum.Enum):
    cuda = "cuda"
    cpu = "cpu"
    mps = "mps"  # Macbook Performance Shaders
    
    @property
    def value_str(self) -> str:
        return str(self.value)

    # noinspection PyTypeChecker,PyMethodParameters
    @common.classproperty
    def choices(cls) -> list[str]:
        return [kind.name for kind in cls]

    def get_device(self, progress_logger: Optional[ProgressLogger] = None) -> torch.device:
        """Return appropriate torch.device based on requested device

        Args:
            progress_logger: Optional[training_utils.ProgressLogger]
                If a progress logger, use that to log messages. Otherwise print.
        Returns:
            torch_device: torch.device
                Device to use for tensor calculations
        """
        def _log(_log_str: str):
            if progress_logger is not None:
                progress_logger(_log_str)
        if self == TorchDeviceKind.cuda:
            if torch.cuda.is_available():
                _log(f"Using {self.name} for torch")
                torch.cuda.empty_cache()
                return torch.device("cuda0")
            else:
                _log(f"{self.name} not available")
                return TorchDeviceKind.cpu.get_device(progress_logger=progress_logger)
        elif self == TorchDeviceKind.mps:
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                _log(f"Using {self.name} for torch")
                return torch.device("mps")
            else:
                _log(f"{TorchDeviceKind.mps} not available")
                return TorchDeviceKind.cpu.get_device(progress_logger=progress_logger)
        else:
            _log(f"Using {self.name} for torch")
            return torch.device("cpu")
