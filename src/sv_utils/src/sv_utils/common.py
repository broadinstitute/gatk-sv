import sys
import os
import importlib
import warnings
import subprocess
import numpy
import pandas
import collections
import psutil
from typing import Text, Any, Union, List, Optional, Iterable, Iterator, Tuple, Generic, Callable, TypeVar
from types import ModuleType


Vector = Union[List, pandas.Series, numpy.array]
Numeric = Union[int, float, numpy.integer, numpy.floating]
TypeT = TypeVar("TypeT")


def deref_iterable_with_sorted_indices(base_list: Iterable[TypeT], sorted_indices: Iterable[int]) -> Iterator[TypeT]:
    """
    Make iterator to values of iterable data structure dereferenced by a sorted iterable of wanted indices.
    Note that this function does not check that sorted_indices are actually sorted.
    Args:
        base_list: Iterable[TypeT]
            Iterable with values to dereference
        sorted_indices: Iterable[int]
            Sorted iterable of integer indices of desired elements of base_list
    Yields:
        dereferenced_vals: TypeT
            Desired elements of base_list
    """
    base_iter = iter(enumerate(base_list))
    for next_sorted_ind in sorted_indices:
        base_ind, base_value = next(base_iter)
        while base_ind != next_sorted_ind:
            base_ind, base_value = next(base_iter)
        yield base_value


def true(size: Union[int, Tuple[int, ...]]) -> numpy.ndarray:
    """
    Syntactic sugar to get a bool array of all true.
    Args:
        size: int or tuple of ints
            size of the array to generate
    Returns:
        true_array: numpy.ndarray
            bool array of all true
    """
    return numpy.ones(size, dtype=bool)


def false(size: Union[int, Tuple[int, ...]]) -> numpy.ndarray:
    """
    Syntactic sugar to get a bool array of all false.
    Args:
        size: int or tuple of ints
            size of the array to generate
    Returns:
        true_array: numpy.ndarray
            bool array of all false
    """
    return numpy.zeros(size, dtype=bool)


class PeekableIter(Generic[TypeT]):
    """
    Iterator that can be "copied" (with .peek) to allow peeking ahead while allowing subsequent iteration over all the
    values in the main iterator (even those that were "peeked").
    """
    __slots__ = ('iterator', 'queue', 'index')

    def __init__(self, iterable: Iterator[TypeT]):
        self.iterator = iterable if isinstance(iterable, Iterator) \
            else iter(iterable)
        self.queue = collections.deque()

    def __next__(self) -> TypeT:
        if self.queue:
            return self.queue.popleft()
        else:
            return next(self.iterator)

    def __iter__(self):
        return self

    def peek(self) -> "PeekIter[TypeT]":
        """ return an iterator over the values in this iterator that doesn't consume the values from this iterator """
        return PeekIter(self)

    def peek_next(self) -> TypeT:
        """ return the next value in this iterator, without consuming it """
        if self.queue:
            return self.queue[0]
        else:
            value = next(self.iterator)
            self.queue.appendleft(value)
            return value

    def has_next(self) -> bool:
        """ overload of Iterator has_next, dealing with peek overhead """
        if self.queue:
            return True
        try:
            self.peek_next()
            return True
        except StopIteration:
            return False

    def put_back(self, value: TypeT):
        """ In case you want to put something back that was consumed instead of peeked. Useful when you usually want
            to consume, but find out you were wrong afterwards """
        self.queue.appendleft(value)


class PeekIter(Generic[TypeT]):
    """
    Helper class for PeekableIter, it does the peaking by always extending the internal queue
    """
    __slots__ = ('iterator', 'queue', 'count', 'queue_iter')

    def __init__(self, parent_iter: PeekableIter[TypeT]):
        self.iterator = parent_iter.iterator
        self.queue = parent_iter.queue
        self.count = len(self.queue)
        self.queue_iter = iter(self.queue)

    def __next__(self) -> TypeT:
        if self.count:
            self.count -= 1
            return next(self.queue_iter)
        else:
            x = next(self.iterator)
            self.queue.append(x)
            return x

    def __iter__(self):
        return self


def command_results(command: Text, exception_on_stderr: bool = True) -> Text:
    """
    Execute shell command. Raise exception if unsuccessful, otherwise return string with output
    """
    sub_p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=True)
    with sub_p.stdout as pipeIn, sub_p.stderr as pipeErr:
        results = pipeIn.read().decode('utf-8')
        err = pipeErr.read().decode('utf-8')
    if err:
        if exception_on_stderr:
            raise RuntimeError('Error executing %s:\n%s' % (command, err[:-1]))
        else:
            warnings.warn(err[:-1])
    return results


def add_exception_context(exception: Exception, context: str):
    """
    Add additional context to a caught exception
    Args:
        exception: Exception
            Exception that was caught
        context: str
            Extra info to add to exception.
    Returns:
        exception: Exception
            Original exeption with annotated context.
    """
    if len(exception.args) == 1 and type(exception.args[0]) is str:
        exception.args = (context, exception.args[0])
    else:
        exception.args = (context,) + exception.args


def static_vars(**kwargs):
    """
    Decorate a function to add static variables. e.g.:
    @static_vars(my_var=6, my_other_var=None)
    def my_func(foo, bar):
        foo += my_func.my_var
        if my_func.my_other_var is not None:
            print(foo + bar)

    Args:
        **kwargs: dict
            keys are variable names, vals are initial values
            Technically they are just added as attributes of decorated function,
            but they fill the roll of static variables

    Returns:
        _decorate: Callable
            function decorator that adds "static" variables
    """
    def _decorate(func):
        for key, val in kwargs.items():
            setattr(func, key, val)
        return func
    return _decorate


# noinspection PyPep8Naming
class classproperty(object):
    """
    Decorate a class member function to make it a property, essentially combining
    @property
    @classmethod
    """
    def __init__(self, getter: Callable[[type], Any]):
        self.getter = getter

    def __get__(self, instance, owner: type) -> Any:
        return self.getter(owner)


def dynamic_import(
        obj_name: Text,
        base: Union[Text, ModuleType, None] = None
) -> Any:
    """
    Import an object by prop_name and return it. Can descend object hierarchy by import_module() or getattr().
    Args:
        obj_name: str
            Name of object in hierarchy, with layers_size separated by '.'
            e.g. "scipy.jax_stats.beta"
        base: str, ModuleType, or None (Default=None)
            Base package / object for import.
            If None, import from global namespace.
            If a str, import from package with that prop_name.
    Returns:
        obj: Any
            Imported object
    """
    if isinstance(base, str):
        base = importlib.import_module(base)
    while '.' in obj_name:
        top_level, obj_name = obj_name.split('.', 1)
        if base is None:
            base = importlib.import_module(top_level)
        elif hasattr(base, top_level):
            base = getattr(base, top_level)
        else:
            # package_name = base if isinstance(base, str) else base.__name__
            base = importlib.import_module('.' + top_level, package=base.__name__)

    if base is None:
        base = importlib.import_module(obj_name)
    elif hasattr(base, obj_name):
        base = getattr(base, obj_name)
    else:
        # package_name = base if isinstance(base, str) else base.__name__
        base = importlib.import_module('.' + obj_name, package=base.__name__)
    return base


def get_cpuinfo() -> (int, int):
    """
    Partially cross-platform method to get physical and logical cpu info.
    Returns:
        num_cores: int
            Number of physical cores in the system
        num_hyperthreads: int
            Number of logical cores in the system
    """
    if sys.platform == 'linux':  # pragma: no cover
        _processor_key = "processor"  # This key should be here
        _core_id_key = "core id"
        core_ids = set()
        num_hyperthreads = 0
        with open("/proc/cpuinfo", 'r') as f_in:
            for line in f_in:
                line = line.strip()
                if not line:
                    continue
                key, value = tuple(word.strip() for word in line.split(':', 1))
                if key == _processor_key:
                    num_hyperthreads += 1
                elif key == _core_id_key:
                    core_ids.add(value)
        # a typical physical processor will have one unique "core id" per physical core, potentially listed multiple
        # times if there are hyperthreads. However emulated linux processors may lack a "core id", in which case assume
        # that the number of physical cores is equal to the number of hyperthreads
        num_cores = len(core_ids) if core_ids else num_hyperthreads
    elif sys.platform == "darwin":  # pragma: no cover
        num_hyperthreads, num_cores = tuple(
            int(word)
            for word in command_results("sysctl -n hw.logicalcpu hw.physicalcpu").split()
        )
    else:  # pragma: no cover
        raise OSError(f"platform '{sys.platform}' is not supported")

    return num_cores, num_hyperthreads


def available_memory(return_gb: bool = True) -> float:
    """
    Return memory available for use by this program in GiB

    Args:
        return_gb: bool
            If True, return available memory in GiB; otherwise return in bytes.
    Returns:
        mem: float
            amount of memory available
    """
    mem = psutil.virtual_memory().total - psutil.virtual_memory().wired \
        if psutil.OSX else psutil.virtual_memory().available
    return mem / 2.0**30 if return_gb else float(mem)


try:
    num_physical_cpus, num_logical_cpus = get_cpuinfo()
except OSError:
    import multiprocessing
    num_physical_cpus = multiprocessing.cpu_count()
    num_logical_cpus = num_physical_cpus
hyperthread_ratio = max(1, num_logical_cpus // num_physical_cpus)
max_allowed_hyperthread = hyperthread_ratio * int(os.environ['NSLOTS']) if 'NSLOTS' in os.environ \
    else max(num_logical_cpus, num_physical_cpus)
max_allowed_physical = int(os.environ.get('NSLOTS', num_physical_cpus))
max_allowed_jobs = max_allowed_hyperthread


def num_jobs_to_use(
        num_jobs: Optional[int] = None,
        required_worker_memory: Optional[Numeric] = 0,
        required_master_memory: Optional[Numeric] = 0,
        require_physical_cpus: bool = False
) -> int:
    """
    Return number of jobs to actually use, based on num_jobs and constraints from available system resources (number of
    physical cpus, number of cpu threads, and memory) and the environment variable NSLOTS (if set).
    Args:
        num_jobs: int or None
            Number of jobs requested to use.
            if num_jobs < 0 or None: try to use NSLOTS if set, otherwise cpu count (physical or threads)
            if num_jobs >= 0: try to use max(num_jobs, 1)
        required_worker_memory: int (Default=0)
            Memory (in GB) that each worker will need to execute a task. If the requested number of jobs will not fit
            into memory, cap the number of created jobs so that the task will finish successfully.
        required_master_memory: int (Default=0)
            Memory (in GB) that the master process will need. If required_worker_memory > 0, this will decrease cap on
            jobs due to memory constraints.
        require_physical_cpus: bool (Default=False)
            If True: cap num_jobs at max_allowed_physical
            If False: cap num_jobs at max_allowed_hyperthread
    Returns:
        num_jobs: int
            Number of jobs to actually use. Will be >= 1
    """
    if required_worker_memory is None:
        required_worker_memory = 0
    if required_master_memory is None:
        required_master_memory = 0
    max_allowed = max_allowed_physical if require_physical_cpus else max_allowed_hyperthread
    num_jobs = max_allowed if num_jobs is None or num_jobs < 0 else min(max(num_jobs, 1), max_allowed)
    if required_worker_memory > 0:
        memory_capped_jobs = (available_memory() - required_master_memory) / required_worker_memory
        num_jobs = min(num_jobs, max(1, int(memory_capped_jobs)))
    elif required_master_memory > available_memory():
        warnings.warn(
            f"Operating in serial because require {required_master_memory} GiB but have {available_memory()} GiB"
        )
        num_jobs = 1
    return num_jobs
