#!/usr/bin/env python

import sys  # just for sys.stdout.flush
import os
import psutil  # used for getting process number only:
import tempfile
import string
import dill
import random
import numpy
from tqdm.auto import tqdm as tqdm
from tqdm import TqdmWarning
# noinspection PyUnresolvedReferences
import multiprocessing
import multiprocessing.pool
from typing import Callable, Union, Optional, Sized, Tuple, Iterator, TypeVar, Collection, Iterable, Any

from sv_utils import common

Numeric = Union[int, float, numpy.integer, numpy.floating]
TaskType = TypeVar("TaskType")
OutputType = TypeVar("OutputType")
ArgsTypes = Tuple[Any, ...]
BasicFuncType = Callable[[TaskType, ArgsTypes], Union[OutputType, Iterable[OutputType]]]
FuncType = Union[
    BasicFuncType,
    Callable[[*TaskType, ArgsTypes], Union[OutputType, Iterable[OutputType]]]
] if isinstance(TaskType, Iterable) else BasicFuncType


tqdm.monitor_interval = 0
# global variables for distributing tasks
_pool_func_args = dict()
_pool_key_len = 10


class Default:
    update_time = 0.5  # seconds
    num_chunk_divs = 20  # chunks per worker per job
    n_jobs = -1  # by default, use all available processes


def get_process_num() -> int:
    """
    Function to return process number of executing worker
    """
    my_pid = os.getpid()
    sibling_pids = [ch.pid for ch in psutil.Process(pid=os.getppid()).children()]
    return sibling_pids.index(my_pid)


def pmap(
        pool: Optional[multiprocessing.Pool],
        func: FuncType,
        tasks: Union[Iterable[TaskType], Iterator[TaskType]],
        num_tasks: Optional[int] = None,
        task_sizes: Collection[float] = None,
        starmap: bool = False,
        flatmap: bool = False,
        ordered: bool = True,
        permute_evaluation: bool = False,
        description: str = 'parallel tasks',
        update_time: Optional[float] = Default.update_time,
        args: Collection = None, kwargs: dict = None,
        chunksize: Optional[int] = None,
        num_chunk_divs: int = Default.num_chunk_divs
) -> Iterator[OutputType]:
    """
    execute func on each task in tasks, using parallel pool
    progress is displayed using tqdm module.
    Typical invocation:
      # NOTE: using 'spawn' context saves a lot of memory
      with multiprocessing.get_context('spawn').Pool(processes=n_jobs) as pool:
          results = [r for r in pmap(pool, func, tasks)]
      assert results == [func(task) for task in tasks]
    INPUTS
        pool: multiprocessing.Pool
        func: function to execute on tasks
        tasks: list, tuple, or iterator yielding data
            each task is the arguments that will be passed to func
        num_tasks: int (Default=None)
            Number of tasks, used for displaying progress.
            If num_tasks is None, it will be set to len(tasks) if available.
        task_sizes: list or tuple (Default=[])
            Task sizes for displaying progress. If empty, each task will be
            assigned size 1
        starmap: bool (Default=False)
            if True, evaluate func(*task)  (i.e. treat task as an iterable)
            if False, evalute func(task)  (i.e. treat task as a single arg)
        flatmap: bool (Default=False)
            if True, treat output as collection of iterables and yield
                successively from each iterable. e.g.
                for iterable in result:
                    for val in iterable:
                        yield val
            if False, just yield results
        ordered: bool (Default=True)
            if True, yield results in order corresponding to tasks
                (i.e. results[n] = func(tasks[n]))
            if False, yield results in order they are completed by parallel pool
        permute_evaluation: bool (Default=False)
            if True, pass tasks to pool in permuted order. This can help
                with load balancing to make evaluation time estimate more
                accurate.
                NOTE: if ordered is True, results will still
                      correspond to the tasks that are passed in.
                ALSO: if tasks were a generator, they will be stored in a list
                      for permutation. This may affect compute / memory usage.
            if False, pass each task to pool in order of tasks collection
        description: str (Default='parallel tasks')
            Descriptive text for progress bar
        update_time: float (Default=0.5)
            Minimum time in seconds for progress bar updates. If None or inf, progress will not be displayed.
        args: list or tuple (Default=None)
            Additional args to pass to func. Constant across tasks:
            results = [func(task, *args) for task in tasks]
        kwargs: dict (Default=None)
            kwargs to pass to func. Constant across tasks:
            results = [func(task, **kwargs) for task in tasks]
        chunksize: int (Default=None)
            chunksize to pass to imap_unordered (number of tasks to pass to each
            worker in a block). If an integer, use that chunksize, if None, then
            attempt to pick an efficient number (num_chunk_divs chunks per
            worker, or chunksize = 1, whichever is larger)
        num_chunk_divs: int (Default=20)
            number of chunks per worker
    OUTPUT:
        results will be yielded either in task or evaluation order, as
        specified specified by the "ordered" keyword.
    """
    sys.stdout.flush()

    # get / manipulate input arguments into final forms used by map
    num_workers, args, kwargs, tasks, task_sizes, num_tasks, num_tasks_str, total, disable, chunksize \
        = _validate_map_params(
            pool, args, kwargs, tasks, task_sizes, num_tasks, update_time, chunksize, num_chunk_divs, permute_evaluation
        )

    if pool is None:
        print('Executing %s on %s tasks serially'
              % (str(func).split()[1], num_tasks_str))
        result_gen = \
            ((t[0], func(*t[1], *args, **kwargs)) for t in tasks) if starmap\
            else ((t[0], func(t[1], *args, **kwargs)) for t in tasks)
    else:
        if chunksize == 1:
            print('Executing %s on %s tasks with %d parallel workers'
                  % (str(func).split()[1], num_tasks_str, num_workers))
        else:
            print(
                'Executing %s on %s tasks (chunksize=%d) with %d parallel '
                'workers' % (str(func).split()[1], num_tasks_str, chunksize,
                             num_workers)
            )
        # translate func to enable re-sorting tasks
        func = _UnorderedMapTranslator(func, args, kwargs, starmap)
        # noinspection PyArgumentList
        result_gen = pool.imap_unordered(func, tasks, chunksize=chunksize)
    sys.stdout.flush()

    num_completed_tasks = 0
    # create info necessary to re-sort tasks
    next_i = 0
    r_dict = dict()
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=TqdmWarning)
        with tqdm(total=total, disable=disable, mininterval=update_time, maxinterval=float('inf'), smoothing=0,
                  desc=description) as progress:
            # loop (with progress)
            for i, r in result_gen:
                if task_sizes:
                    num_completed_tasks += 1
                    progress.set_postfix(refresh=False, tasks=f"{num_completed_tasks}/{num_tasks}")
                    progress.update(task_sizes[i])
                else:
                    progress.update(1)

                if ordered:
                    if i == next_i:
                        # get next wanted result, yield it
                        if flatmap:
                            yield from r
                        else:
                            yield r
                        next_i += 1
                        while next_i in r_dict:
                            if flatmap:
                                yield from r_dict.pop(next_i)
                            else:
                                yield r_dict.pop(next_i)
                            next_i += 1
                    else:
                        # this result is not wanted yet. store for later
                        r_dict[i] = r
                else:
                    if flatmap:
                        yield from r
                    else:
                        yield r


def _validate_map_params(
        pool: Optional[multiprocessing.Pool],
        args: Collection,
        kwargs: dict,
        tasks: Union[Iterable, Iterator],
        task_sizes: Iterable[float],
        num_tasks: Optional[int],
        update_time: Optional[float],
        chunksize: Optional[int],
        num_chunk_divs: int,
        permute_evaluation: bool
) -> Tuple[int, tuple, dict, Union[Iterable, Iterator], tuple, int, str, float, bool, int]:
    """
    Get / manipulate input arguments into final forms used by pmap
    Args:
        pool: multiprocessing.Pool
        args: list or tuple (Default=None)
            Additional args to pass to func. Constant across tasks:
            results = [func(task, *args) for task in tasks]
        kwargs: dict (Default=None)
            kwargs to pass to func. Constant across tasks:
            results = [func(task, **kwargs) for task in tasks]
        tasks: list, tuple, or iterator yielding data
            each task is the arguments that will be passed to func
        task_sizes: list or tuple (Default=[])
            Task sizes for displaying progress. If empty, each task will be
            assigned size 1
        num_tasks: int (Default=None)
            Number of tasks, used for displaying progress.
            If num_tasks is None, it will be set to len(tasks) if available.
        update_time: float (Default=0.5)
            Minimum time in seconds for progress bar updates
        chunksize: int (Default=None)
            chunksize to pass to imap_unordered (number of tasks to pass to each
            worker in a block). If an integer, use that chunksize, if None, then
            attempt to pick an efficient number (num_chunk_divs chunks per
            worker, or chunksize = 1, whichever is larger)
        num_chunk_divs: int (Default=20)
            number of chunks per worker
        permute_evaluation: bool (Default=False)
            if True, pass tasks to pool in permuted order. This can help
                with load balancing to make evaluation time estimate more
                accurate.
                NOTE: if ordered is True, results will still
                      correspond to the tasks that are passed in.
                ALSO: if tasks were a generator, they will be stored in a list
                      for permutation. This may affect compute / memory usage.
            if False, pass each task to pool in order of tasks collection
    Returns:
        num_workers: int
            number of workers in pool
        args: list or tuple (Default=None)
            Convert to tuple if was None
        kwargs: dict (Default=None)
            Convert to dict if was None
        tasks: list, tuple, or iterator yielding data
            each task is the arguments that will be passed to func
        task_sizes: tuple[float]
            Convert to tuple
        num_tasks: int
            Calculate if was None and enough info is present:
                Will be assigned non-None if len(tasks) exists or permute == True
                or task_sizes exists.
        num_tasks_str: str
            A description of num_tasks (specifying "unknown number" if num_tasks
            is None)
        total: float
            if task_sizes is non-empty: sum(task_sizes)
            otherwise, num_tasks
        disable: bool
            Whether to use the progress meter.
        chunksize: int:
            Final chunksize to use.
    """
    if task_sizes is None:
        task_sizes = tuple()
    if args is None:
        args = tuple()
    if kwargs is None:
        kwargs = dict()
    # get parameters for progress meter
    disable = update_time is None or update_time >= float('inf')
    tasks = enumerate(tasks)
    if permute_evaluation:
        tasks = list(tasks)
        random.shuffle(tasks)
        if not disable and num_tasks is None:
            num_tasks = len(tasks)
    if task_sizes:
        task_sizes = tuple(task_sizes)
        total = sum(task_sizes)
        if num_tasks is None:
            num_tasks = len(task_sizes)
    else:
        if hasattr(tasks, '__len__'):
            num_tasks = len(tasks)
        total = num_tasks

    # noinspection PyProtectedMember
    num_workers = 1 if pool is None else len(pool._pool)
    if num_tasks is None:
        num_tasks_str = 'an unknown number of'
        if chunksize is None or chunksize < 1:
            chunksize = 1
    else:
        num_tasks_str = '%d' % num_tasks
        if chunksize is None or chunksize < 1:
            if task_sizes is not None and task_sizes:
                size_ratio = max(task_sizes) / min(s for s in task_sizes if s > 0)
                if numpy.isfinite(size_ratio) and size_ratio > num_chunk_divs:
                    num_chunk_divs = int(round(size_ratio))

            chunksize = max(num_tasks // (num_chunk_divs * num_workers), 1)
    assert type(chunksize) is int and chunksize > 0, \
        'chunksize must be a positive integer'

    return num_workers, args, kwargs, tasks, task_sizes, num_tasks, num_tasks_str, total, disable, chunksize


class _UnorderedMapTranslator(object):
    """
    Class for translating functions to work with pmap
    """
    def __init__(self, func, args=None, kwargs=None, starmap=False):
        if args is None:
            args = []
        if kwargs is None:
            kwargs = {}
        self._starmap = starmap
        self._master_pid = os.getpid()
        self._pickle_args(func, args, kwargs)

    def __call__(self, tup):
        if self._key not in _pool_func_args:
            self._load_func_args()
        func, args, kwargs = _pool_func_args[self._key]

        if self._starmap:
            return tup[0], func(*tup[1], *args, **kwargs)
        else:
            return tup[0], func(tup[1], *args, **kwargs)

    def __del__(self):
        if os.getpid() == self._master_pid:
            # note: unfortunately workers can't free up memory because they are
            # constantly being recycled, so they never know when the memory is
            # no longer needed. Fortunately the re-use of keys means that a new
            # job will overwrite the old memory. And closing the pool will free
            # all the workers memory. So this is not too bad.
            global _pool_func_args
            if self._key in _pool_func_args:
                _pool_func_args.pop(self._key)
            # master removes _key_file
            try:
                os.remove(self._key_file())
            except FileNotFoundError:
                pass

    def _pickle_args(self, func, args, kwargs):
        self._key = ''.join(
            random.SystemRandom().choices(
                string.ascii_letters + string.digits, k=_pool_key_len
            )
        )
        global _pool_func_args
        _pool_func_args[self._key] = (func, args, kwargs)
        with open(self._key_file(), 'wb') as f_out:
            dill.dump((func, args, kwargs), f_out)

    def _load_func_args(self):
        with open(self._key_file(), 'rb') as f_in:
            func, args, kwargs = dill.load(f_in)
        global _pool_func_args
        _pool_func_args[self._key] = (func, args, kwargs)
        self._is_master = False

    def _key_file(self):
        return os.path.join(tempfile.gettempdir(), self._key)


Pool = multiprocessing.get_context('spawn').Pool
multiprocessing.pool.Pool.pmap = pmap


def parmap(
        func: FuncType,
        tasks: Union[Iterable[TaskType], Iterator[TaskType]],
        num_tasks: Optional[int] = None,
        task_sizes: Iterable[float] = None,
        starmap: bool = False,
        flatmap: bool = False,
        ordered: bool = True,
        permute_evaluation: bool = False,
        description: str = 'parallel tasks',
        update_time: Optional[float] = Default.update_time,
        args: Collection = None, kwargs: dict = None,
        chunksize: Optional[int] = None,
        num_chunk_divs: int = Default.num_chunk_divs,
        n_jobs: int = Default.n_jobs,
        required_worker_memory: Optional[Numeric] = None,
        required_master_memory: Optional[Numeric] = None,
        require_physical_cpus: bool = False
) -> Iterator[OutputType]:
    """
    execute func on each task in tasks, using parallel pool
    progress is displayed using tqdm module. helper module to avoid need to manually create pool before calling pmap
    Typical invocation:
      results = [r for r in parmap(func, tasks)]
      assert results == [func(task) for task in tasks]
    INPUTS
        pool: multiprocessing.Pool
        func: function to execute on tasks
        tasks: list, tuple, or iterator yielding data
            each task is the arguments that will be passed to func
        num_tasks: int (Default=None)
            Number of tasks, used for displaying progress.
            If num_tasks is None, it will be set to len(tasks) if available.
        task_sizes: list or tuple (Default=[])
            Task sizes for displaying progress. If empty, each task will be
            assigned size 1
        starmap: bool (Default=False)
            if True, evaluate func(*task)  (i.e. treat task as an iterable)
            if False, evalute func(task)  (i.e. treat task as a single arg)
        flatmap: bool (Default=False)
            if True, treat output as collection of iterables and yield
                successively from each iterable. e.g.
                for iterable in result:
                    for val in iterable:
                        yield val
            if False, just yield results
        ordered: bool (Default=True)
            if True, yield results in order corresponding to tasks
                (i.e. results[n] = func(tasks[n]))
            if False, yield results in order they are completed by parallel pool
        permute_evaluation: bool (Default=False)
            if True, pass tasks to pool in permuted order. This can help
                with load balancing to make evaluation time estimate more
                accurate.
                NOTE: if ordered is True, results will still
                      correspond to the tasks that are passed in.
                ALSO: if tasks were a generator, they will be stored in a list
                      for permutation. This may affect compute / memory usage.
            if False, pass each task to pool in order of tasks collection
        description: str (Default='parallel tasks')
            Descriptive text for progress bar
        update_time: float (Default=0.5)
            Minimum time in seconds for progress bar updates.  If None or inf, progress will not be displayed.
        args: list or tuple (Default=None)
            Additional args to pass to func. Constant across tasks:
            results = [func(task, *args) for task in tasks]
        kwargs: dict (Default=None)
            kwargs to pass to func. Constant across tasks:
            results = [func(task, **kwargs) for task in tasks]
        chunksize: int (Default=None)
            chunksize to pass to imap_unordered (number of tasks to pass to each
            worker in a block). If an integer, use that chunksize, if None, then
            attempt to pick an efficient number (num_chunk_divs chunks per
            worker, or chunksize = 1, whichever is larger)
        num_chunk_divs: int (Default=20)
            number of chunks per worker
        required_worker_memory: float (Default=0)
            Used to constrain number of parallel workers based on available
            memory. Amount of memory needed by each worker process.
        required_master_memory: float (Default=0)
            Used to constrain number of parallel workers based on available
            memory. Amount of memory needed by master process.
        require_physical_cpus: bool (Default=False)
            If True, limit number of jobs to number of physical cores in system, not number of hyperthreads.
    OUTPUT:
        results will be yielded either in task or evaluation order, as
        specified specified by the "ordered" keyword.
    """
    n_jobs = common.num_jobs_to_use(
        n_jobs,
        required_master_memory=required_master_memory,
        required_worker_memory=required_worker_memory,
        require_physical_cpus=require_physical_cpus
    )
    if num_tasks is None:
        if isinstance(tasks, Sized):
            num_tasks = len(tasks)
        elif task_sizes is not None and isinstance(task_sizes, Sized):
            num_tasks = len(task_sizes)
    if num_tasks is not None:
        n_jobs = min(n_jobs, num_tasks)

    if n_jobs <= 1:
        for result in pmap(
            pool=None, func=func, tasks=tasks, num_tasks=num_tasks,
            task_sizes=task_sizes, starmap=starmap, flatmap=flatmap,
            ordered=ordered, permute_evaluation=permute_evaluation,
            description=description, update_time=update_time, args=args,
            kwargs=kwargs, chunksize=chunksize, num_chunk_divs=num_chunk_divs
        ):
            yield result
    else:
        # syntactic sugar causes problems with pycov
        # with Pool(processes=n_jobs) as pool:
        pool = Pool(processes=n_jobs)
        try:
            for result in pmap(
                pool=pool, func=func, tasks=tasks, num_tasks=num_tasks,
                task_sizes=task_sizes, starmap=starmap, flatmap=flatmap,
                ordered=ordered, permute_evaluation=permute_evaluation,
                description=description, update_time=update_time, args=args,
                kwargs=kwargs, chunksize=chunksize, num_chunk_divs=num_chunk_divs
            ):
                yield result
        finally:
            pool.close()
            pool.join()
