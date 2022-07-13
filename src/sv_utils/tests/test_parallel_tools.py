#!/usr/bin/env python

import numpy
import numpy.random
import time
import pytest
import random
import string

from typing import Iterator, Tuple, Union

from sv_utils import common, parallel_tools


BasicTask = numpy.ndarray
StarTask = Tuple[BasicTask, str]
Task = Union[BasicTask, StarTask]


@pytest.fixture(scope='session')
def pool(n_jobs=common.num_physical_cpus):
    """
    yield parallel pool for executing tests
    """
    with parallel_tools.Pool(processes=n_jobs) as par_pool:
        yield par_pool


def _make_star_tasks(num_tasks: int = 1000) -> Iterator[StarTask]:
    """
    generate tasks for demo that are inherently multi-argument
    """
    chars = string.ascii_letters + string.digits
    num_chars = len(chars)
    for task_num in range(num_tasks):
        task_name = chars[task_num % num_chars]
        tot = task_num // num_chars
        while tot:
            task_name += chars[tot % num_chars]
            tot //= num_chars
        num_rand = numpy.random.randint(3, 10)
        yield numpy.random.randn(num_rand), task_name


def _make_basic_tasks() -> Iterator[BasicTask]:
    """
    generate tasks for demo that are single-argument
    """
    for vals, name in _make_star_tasks():
        yield vals


def _get_task_size(task: Task) -> float:
    """
    Return estimate of task difficulty
    """
    if isinstance(task, Tuple):
        vals = task[0]  # it's a StarTask
    else:
        vals = task  # it's a BasicTask
    return round(sum(abs(vals)), ndigits=6)


def _star_func(v, k="no star", extra_val_1="extra 1", extra_val_2="extra 2", sleep_time=None) -> str:
    """
    evaluate input values returning string
    principle virtues of this function are:
     1) simple
     2) can be a star function or a single-valued function
     3) fairly unique mapping from inputs to outputs
     4) if sleep_time is not None, should be possible to measure
        speed-up from parallel execution
    """
    if sleep_time is not None:
        if not sleep_time:
            sleep_time = 1.0
        sleep_time *= sum(abs(v))
        time.sleep(sleep_time)
    return "%s, %s, %s : %f" % (k, extra_val_1, extra_val_2, sum(abs(v)))


def _is_multiple_of_5(int_list: numpy.ndarray) -> numpy.ndarray:
    return int_list.compress(int_list % 5 == 0)


def _run_accuracy_trial(pool, starmap=False, args=None, kwargs=None,
                        permute_evaluation=False, ordered=True,
                        chunksize=None):
    """
    For given set of input parameters test that parallel execution
    returns correct results by comparing to serial
    """
    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}
    # get tasks and serial results
    if starmap:
        tasks = list(_make_star_tasks())
        serial_results = \
            [_star_func(*task, *args, **kwargs) for task in tasks]
    else:
        tasks = list(_make_basic_tasks())
        serial_results = \
            [_star_func(task, *args, **kwargs) for task in tasks]
    # get parallel results
    parallel_results = list(parallel_tools.pmap(
        pool, _star_func, tasks, args=args, kwargs=kwargs, update_time=None,
        starmap=starmap, ordered=ordered, chunksize=chunksize,
        permute_evaluation=permute_evaluation
    ))
    # check equality
    if ordered:
        assert parallel_results == serial_results
    else:
        assert sorted(parallel_results) == sorted(serial_results)


def _random_word(word_len=10):
    """
    Generate a random "word" of printable characters
    """
    return ''.join(random.choices(string.printable, k=word_len))


def test_accuracy(pool):
    """
    Test accuracy of parallel execution under a variety of use cases
    """
    # test basic imap functionality
    _run_accuracy_trial(pool)

    # test imap with extra positional args
    args = [_random_word() for __ in range(3)]
    _run_accuracy_trial(pool, args=args)

    # test imap with extra keyword args
    kwargs = {'extra_val_2': _random_word()}
    _run_accuracy_trial(pool, kwargs=kwargs)

    # test imap with extra args and keyword args
    args = [_random_word()]
    kwargs = {'extra_val_2': _random_word()}
    _run_accuracy_trial(pool, args=args, kwargs=kwargs)

    # test basic starmap functionality
    _run_accuracy_trial(pool, starmap=True)

    # test starmap with extra positional args
    args = [_random_word() for __ in range(2)]
    _run_accuracy_trial(pool, starmap=True, args=args)

    # test starmap with extra keyword args
    kwargs = {'extra_val_2': _random_word()}
    _run_accuracy_trial(pool, starmap=True, kwargs=kwargs)

    # test starmap with extra args and keyword args
    args = [_random_word()]
    kwargs = {'extra_val_2': _random_word()}
    _run_accuracy_trial(pool, starmap=True, args=args, kwargs=kwargs)

    # test permuting evaluation order
    _run_accuracy_trial(pool, permute_evaluation=True)

    # test unordered results
    _run_accuracy_trial(pool, permute_evaluation=True, ordered=False)

    # test with chunksize > 1
    _run_accuracy_trial(pool, chunksize=10)

    # test with chunksize == 1
    _run_accuracy_trial(pool, chunksize=1)


def test_performance(pool, capsys, test_time=0.5):
    """
    Check that parallel execution is faster, wait bar can draw
    """
    # get tasks as generators to test that functionality
    tasks = list(_make_star_tasks())
    task_sizes = [_get_task_size(task) for task in tasks]

    num_workers = len(pool._pool)
    # make serial tasks go faster to avoid waiting a long time
    t_s = (test_time / 2.0) / sum(task_sizes)
    t_p = t_s * num_workers

    t0 = time.time()
    serial_results = [_star_func(*task, sleep_time=t_s) for task in tasks]
    t1 = time.time()
    with capsys.disabled():
        parallel_results = list(parallel_tools.pmap(
            pool, _star_func, tasks, task_sizes=task_sizes,
            starmap=True, kwargs={'sleep_time': t_p}, update_time=0.0,
            description='This should draw a progress bar'
        ))
        t2 = time.time()
        assert serial_results == parallel_results
    t_serial = t1 - t0
    t_parallel = t2 - t1
    speed_up = (t_serial / t_s) / (t_parallel / t_p)
    assert speed_up > 0.67 * num_workers


def test_parmap():
    tasks = list(_make_star_tasks())
    serial_results = [_star_func(*task) for task in tasks]
    parmap_serial_results = list(
        parallel_tools.parmap(_star_func, tasks, starmap=True, n_jobs=1)
    )
    assert parmap_serial_results == serial_results, 'parmap serial results don''t match serial'
    parallel_results = list(
        parallel_tools.parmap(_star_func, tasks, starmap=True)
    )
    assert parallel_results == serial_results, \
        'parmap results don''t match serial'


def test_flatmap(max_check=10000, num_tasks=16):
    serial_results = _is_multiple_of_5(numpy.arange(max_check))
    tasks = [
        numpy.arange(t * max_check // num_tasks, (t + 1) * max_check // num_tasks)
        for t in range(num_tasks)
    ]
    parallel_results = list(
        parallel_tools.parmap(
            _is_multiple_of_5, tasks=tasks, flatmap=True
        )
    )
    assert numpy.array_equal(serial_results, parallel_results), "parmap flatmap results don't match serial"

    parallel_results = list(
        parallel_tools.parmap(
            _is_multiple_of_5, tasks=tasks, flatmap=True, permute_evaluation=True
        )
    )
    assert numpy.array_equal(serial_results, parallel_results), \
        "parmap flatmap results don't match serial when permuting evaluation"


def test_empty():
    # check that empty tasks return empty results
    assert not list(parallel_tools.parmap(_star_func, [])), "empty tasks list returned non-empty results"
    assert not list(parallel_tools.parmap(_star_func, (_x for _x in range(0)))), \
        "empty tasks generator returned non-empty results"
    assert not list(parallel_tools.parmap(_star_func, [], flatmap=True)), "empty tasks list returned non-empty results"
    assert not list(parallel_tools.parmap(_star_func, (_x for _x in range(0)), flatmap=True)), \
        "empty tasks generator returned non-empty results"
