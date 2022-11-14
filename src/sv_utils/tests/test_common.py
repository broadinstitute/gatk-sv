import os
import numpy
import string
import pytest
import warnings

from sv_utils import common


@pytest.fixture(scope='function')
def tmpfile(tmpdir, filename=''):
    """
    yield temporary file for testing that will be deleted when the test function
    ends
    """
    if not filename:
        filename = 'tmp_' + ''.join(
            numpy.random.choice(list(string.ascii_letters), size=(6,))
        )
    filename = os.path.join(tmpdir, filename)
    yield filename
    if os.path.isfile(filename):
        os.remove(filename)


def _random_text(text_len: int = 5):
    return ''.join(numpy.random.choice(list(string.ascii_letters), size=text_len))


def test_deref_iterable_with_sorted_indices(num_values_test: int = 100, num_indices_test: int = 10):
    values = [_random_text() for _ in range(num_values_test)]
    sorted_indices = numpy.sort(numpy.random.permutation(num_values_test)[:num_indices_test])
    val_iter = (v for v in values)
    derefed_iter = common.deref_iterable_with_sorted_indices(val_iter, sorted_indices)
    derefed_arr = [values[i] for i in sorted_indices]
    assert numpy.array_equal(list(derefed_iter), derefed_arr)


def test_true_and_false(num_bools: int = 100):
    true_arr = common.true(num_bools)
    assert len(true_arr) == num_bools
    assert true_arr.all()
    false_arr = common.false(num_bools)
    assert len(false_arr) == num_bools
    assert (~false_arr).all()


def test_peekable_iter(list_len=5, data_len=1000, num_checks=10):
    peekable_iter = common.PeekableIter(numpy.arange(list_len))
    n = -1
    while peekable_iter.has_next():
        n += 1
        v_peek = peekable_iter.peek_next()
        assert v_peek == n
        v = next(peekable_iter)
        assert v == n
        peekable_iter.put_back(v)
        forward_vals = list(peekable_iter.peek())
        assert numpy.array_equal(forward_vals, numpy.arange(n, list_len))
        v2 = next(peekable_iter)
        assert v2 == n

    num_check_array = numpy.concatenate(
        ([0], numpy.logspace(0, numpy.log10(1000), num_checks - 1, dtype=int))
    )
    for num_check in num_check_array:
        raw_data = numpy.random.choice(list(string.printable), size=data_len)
        it = common.PeekableIter(raw_data)
        if num_check > 0:
            peak_it = it.peek()
            front = [next(peak_it) for _ in range(num_check)]
            assert front == raw_data[:num_check].tolist(), \
                "Peeked data doesn't match original data"
        assert list(it) == raw_data.tolist(), \
            "Data from PeekableIter doesn't match original data"


def test_command_results():
    # other tests will test valid commands where the output of the command is
    # expected. Here, just ensure that bad commands produce errors and good
    # commands produce non-empty
    bad_command = "Presumably this isn't a valid command"
    with pytest.raises(RuntimeError):
        common.command_results(bad_command)
    with warnings.catch_warnings():
        # this should just throw a warning and produce no result
        warnings.simplefilter("ignore")
        assert len(common.command_results(bad_command, exception_on_stderr=False)) == 0
    tests_path = os.path.dirname(os.path.realpath(__file__))
    ls_results = common.command_results(f"ls {tests_path}")
    assert os.path.basename(__file__) in ls_results


def test_add_exception_context():
    # Test exceptions with 0, 1, and 2 arguments passed to them, to make sure that the context is added sensibly
    # Try a few error types. No reason to believe any error type will cause problems, but just in case ...
    with pytest.raises(RuntimeError) as exception_info:
        # check it's still a RuntimeError
        try:
            raise RuntimeError()
        except RuntimeError as error:
            common.add_exception_context(error, "context")
            raise
    assert exception_info.value.args == ("context",)

    with pytest.raises(ValueError) as exception_info:
        # check it's still a ValueError
        try:
            raise ValueError("Foo")
        except ValueError as error:
            common.add_exception_context(error, "context")
            raise
    assert exception_info.value.args == ("context", "Foo")

    with pytest.raises(IOError) as exception_info:
        # check it's still an IOError
        try:
            raise IOError("Foo", "Bar")
        except IOError as error:
            common.add_exception_context(error, "context")
            raise
    assert exception_info.value.args == ("context", "Foo", "Bar")


def test_static_vars():
    @common.static_vars(my_static_var=42)
    def _func_with_static(new_val: float):
        old_val = _func_with_static.my_static_var
        _func_with_static.my_static_var = new_val
        return old_val

    # check that the initial stored value is good
    assert _func_with_static(3.14159) == 42
    # check that the value changed
    assert _func_with_static(0) == 3.14159


def test_classproperty():
    class Foo:
        def __init__(self, foo):
            self.foo = foo

        # noinspection PyMethodParameters
        @common.classproperty
        def default_foo(cls: type):
            return cls(42)

    assert Foo.default_foo.foo == 42


def test_dynamic_import():
    # test that dynamic_import can import a function, as usual
    acos_func = common.dynamic_import("numpy.arccos")
    pi_1 = acos_func(-1)
    # test by specifying base as a string
    pi_2 = common.dynamic_import("pi", base="numpy")
    assert pi_1 == pi_2
    # test by specifying the base as a module
    acos_func_2 = common.dynamic_import("arccos", base=numpy)
    assert acos_func(0.5) == acos_func_2(0.5)
    # finally ensure that if you try to dynamically import nonsense, you get a sensible error
    with pytest.raises(ModuleNotFoundError):
        common.dynamic_import("numpy.foo")


def test_cpu_and_memory_info():
    # hard to "test" this. Just make sure it returns stuff that isn't obviously crazy
    num_cores, num_hyperthreads = common.get_cpuinfo()
    assert num_hyperthreads >= num_cores > 0
    mem_bytes = common.available_memory(return_gb=False)
    mem_gb = common.available_memory(return_gb=True)
    assert mem_bytes > mem_gb > 0
    assert common.num_jobs_to_use() > 0
