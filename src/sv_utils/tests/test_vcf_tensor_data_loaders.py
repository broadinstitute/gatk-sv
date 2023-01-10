import pytest
from pathlib import Path


class Default:
    num_samples: int = 100
    num_variants: int = 1000


@pytest.fixture(scope='session')
def test_parquet_file(tmpdir) -> Path:
    """
    yield parallel pool for executing tests
    """
    temp_out_dir = tmpdir.mkdir("test_vcf_tensor_data_loaders")

    yield parquet_file