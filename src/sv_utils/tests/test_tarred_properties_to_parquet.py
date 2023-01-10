import py
from pathlib import Path
import json
import pytest
import numpy
import pandas
from sv_utils import common, genomics_io
from gq_recalibrator import tarred_properties_to_parquet


class Default:
    random_seed = 0
    num_samples: int = 20
    num_variants: int = 100
    num_contigs: int = 2
    max_int_training_value = 100000
    min_int_training_value = -max_int_training_value
    num_int_training_properties: int = 2
    num_bool_training_properties: int = 2
    column_levels = tarred_properties_to_parquet.Default.column_levels


class Keys:
    variant_id = genomics_io.Keys.id
    contig = genomics_io.Keys.contig
    begin = genomics_io.Keys.begin
    end = genomics_io.Keys.end
    svlen = genomics_io.Keys.svlen
    svtype = genomics_io.Keys.svtype
    property = tarred_properties_to_parquet.Keys.property
    sample_id = tarred_properties_to_parquet.Keys.sample_id


PropertySummary = tarred_properties_to_parquet.PropertySummary
PropertiesSummary = tarred_properties_to_parquet.PropertiesSummary


@pytest.fixture(scope="function")
def random_generator(
        random_seed: int = Default.random_seed
) -> numpy.random.Generator:
    yield common.init_generator(generator_init=random_seed)


@pytest.fixture(scope="function")
def sample_ids(num_samples: int = Default.num_samples) -> tuple[str]:
    yield tuple(f"sample_{sample_num}" for sample_num in range(num_samples))


@pytest.fixture(scope="function")
def properties_dataframe(
        random_generator: numpy.random.Generator,
        sample_ids: tuple[str],
        num_variants: int = Default.num_variants,
        num_int_training_properties: int = Default.num_int_training_properties,
        num_bool_training_properties: int = Default.num_bool_training_properties,
        num_contigs: int = Default.num_contigs
) -> pandas.DataFrame:
    variant_locs = _generate_random_sv_locations(
        random_generator=random_generator, num_variants=num_variants, num_contigs=num_contigs
    )
    num_int_variant_properties = num_int_training_properties // 2
    num_bool_variant_properties = num_bool_training_properties // 2
    num_variant_properties = num_int_variant_properties + num_bool_variant_properties

    variant_training_properties = _generate_random_training_properties(
        random_generator=random_generator, num_variants=num_variants,
        int_training_properties_range=range(0, num_int_variant_properties),
        bool_training_properties_range=range(num_int_variant_properties, num_variant_properties),
        sample_id=None
    )

    num_int_sample_properties = num_int_training_properties - num_int_variant_properties
    end_int_sample_properties = num_variant_properties + num_int_sample_properties
    num_bool_sample_properties = num_bool_training_properties - num_bool_variant_properties
    end_bool_sample_properties = end_int_sample_properties + num_bool_sample_properties

    properties_dataframe = pandas.concat(
        (variant_locs, variant_training_properties) + tuple(
            _generate_random_training_properties(
                random_generator=random_generator, num_variants=num_variants,
                int_training_properties_range=range(num_variant_properties,
                                                    end_int_sample_properties),
                bool_training_properties_range=range(end_int_sample_properties,
                                                     end_bool_sample_properties),
                sample_id=sample_id
            ) for sample_id in sample_ids
        ),
        axis=1
    )
    return properties_dataframe


def _generate_random_sv_locations(
        random_generator: numpy.random.Generator,
        num_variants: int,
        num_contigs: int
) -> pandas.DataFrame:
    # generate independent properties
    base_properties = pandas.DataFrame(
        {
            Keys.variant_id: [f"variant_{n}" for n in range(num_variants)],
            Keys.contig: random_generator.choice(
                [f"chr{n}" for n in range(num_contigs)], size=(num_variants,)
            ),
            Keys.begin: random_generator.integers(low=1, high=10**6, size=(num_variants,)),
            Keys.svlen:
                50 + random_generator.exponential(scale=10000, size=(num_variants,))
                .round().astype(int),
            Keys.svtype: random_generator.choice(
                ["INS", "DEL", "DUP", "INV"], size=(num_variants,)
            )
        }
    )
    # generate END depending on BEGIN and SVTYPE+SVLEN
    base_properties[Keys.end] = numpy.where(
        base_properties[Keys.svtype] == "INS",
        base_properties[Keys.begin],
        base_properties[Keys.begin] + base_properties[Keys.svlen]
    )
    # sort properties
    base_properties = genomics_io.sort_intervals_table(base_properties)
    # convert column names to multi-sample format
    # noinspection PyTypeChecker
    base_properties.columns = pandas.MultiIndex.from_tuples(
        [(None, column) for column in base_properties.columns],
        names=Default.column_levels
    )
    return base_properties


def _generate_random_training_properties(
        random_generator: numpy.random.Generator,
        num_variants: int,
        int_training_properties_range: range,
        bool_training_properties_range: range,
        sample_id: str | None
) -> pandas.DataFrame:
    training_properties = pandas.DataFrame(
        {
            **{
                (sample_id, f"int_train_{property_num}"):
                    random_generator.integers(-10**6, 10**6, (num_variants,))
                for property_num in int_training_properties_range
            },
            **{
                (sample_id, f"bool_train_{property_num}"):
                    random_generator.integers(0, 2, (num_variants,), dtype=bool)
                for property_num in bool_training_properties_range
            }
        }
    )
    # convert column names to multi-sample format
    training_properties.columns = pandas.MultiIndex.from_tuples(
        training_properties.columns, names=Default.column_levels
    )
    return training_properties


# noinspection PyProtectedMember,PyUnresolvedReferences
@pytest.fixture(scope="function")
def temp_path(tmpdir: py._path.local.LocalPath) -> Path:
    """ explicitly convert to Path """
    yield Path(tmpdir)


@pytest.fixture(scope="function")
def tarred_properties_to_parquet_dir(temp_path: Path) -> Path:
    """ Yield consistent temp folder for files related to these tests """
    properties_dir = temp_path / "test_tarred_properties_to_parquet"
    properties_dir.mkdir(parents=True, exist_ok=True)
    yield properties_dir


class PropertySummaryEncoder:
    def default(self, obj):
        if isinstance(obj, tarred_properties_to_parquet.PropertySummary):
            return obj.__dict__
        return json.JSONEncoder.default(self, obj)


@pytest.fixture(scope="function")
def properties_tsvs_folder(
        tarred_properties_to_parquet_dir,
        properties_dataframe: pandas.DataFrame,
        sample_ids: tuple[str]
) -> Path:
    """ Yield consistent temp folder for files related to these tests """
    tsvs_folder = tarred_properties_to_parquet_dir / "properties"
    tsvs_folder.mkdir()
    property_names = set(properties_dataframe.columns.get_level_values(level=Keys.property))
    properties_summary: PropertiesSummary = PropertiesSummary({})
    for property_name in property_names:
        property_values = properties_dataframe.loc[:, (slice(None), property_name)]
        # pandas won't write multi-index TSV even if not writing header / index. so manually flatten
        # columns:
        tarred_properties_to_parquet.flatten_columns(property_values)
        tsv_file = tsvs_folder / property_name
        genomics_io.pandas_to_tsv(tsv_file, property_values, write_index=False, write_header=False)
        properties_summary[property_name] = PropertiesSummary(
            _get_property_summary(property_values)
        )

    sample_ids_file = tsvs_folder / "sample_ids.list"
    with open(sample_ids_file, 'w') as f_out:
        for sample_id in sample_ids:
            f_out.write(f"{sample_id}\n")

    properties_summary_file = tsvs_folder / "properties_summary.json"
    with open(properties_summary_file, 'w') as f_out:
        json.dump()

    yield tsvs_folder


@pytest.fixture(scope="function")
def properties_tsvs_tar_file(
        properties_tsvs_folder: Path
) -> Path:
    """
    yield parallel pool for executing tests
    """
    yield tarred_properties_to_parquet.archive_to_tar(
        folder_to_archive=properties_tsvs_folder, remove_files=False
    )


@pytest.fixture(scope="function")
def properties_parquet_files(
        properties_tsvs_tar_file: Path
) -> (Path, Path):
    tarred_parquet_file = properties_tsvs_tar_file.with_suffix(".pq.tar")
    scaling_json = properties_tsvs_tar_file.parent.joinpath("property-scales.json")
    tarred_properties_to_parquet.tarred_properties_to_parquet(
        input_tar=properties_tsvs_tar_file,
        output_path=tarred_parquet_file,
        properties_scaling_json=scaling_json,
        remove_input_tar=False
    )
    yield tarred_parquet_file, scaling_json


@pytest.fixture(scope="function")
def properties_tarred_parquet_file(
        properties_parquet_files: tuple[Path, Path]
) -> Path:
    yield properties_parquet_files[0]


@pytest.fixture(scope="function")
def properties_scaling_json(
        properties_parquet_files: tuple[Path, Path]
) -> Path:
    yield properties_parquet_files[1]


def test_conversion_faithful(
        properties_dataframe: pandas.DataFrame,
        properties_tarred_parquet_file: Path
):
    # load as parquet files as dask DataFrame and then immediately compute (convert to pandas)
    # because this is small and can easily fit into memory
    loaded_properties_dataframe = tarred_properties_to_parquet.parquet_to_df(
        input_path=properties_tarred_parquet_file, remove_input_tar=False
    ).compute()
    assert properties_dataframe.equals(loaded_properties_dataframe)
