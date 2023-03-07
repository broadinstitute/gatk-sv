import tarfile
from pathlib import Path
import json
import pytest
import numpy
import pandas
from sv_utils import common, genomics_io
from gq_recalibrator import tarred_properties_to_parquet
import common_test_utils

"""
This test suite looks really complicated, so a brief explanation is in order.
The basic idea is
1) To create a random data set that could represent training data for actual SVs.
   We do this instead of storing one in a testing folder because it may be large, and the
   requirements to test may change over time. It's a pain to create and maintain this by hand, and
   it stresses git.
2) Having created the random data set, we create corresponding TSVs
3) Then use the tarred_properties_to_parquet program to convert back to properties.
4) Finally we check that the loaded properties are equivalent to the randomly-generated ones

We do 1-3 with fixtures that can be re-used for other tests (e.g training filters). Hence most of
this test suite is fixtures.

Some notes:
* In making this test suite I encountered bugs related to there only being one parquet file in the
  tarred parquet data set. Typical data sets will have multiple parquet files. It's therefore
  necessary to test both cases.
* Opening the tarred gzipped TSVs could cause problems if they were already decompressed, so that's
  a second reason that they must be opened twice (to test that it works)
"""


class Default:
    random_seed = 0
    num_samples: int = 10
    num_variants: int = 50
    num_contigs: int = 2
    max_int_training_value = 100000
    min_int_training_value = -max_int_training_value
    num_variant_training_properties_per_type: int = 1
    num_sample_training_properties_per_type: int = 1
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
    property_type = tarred_properties_to_parquet.Keys.property_type
    num_rows = tarred_properties_to_parquet.Keys.num_rows
    num_samples = tarred_properties_to_parquet.Keys.num_samples
    codes = tarred_properties_to_parquet.Keys.codes


PropertyType = tarred_properties_to_parquet.PropertyType
PropertySummary = tarred_properties_to_parquet.PropertySummary
PropertiesSummary = tarred_properties_to_parquet.PropertiesSummary


@pytest.fixture(scope="session")
def random_generator(
        random_seed: int = Default.random_seed
) -> numpy.random.Generator:
    yield common.init_generator(generator_init=random_seed)


@pytest.fixture(scope="session")
def sample_ids(num_samples: int = Default.num_samples) -> tuple[str]:
    yield tuple(f"sample_{sample_num}" for sample_num in range(num_samples))


@pytest.fixture(scope="session")
def properties_dataframe_and_summary(
        random_generator: numpy.random.Generator,
        sample_ids: tuple[str],
        num_variants: int = Default.num_variants,
        num_variant_training_properties_per_type: int
        = Default.num_variant_training_properties_per_type,
        num_sample_training_properties_per_type: int
        = Default.num_sample_training_properties_per_type,
        num_contigs: int = Default.num_contigs
) -> (pandas.DataFrame, PropertiesSummary):
    variant_locs, properties_summary = _generate_random_sv_locations(
        random_generator=random_generator, num_variants=num_variants, num_contigs=num_contigs
    )

    variant_training_properties = _generate_random_training_properties(
        random_generator=random_generator, num_variants=num_variants,
        training_properties_start_index=0,
        num_training_properties_per_type=num_variant_training_properties_per_type,
        sample_id=None, num_samples=0, properties_summary=properties_summary
    )

    properties_dataframe = pandas.concat(
        (variant_locs, variant_training_properties) + tuple(
            _generate_random_training_properties(
                random_generator=random_generator, num_variants=num_variants,
                training_properties_start_index=len(variant_training_properties.columns),
                num_training_properties_per_type=num_sample_training_properties_per_type,
                sample_id=sample_id, num_samples=len(sample_ids),
                properties_summary=properties_summary
            ) for sample_id in sample_ids
        ),
        axis=1
    )
    return properties_dataframe, properties_summary


@pytest.fixture(scope="session")
def properties_dataframe(
    properties_dataframe_and_summary: tuple[pandas.DataFrame, PropertiesSummary]
) -> pandas.DataFrame:
    return properties_dataframe_and_summary[0]


@pytest.fixture(scope="session")
def one_hot_properties_dataframe(
        properties_dataframe_and_summary: tuple[pandas.DataFrame, PropertiesSummary]
) -> pandas.DataFrame:
    yield tarred_properties_to_parquet.categorical_properties_to_one_hot(
        properties_dataframe=properties_dataframe_and_summary[0],
        properties_summary=properties_dataframe_and_summary[1]
    )


@pytest.fixture(scope="session")
def properties_summary(
        properties_dataframe_and_summary: tuple[pandas.DataFrame, PropertiesSummary]
) -> PropertiesSummary:
    return properties_dataframe_and_summary[1]


def _generate_random_sv_locations(
        random_generator: numpy.random.Generator,
        num_variants: int,
        num_contigs: int
) -> (pandas.DataFrame, PropertiesSummary):
    # make initial properties_summary for dtype consistency, then add codes
    properties_summary = PropertiesSummary({
        Keys.variant_id: PropertySummary(
            property_type=PropertyType.String, num_rows=num_variants, num_samples=0
        ),
        Keys.contig: PropertySummary(
            property_type=PropertyType.String, num_rows=num_variants, num_samples=0,
            codes=[f"chr{n}" for n in range(num_contigs)]
        ),
        Keys.begin: PropertySummary(
            property_type=PropertyType.int, num_rows=num_variants, num_samples=0
        ),
        Keys.svlen: PropertySummary(
            property_type=PropertyType.int, num_rows=num_variants, num_samples=0
        ),
        Keys.svtype: PropertySummary(
            property_type=PropertyType.String, num_rows=num_variants, num_samples=0,
            codes=["INS", "DEL", "DUP", "INV"]
        ),
        Keys.end: PropertySummary(
            property_type=PropertyType.int, num_rows=num_variants, num_samples=0
        )
    })
    # generate independent properties
    base_properties = pandas.DataFrame(
        {
            Keys.variant_id: [f"variant_{n}" for n in range(num_variants)],
            Keys.contig: random_generator.choice(
                properties_summary[Keys.contig].codes, size=(num_variants,)
            ),
            Keys.begin: random_generator.integers(
                low=1, high=10**6, size=(num_variants,),
                dtype=properties_summary[Keys.begin].dtype
            ),
            Keys.svlen:
                50 + random_generator.exponential(
                    scale=10000, size=(num_variants,),
                ).round().astype(properties_summary[Keys.svlen].dtype),
            Keys.svtype: random_generator.choice(
                properties_summary[Keys.svtype].codes, size=(num_variants,)
            )
        }
    )
    # make id of type "string[pyarrow]"
    base_properties[Keys.variant_id] = base_properties[Keys.variant_id].astype("string[pyarrow]")
    # make contig and svtype categorical
    base_properties[Keys.contig] = base_properties[Keys.contig].astype(
        properties_summary[Keys.contig].dtype
    )
    base_properties[Keys.svtype] = base_properties[Keys.svtype].astype(
        properties_summary[Keys.svtype].dtype
    )
    # generate END depending on BEGIN and SVTYPE+SVLEN
    base_properties[Keys.end] = numpy.where(
        base_properties[Keys.svtype] == "INS",
        base_properties[Keys.begin],
        base_properties[Keys.begin] + base_properties[Keys.svlen]
    ).astype(properties_summary[Keys.end].dtype)
    # sort properties
    base_properties = genomics_io.sort_intervals_table(base_properties)
    # convert column names to multi-sample format
    # noinspection PyTypeChecker
    base_properties.columns = pandas.MultiIndex.from_tuples(
        [(None, column) for column in base_properties.columns],
        names=Default.column_levels
    )

    return base_properties, properties_summary


def _generate_random_training_properties(
        random_generator: numpy.random.Generator,
        num_variants: int,
        training_properties_start_index: int,
        num_training_properties_per_type: int,
        sample_id: str | None,
        num_samples: int,
        properties_summary: PropertiesSummary
) -> pandas.DataFrame:
    training_properties = pandas.DataFrame(
        {
            (sample_id, property_name): _gen_random_property_values(
                random_generator=random_generator,
                property_name=property_name,
                property_type=property_type,
                num_values=num_variants,
                num_samples=num_samples,
                properties_summary=properties_summary
            )
            for property_name, property_type in
            (
                (f"train-{property_type.name}-{property_num}", property_type)
                for property_num, property_type in enumerate(
                    (
                        property_type
                        for property_type in PropertyType
                        for __ in range(num_training_properties_per_type)
                    ),
                    start=training_properties_start_index
                )
            )
        }
    )

    # convert column names to multi-sample format
    training_properties.columns = pandas.MultiIndex.from_tuples(
        training_properties.columns, names=Default.column_levels
    )
    return training_properties


def _gen_random_property_values(
        random_generator: numpy.random.Generator,
        property_name: str,
        property_type: PropertyType,
        num_values: int,
        num_samples: int,
        properties_summary: PropertiesSummary
) -> pandas.Series:
    if property_name in properties_summary:
        property_summary = properties_summary[property_name]
    else:
        # need to add new PropertySummary
        codes = _get_random_codes(
            random_generator=random_generator, property_name=property_name,
            property_type=property_type
        ) if property_type.has_codes else []
        property_summary = PropertySummary(
            property_type=property_type, num_rows=num_values, num_samples=num_samples, codes=codes
        )
        properties_summary[property_name] = property_summary

    if property_type.is_floating_type:
        values = random_generator.normal(
            loc=0.0, scale=10.0, size=(num_values,)
        ).astype(property_summary.dtype)
    elif property_summary.is_integer_type:
        iinfo = numpy.iinfo(property_summary.dtype)
        values = random_generator.integers(
            max(-10**6, iinfo.min), min(10**6, iinfo.max), size=(num_values,),
            dtype=property_summary.dtype
        )
    elif property_summary.is_bool_type:
        values = random_generator.integers(0, 2, (num_values,), dtype=bool)
    else:
        values = random_generator.choice(property_summary.codes, size=(num_values,), replace=True)
    return pandas.Series(values, name=property_name, dtype=property_summary.dtype)


def _get_random_codes(
        random_generator: numpy.random.Generator, property_name: str, property_type: PropertyType
) -> list[str]:
    num_codes = random_generator.integers(2, 6)
    # codes don't have to be sorted, make them randomly ordered
    codes = [f"{property_name}={code_index}"
             for code_index in random_generator.permutation(range(num_codes))]
    if property_type == PropertyType.String:
        return codes
    # we're making a StringSet property, so actual code dict is some observed subset of possible
    # combinations of these base codes, specified as a comma-separated list
    combined_codes = []
    while len(combined_codes) < min(2 ** num_codes - 1, 10):
        num_combine = random_generator.integers(0, num_codes)
        new_code = ','.join(
            sorted(
                random_generator.choice(codes, num_combine, replace=False)
            )
        )
        if new_code not in combined_codes:
            combined_codes.append(new_code)
    return combined_codes


@pytest.fixture(scope="session")
def temp_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """ explicitly convert to Path """
    yield tmp_path_factory.mktemp("test_tarred_properties_to_parquet")


@pytest.fixture(scope="session")
def tarred_properties_to_parquet_dir(temp_path: Path) -> Path:
    """ Yield consistent temp folder for files related to these tests """
    properties_dir = temp_path / "test_tarred_properties_to_parquet"
    properties_dir.mkdir(parents=True, exist_ok=True)
    yield properties_dir


class PropertySummaryEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, PropertySummary):
            return {
                getattr(Keys, _slot): getattr(obj, _slot) for _slot in obj.__slots__
                if _slot != "__weakref__"
            }
        elif isinstance(obj, PropertyType):
            return str(obj.value)
        return json.JSONEncoder.default(self, obj)


@pytest.fixture(scope="session")
def properties_tsvs_folder(
        tarred_properties_to_parquet_dir,
        properties_dataframe: pandas.DataFrame,
        properties_summary: PropertiesSummary,
        sample_ids: tuple[str]
) -> Path:
    """ Yield consistent temp folder for files related to these tests """
    tsvs_folder = tarred_properties_to_parquet_dir / "properties"
    tsvs_folder.mkdir()
    property_names = set(properties_dataframe.columns.get_level_values(level=Keys.property))
    for property_name in property_names:
        property_values = properties_dataframe.loc[:, (slice(None), property_name)]
        property_summary = properties_summary[property_name]
        if property_summary.codes:
            # the TSVs for encodings are the ordinal encoding according to the codes list, so need
            # to get those ordinal codes
            property_values = pandas.DataFrame(
                {
                    column: property_values[column].cat.codes.astype(
                        property_summary.encoding_dtype
                    )
                    for column in property_values.columns
                },
                columns=property_values.columns
            )
        elif property_summary.is_bool_type:
            property_values = pandas.DataFrame(
                {
                    column: property_values[column].astype(property_summary.encoding_dtype)
                    for column in property_values.columns
                },
                columns=property_values.columns
            )

        # pandas won't write multi-index TSV even if not writing header / index. so manually
        # flatten columns:
        tarred_properties_to_parquet.flatten_columns(property_values)
        tsv_file = tsvs_folder / f"{property_name}.tsv.gz"
        genomics_io.pandas_to_tsv(tsv_file, property_values, write_index=False, write_header=False)

    sample_ids_file = tsvs_folder / "sample_ids.list"
    with open(sample_ids_file, 'w') as f_out:
        for sample_id in sample_ids:
            f_out.write(f"{sample_id}\n")

    properties_summary_file = tsvs_folder / "properties_summary.json"
    with open(properties_summary_file, 'w') as f_out:
        json.dump(properties_summary, f_out, cls=PropertySummaryEncoder)

    yield tsvs_folder


@pytest.fixture(scope="session")
def properties_tsvs_tar_file(
        properties_tsvs_folder: Path
) -> Path:
    """
    yield parallel pool for executing tests
    """
    yield tarred_properties_to_parquet.archive_to_tar(
        folder_to_archive=properties_tsvs_folder, remove_files=True
    )


@pytest.fixture(scope="session")
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
    with tarfile.TarFile(tarred_parquet_file, 'r') as tar_in:
        assert len(tar_in.getnames()) == 1, "Expected only one parquet file in this archive"
    yield tarred_parquet_file, scaling_json


@pytest.fixture(scope="session")
def properties_parquet_files_multiple_partitions(
        properties_tsvs_tar_file: Path
) -> (Path, Path):
    tarred_parquet_file = properties_tsvs_tar_file.with_suffix(".multiple.pq.tar")
    scaling_json = properties_tsvs_tar_file.parent.joinpath("multiple-property-scales.json")
    tarred_properties_to_parquet.tarred_properties_to_parquet(
        input_tar=properties_tsvs_tar_file,
        output_path=tarred_parquet_file,
        properties_scaling_json=scaling_json,
        remove_input_tar=False,
        dask_partition_size_mb=0.001
    )
    with tarfile.TarFile(tarred_parquet_file, 'r') as tar_in:
        assert len(tar_in.getnames()) > 1, "Expected multiple parquet files in this archive"
    yield tarred_parquet_file, scaling_json


@pytest.fixture(scope="session")
def properties_tarred_parquet_file(
        properties_parquet_files: tuple[Path, Path]
) -> Path:
    yield properties_parquet_files[0]


@pytest.fixture(scope="session")
def properties_scaling_json(
        properties_parquet_files: tuple[Path, Path]
) -> Path:
    yield properties_parquet_files[1]


def test_conversion_faithful(
        one_hot_properties_dataframe: pandas.DataFrame,
        properties_tarred_parquet_file: Path
):
    # load as parquet files as dask DataFrame and then immediately compute (convert to pandas)
    # because this is small and can easily fit into memory
    loaded_properties_dataframe = tarred_properties_to_parquet.parquet_to_df(
        input_path=properties_tarred_parquet_file, remove_input_tar=False
    ).compute()
    common_test_utils.assert_dataframes_equal(
        one_hot_properties_dataframe, loaded_properties_dataframe,
        "original properties_dataframe != loaded_properties_dataframe",
        check_column_order=False
    )


def test_multiple_partitions_equivalent_to_one(
        one_hot_properties_dataframe: pandas.DataFrame,
        properties_scaling_json: Path,
        properties_parquet_files_multiple_partitions: tuple[Path, Path]
):
    multiple_partitions_tarred_parquet, multiple_partitions_json \
        = properties_parquet_files_multiple_partitions
    with open(properties_scaling_json, 'r') as f_in:
        single_json_obj = json.load(f_in)
    with open(multiple_partitions_json, 'r') as f_in:
        multiple_json_obj = json.load(f_in)
    # floating baseline and scale use median-of-medians, so don't expect exact agreement,
    # especially on such a small data set
    # 1) all keys should be the same
    assert single_json_obj.keys() == multiple_json_obj.keys()
    # 2) all non-float key,value pairs should be the same
    assert {
        key: value for key, value in single_json_obj.items()
        if not any(f"train-{_t}-" in key for _t in ("float", "double"))
    } == {
        key: value for key, value in multiple_json_obj.items()
        if not any(f"train-{_t}-" in key for _t in ("float", "double"))
    }

    # both dataframes should be the same when computed
    multiple_dataframe = tarred_properties_to_parquet.parquet_to_df(
        input_path=multiple_partitions_tarred_parquet, remove_input_tar=False
    ).compute()
    common_test_utils.assert_dataframes_equal(
        one_hot_properties_dataframe, multiple_dataframe,
        "single-partition properties_dataframe != multiple-partition properties_dataframe",
        check_column_order=False
    )

