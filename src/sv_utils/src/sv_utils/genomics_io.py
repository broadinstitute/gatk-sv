#!/usr/bin/env python
import io
import os
import sys
import ast
import types
import warnings
from enum import Enum
from types import MappingProxyType
from typing import Text, Union, Tuple, Mapping, Optional, Any, Dict, Callable, Sequence, List, Collection, ValuesView, \
    Iterator, Set, Iterable, TypeVar
import numpy
import pandas
import pandas.core.dtypes.base
import pandas.core.dtypes.cast
import pysam

from sv_utils import common

NA = pandas.NA  # missing field_type, allows for nullable integer
NAType = type(NA)
NoneType = type(None)
dtype = Union[numpy.dtype, pandas.core.dtypes.base.ExtensionDtype]
MissingType = Union[str, float, NAType]
EncodedVcfField = Union[int, float, str, NAType]
VcfEncodingMapper = Mapping[str, EncodedVcfField]
VariantPropertiesSource = Union[pysam.VariantRecord, pysam.libcbcf.VariantRecordInfo]
VcfPropSpecifier = Optional[Union[Dict[str, str], Sequence[str]]]
VariantPropertyExtractor = Callable[[pysam.VariantRecord], Optional[EncodedVcfField]]
GenotypePropertyExtractor = Callable[[pysam.libcbcf.VariantRecordSample], Optional[EncodedVcfField]]
Genotype = Optional[Tuple[Optional[int], ...]]


# monkey-patch warnings to just display the warning
# noinspection PyUnusedLocal
def _showwarning(message, category=UserWarning, filename='', lineno=-1, file='', line=''):
    message = f"\n{category.__name__}  {filename}:{lineno}  {message}" if lineno > 0 \
        else f"\n{category.__name__}  {message}"
    print(message)
    return message


warnings.showwarning = _showwarning
warnings.formatwarning = _showwarning


class ErrorAction(Enum):
    """ Simple Enum to control behavior when a problem is identified """
    Ignore = 0
    Warn = 1
    RaiseException = 2

    def handle_error(self, message):
        if self == ErrorAction.Warn:
            warnings.warn(message, stacklevel=1)
        elif self == ErrorAction.RaiseException:
            # remove last frame from traceback so that the exception points to the actual problem, and not this
            # wrapper
            # noinspection PyUnresolvedReferences,PyProtectedMember
            back_frame = sys._getframe(1)
            back_traceback = types.TracebackType(tb_next=None, tb_frame=back_frame, tb_lasti=back_frame.f_lasti,
                                                 tb_lineno=back_frame.f_lineno)
            raise ValueError(message).with_traceback(back_traceback)


class VcfKeys:  # note: for convenience we use .lower() before processing to make all fields lower-case
    id = "id"
    chrom = "chrom"
    pos = "pos"
    stop = "stop"
    qual = "qual"
    ref = "ref"
    alts = "alts"
    filter = "filter"
    rlen = "rlen"
    rid = "rid"
    svtype = "SVTYPE"
    svlen = "SVLEN"
    bnd_contig_2 = "CHR2"
    bnd_end_2 = "END2"
    source = "SOURCE"
    alleles = "alleles"
    gt = "GT"
    gq = "GQ"
    end = "END"  # Note: don't use this, use stop instead
    cn = "CN"  # Note: By the end of CleanVcf, CN/CNQ are the same as RD_CN/RD_GQ
    cnq = "CNQ"
    rd_cn = "RD_CN"
    rd_gq = "RD_GQ"


class Keys:
    id = "id"
    contig = "contig"
    begin = "begin"
    end = "end"
    gt = VcfKeys.gt.lower()
    gq = VcfKeys.gq.lower()
    allele_frequency = "af"
    qual = "qual"
    ref = "ref"
    alts = "alts"
    filter = "filter"
    alleles = "alleles"
    rlen = "rlen"
    rid = "rid"
    svtype = VcfKeys.svtype.lower()
    svlen = VcfKeys.svlen.lower()
    bnd_contig_2 = VcfKeys.bnd_contig_2.lower()
    bnd_end_2 = VcfKeys.bnd_end_2.lower()
    source = VcfKeys.source.lower()
    other_contig = "other_contig"
    other_begin = "other_begin"
    other_end = "other_end"
    sample_ids = "sample_ids"
    sample_id = "sample_id"
    property = "property"
    allele_count = "ac"
    cn = VcfKeys.cn.lower()
    cnq = VcfKeys.cnq.lower()
    rd_cn = VcfKeys.rd_cn.lower()
    rd_gq = VcfKeys.rd_gq.lower()


VcfMappingClasses = (pysam.libcbcf.VariantRecordSamples, pysam.libcbcf.VariantRecordSample)


class VcfFieldTypes:
    string = "String"
    integer = "Integer"
    float = "Float"


class BedKeys:
    chrom = "chrom"
    chrom_start = "chromStart"
    chrom_end = "chromEnd"
    other_chrom = "otherChrom"
    other_start = "otherStart"
    other_end = "otherEnd"
    name = "name"


class VaporKeys:
    gt = "gt"
    gq = "gq"
    qs = "qs"
    gs = "gs"
    read_scores = "read_scores"


class Default:
    int_type = numpy.int32
    float_type = numpy.float32
    missing_value = "MISSING"
    location_vcf_columns = MappingProxyType({
        VcfKeys.id: Keys.id, VcfKeys.chrom: Keys.contig, VcfKeys.pos: Keys.begin, VcfKeys.stop: Keys.end
    })
    basic_vcf_columns = MappingProxyType({
        **location_vcf_columns,
        **{VcfKeys.qual: Keys.qual, VcfKeys.ref: Keys.ref, VcfKeys.alts: Keys.alts, VcfKeys.filter: Keys.filter}
    })
    load_bed_columns = MappingProxyType({
        BedKeys.chrom: Keys.contig, BedKeys.chrom_start: Keys.begin, BedKeys.chrom_end: Keys.end, BedKeys.name: Keys.id,
        BedKeys.other_chrom: Keys.other_contig, BedKeys.other_start: Keys.other_begin,
        BedKeys.other_end: Keys.other_end
    })
    save_bed_columns = MappingProxyType(
        {val: key for key, val in load_bed_columns.items()}
    )
    genome_origin = 1
    bed_origin = 0
    vcf_origin = 1
    encoding = "utf-8"
    drop_missing_columns = False
    drop_trivial_multi_index = False
    cast_whole_nums_to_int = True
    category_prop_getter_ordered = False
    missing_samples_action = ErrorAction.RaiseException
    missing_properties_action = ErrorAction.RaiseException
    header_start = '#'  # most files use this as their header indicator
    log_progress = True
    location_columns = frozenset({Keys.begin, Keys.end, Keys.bnd_end_2, Keys.other_begin, Keys.other_end})
    use_copy_number = False
    use_cn = False  # Note: By the end of CleanVcf, CN/CNQ is identical to RD_CN/RD_GQ when it's present


def _number_more_than_1(vcf_number: str) -> bool:
    return vcf_number > 1 if isinstance(vcf_number, int) else True


def _get_genotype_extractor(vcf_prop_name: str, header: pysam.VariantHeader) -> GenotypePropertyExtractor:
    header_field = header.formats.get(vcf_prop_name)
    field_type, number = header_field.type, header_field.number
    if field_type == VcfFieldTypes.string or _number_more_than_1(number):
        encode_dict = CategoryPropertyCollator.get_encode_dict(vcf_prop_name)

        def extractor(genotype_record: pysam.libcbcf.VariantRecordSample) -> Optional[EncodedVcfField]:
            return CategoryPropertyCollator.encode(
                genotype_record.get(vcf_prop_name, None),
                encode_dict
            )
    else:
        def extractor(genotype_record: pysam.libcbcf.VariantRecordSample) -> Optional[EncodedVcfField]:
            return genotype_record.get(vcf_prop_name, None)
    return extractor


def _get_variant_extractor(
        vcf_prop_name: str, header: pysam.VariantHeader, missing_value: Any = Default.missing_value
) -> VariantPropertyExtractor:
    if vcf_prop_name in header.info:
        # get from info
        header_field = header.info.get(vcf_prop_name)
        field_type, number = header_field.type, header_field.number
        if field_type == VcfFieldTypes.string or _number_more_than_1(number):
            encode_dict = CategoryPropertyCollator.get_encode_dict(vcf_prop_name)

            def extractor(variant_record: pysam.VariantRecord) -> Optional[EncodedVcfField]:
                return CategoryPropertyCollator.encode(
                    variant_record.info.get(vcf_prop_name, None),
                    encode_dict
                )
        else:
            def extractor(variant_record: pysam.VariantRecord) -> Optional[EncodedVcfField]:
                return variant_record.info.get(vcf_prop_name, None)
    elif vcf_prop_name == VcfKeys.filter:
        # always a tuple of strings
        encode_dict = CategoryPropertyCollator.get_encode_dict(Keys.filter)

        def extractor(variant_record: pysam.VariantRecord) -> Optional[EncodedVcfField]:
            filters = tuple(variant_record.filter)
            return CategoryPropertyCollator.encode(filters if filters else missing_value, encode_dict)
    else:
        # get property directly from record. These are fields defined in the vcf spec
        if vcf_prop_name in {VcfKeys.pos, VcfKeys.stop, VcfKeys.rlen}:
            field_type = VcfFieldTypes.integer
        elif vcf_prop_name in {VcfKeys.qual}:
            field_type = VcfFieldTypes.float
        elif vcf_prop_name in {VcfKeys.id, VcfKeys.chrom, Keys.contig, VcfKeys.alts, VcfKeys.alleles, VcfKeys.ref}:
            # note: pysam allows querying chrom as contig, so check for both here in case user overrides default
            # properties
            field_type = VcfFieldTypes.string
        else:
            raise ValueError(f"Unknown vcf property '{vcf_prop_name}'")

        if field_type == VcfFieldTypes.string:
            encode_dict = CategoryPropertyCollator.get_encode_dict(vcf_prop_name)

            def extractor(variant_record: pysam.VariantRecord) -> Optional[EncodedVcfField]:
                return CategoryPropertyCollator.encode(
                    getattr(variant_record, vcf_prop_name, None),
                    encode_dict
                )
        else:
            def extractor(variant_record: pysam.VariantRecord) -> Optional[EncodedVcfField]:
                return getattr(variant_record, vcf_prop_name, None)
    return extractor


class VcfPropertyCollator:
    __slots__ = ("vcf_prop_name", "column_name")
    column_name: Tuple[Optional[Text], Text]

    def __init__(self, vcf_prop_name: str, column_name: Tuple[Optional[Text], Text]):
        self.vcf_prop_name = vcf_prop_name
        self.column_name = column_name

    @property
    def table_prop_name(self) -> str:
        return self.column_name[1]

    @property
    def sample_id(self) -> Optional[str]:
        return self.column_name[0]

    def get_dtype(self, values: Sequence[EncodedVcfField]) -> dtype:
        raise NotImplementedError("class derived from VcfPropertyCollator should implement .dtype()")

    def get_pandas_series(self, values: Sequence[EncodedVcfField]) -> pandas.Series:
        try:
            return pandas.Series(values, dtype=self.get_dtype(values), name=self.column_name)
        except Exception as err:
            common.add_exception_context(err, f"making series for {self.column_name}")
            raise

    @staticmethod
    def get_field_type_and_number(vcf_prop_name: str, header: pysam.VariantHeader) -> (str, int):
        if vcf_prop_name in header.info:
            # get from info
            header_field = header.info.get(vcf_prop_name)
            field_type, number = header_field.type, header_field.number
        elif vcf_prop_name in header.formats:
            header_field = header.formats.get(vcf_prop_name)
            field_type, number = header_field.type, header_field.number
        elif vcf_prop_name == VcfKeys.filter:
            field_type, number = VcfFieldTypes.string, -1
        else:
            # get property directly from record. These are fields defined in the vcf spec
            if vcf_prop_name in {VcfKeys.pos, VcfKeys.stop, VcfKeys.rlen}:
                field_type, number = VcfFieldTypes.integer, 1
            elif vcf_prop_name in {VcfKeys.qual}:
                field_type, number = VcfFieldTypes.float, 1
            elif vcf_prop_name in {VcfKeys.id, VcfKeys.chrom, Keys.contig, VcfKeys.alts, VcfKeys.alleles, VcfKeys.ref}:
                # note: pysam allows querying chrom as contig, so check for both here in case user overrides default
                # properties
                field_type, number = VcfFieldTypes.string, -1
            else:
                raise ValueError(f"Unknown vcf property '{vcf_prop_name}'")
        return field_type, number

    @staticmethod
    def get_collator(
            vcf_prop_name: str,
            column_name: Tuple[Optional[Text], Text],
            header: pysam.VariantHeader,
            missing_value: str = Default.missing_value,
            int_type: type = Default.int_type,
            float_type: type = Default.float_type,
            cast_whole_nums_to_int: bool = Default.cast_whole_nums_to_int,
            ordered: bool = Default.category_prop_getter_ordered
    ) -> "VcfPropertyCollator":
        field_type, number = VcfPropertyCollator.get_field_type_and_number(vcf_prop_name=vcf_prop_name, header=header)

        if field_type == VcfFieldTypes.integer:
            if number == 1:
                return IntPropertyCollator(vcf_prop_name=vcf_prop_name, column_name=column_name, int_type=int_type)
        if field_type == VcfFieldTypes.float:
            if number == 1:
                return FloatPropertyCollator(vcf_prop_name=vcf_prop_name, column_name=column_name,
                                             float_type=float_type, cast_whole_nums_to_int=cast_whole_nums_to_int)
        return CategoryPropertyCollator(vcf_prop_name=vcf_prop_name, column_name=column_name,
                                        missing_value=missing_value, ordered=ordered)


class CategoryPropertyCollator(VcfPropertyCollator):
    __slots__ = ("encode_dict", "missing_value", "ordered")
    encode_dict: Dict[Any, int]
    missing_value: Any
    ordered: bool

    def __init__(self, vcf_prop_name: str, column_name: Tuple[Optional[Text], Text],
                 missing_value: Any = Default.missing_value, ordered: bool = Default.category_prop_getter_ordered):
        super().__init__(vcf_prop_name=vcf_prop_name, column_name=column_name)
        # create global encode dict (to keep common encoding in case later another VCF is read in)
        self.encode_dict = CategoryPropertyCollator.get_encode_dict(self.vcf_prop_name)
        self.missing_value = missing_value
        self.ordered = ordered

    @staticmethod
    def get_encode_dict(encode_name: str) -> Dict:
        global_name = "_encode_" + encode_name
        if global_name in globals():
            encode_dict = globals()[global_name]
        else:
            encode_dict = {}
            globals()[global_name] = encode_dict
        return encode_dict

    @staticmethod
    def encode(raw_value: Any, code_dict: Dict[Any, int]) -> int:
        code = code_dict.get(raw_value, None)
        if code is None:
            code = len(code_dict)
            code_dict[raw_value] = code
        return code

    @property
    def categories(self) -> List[Any]:
        unique_values = [None] * len(self.encode_dict)
        for unique_value, code in self.encode_dict.items():
            unique_values[code] = self.missing_value if unique_value is None else unique_value
        return unique_values

    @staticmethod
    def _is_compressable(num_categories: int, num_values: int) -> bool:
        return num_categories < num_values // 2

    def get_dtype(self, values: Sequence[EncodedVcfField]) -> dtype:
        categories = self.categories
        if CategoryPropertyCollator._is_compressable(len(categories), len(values)):
            # this data is compressable, use categorical.
            return pandas.CategoricalDtype(categories=self.categories, ordered=self.ordered)
        else:
            return numpy.dtype("object")

    def get_pandas_series(self, values: Sequence[EncodedVcfField]) -> pandas.Series:
        try:
            categories = self.categories
            if CategoryPropertyCollator._is_compressable(len(categories), len(values)):
                # this data is compressable, use categorical
                try:
                    return pandas.Series(
                        pandas.Categorical.from_codes(values, categories=categories, ordered=self.ordered),
                        name=self.column_name
                    )
                except Exception as exception:
                    common.add_exception_context(
                        exception, f"Error compressing {len(values)} values into categories: {categories}"
                    )
                    raise
            else:
                # expected a categorical, but it's better to be just an object series
                return pandas.Series(
                    [categories[value] for value in values], name=self.column_name, dtype="object"
                )
        except Exception as err:
            common.add_exception_context(err,
                                         f"making series for sample {self.sample_id}, property {self.table_prop_name}")
            raise


class IntPropertyCollator(VcfPropertyCollator):
    __slots__ = ("int_type",)
    int_type: type

    def __init__(self, vcf_prop_name: str, column_name: Tuple[Optional[Text], Text],
                 int_type: type = Default.int_type):
        # ignore missing_value: for ints, we need to use NA
        super().__init__(vcf_prop_name=vcf_prop_name, column_name=column_name)
        self.int_type = int_type

    @staticmethod
    def is_integer_values(values: Union[Sequence[EncodedVcfField], numpy.ndarray, pandas.Series]) -> bool:
        return all(value is None or pandas.api.types.is_integer(value) for value in values)

    @staticmethod
    def make_nullable_dtype(
            dt: Union[numpy.dtype, pandas.api.extensions.ExtensionDtype]
    ) -> pandas.api.extensions.ExtensionDtype:
        """ convert to pandas dtype which supports pandas.NA (same class names with upper-case) """
        return pandas.core.dtypes.cast.pandas_dtype(
            dt.name.replace('i', 'I').replace('u', 'U')
        )

    @staticmethod
    def range_to_min_int_dtype(
            min_value: int, max_value: int, nullable: bool
    ) -> Union[numpy.dtype, pandas.api.extensions.ExtensionDtype]:
        # noinspection PyTypeChecker
        min_type: numpy.dtype = numpy.min_scalar_type(
            int(min(min_value, -max_value) if min_value < 0 else max_value)
        )

        return IntPropertyCollator.make_nullable_dtype(min_type) if nullable else min_type

    @staticmethod
    def min_int_dtype(
            values: Union[Sequence[EncodedVcfField], numpy.ndarray, pandas.Series],
            prop_name: str,
            int_type: type,
            location_columns: Collection[str] = Default.location_columns
    ) -> dtype:
        try:
            _has_missing = any(pandas.isna(value) for value in values) if isinstance(values, Tuple) \
                else pandas.isna(values).any()
        except AttributeError as attribute_error:
            common.add_exception_context(attribute_error, f"values={values}")
            raise

        if _has_missing:
            min_val = min((value for value in values if not pandas.isna(value)), default=0)
            max_val = max((value for value in values if not pandas.isna(value)), default=0)
        else:
            min_val = min(values, default=0)
            max_val = max(values, default=0)

        if _is_location_column(prop_name, location_columns=location_columns):
            # ensure that begin and end can handle addition / subtraction with typical-sized other values to avoid
            # potential overflows in downstream computations
            min_val = min(min_val, numpy.iinfo(int_type).min)
            max_val = max(max_val, numpy.iinfo(int_type).max)

        return IntPropertyCollator.range_to_min_int_dtype(min_val, max_val, nullable=_has_missing)

    def get_dtype(self, values: Sequence[EncodedVcfField]) -> dtype:
        return IntPropertyCollator.min_int_dtype(values, self.table_prop_name, self.int_type)

    def get_pandas_series(self, values: Sequence[EncodedVcfField]) -> pandas.Series:
        try:
            return pandas.Series(values, dtype=self.get_dtype(values), name=self.column_name)
        except TypeError:
            # Sometimes variants have multiple values for an "integer" property (e.g. svlen). Unfortunately, in *MANY*
            # places in the codebase, we must make decisions about svlen, etc as though they were simple integers. So if
            # this error occurs, try ONE TIME to coerce all values to simple integers and re-calculate min/max values
            values = [value if isinstance(value, (int, NoneType)) else max(value, key=abs) for value in values]
            try:
                return pandas.Series(values, dtype=self.get_dtype(values), name=self.column_name)
            except Exception as err:
                common.add_exception_context(
                    err, f"making series for sample {self.sample_id}, property {self.table_prop_name}"
                )
                raise


class FloatPropertyCollator(VcfPropertyCollator):
    __slots__ = ("float_type", "cast_whole_nums_to_int")
    float_type: type
    cast_whole_nums_to_int: bool

    def __init__(self, vcf_prop_name: str, column_name: Tuple[Optional[Text], Text],
                 float_type: type = Default.float_type, cast_whole_nums_to_int: bool = Default.cast_whole_nums_to_int):
        # ignore missing_value: for ints, we need to use NA
        super().__init__(vcf_prop_name=vcf_prop_name, column_name=column_name)
        self.float_type = float_type
        self.cast_whole_nums_to_int = cast_whole_nums_to_int

    @staticmethod
    def _float_values_can_be_cast_to_int(float_values: Iterable[float]) -> bool:
        return all(value is None or numpy.isnan(value) or value % 1 == 0 for value in float_values)

    @staticmethod
    def min_numeric_dtype(
            values: Sequence[EncodedVcfField],
            prop_name: str,
            cast_whole_nums_to_int: bool = Default.cast_whole_nums_to_int,
            int_type: type = Default.int_type,
            float_type: type = Default.float_type
    ) -> dtype:
        if IntPropertyCollator.is_integer_values(values) or \
                (cast_whole_nums_to_int and FloatPropertyCollator._float_values_can_be_cast_to_int(values)):
            return IntPropertyCollator.min_int_dtype(values, prop_name=prop_name, int_type=int_type)
        else:
            return numpy.dtype(float_type)

    def get_dtype(self, values: Sequence[EncodedVcfField], int_type: type = Default.int_type) -> dtype:
        return FloatPropertyCollator.min_numeric_dtype(
            values, self.table_prop_name, cast_whole_nums_to_int=self.cast_whole_nums_to_int,
            float_type=self.float_type, int_type=int_type
        )

    def get_pandas_series(self, values: Sequence[EncodedVcfField]) -> pandas.Series:
        try:
            if self.cast_whole_nums_to_int and FloatPropertyCollator._float_values_can_be_cast_to_int(values):
                return pandas.Series(
                    values,
                    dtype=IntPropertyCollator.min_int_dtype(values, prop_name=self.table_prop_name,
                                                            int_type=Default.int_type),
                    name=self.column_name
                )
            return pandas.Series(values, dtype=numpy.dtype(self.float_type), name=self.column_name)
        except Exception as err:
            common.add_exception_context(
                err, f"making series for sample {self.sample_id}, property {self.table_prop_name}"
            )
            raise


def _is_location_column(
        column_name: Optional[str],
        location_columns: Collection[str] = Default.location_columns
) -> bool:
    return column_name in location_columns


def _literal_eval(val: Any) -> Any:
    return val if isinstance(val, (float, int)) \
        else val if (isinstance(val, str) and not (val.startswith('(') and val.endswith(')'))) \
        else ast.literal_eval(val)


PandasType = TypeVar("PandasType", pandas.Series, pandas.DataFrame)


def make_nullable(pandas_obj: PandasType) -> PandasType:
    if isinstance(pandas_obj, pandas.Series):
        return pandas_obj.astype(IntPropertyCollator.make_nullable_dtype(pandas_obj.dtype))
    else:
        return pandas_obj.astype(
            {
                col_name: IntPropertyCollator.make_nullable_dtype(col_dtype)
                for col_name, col_dtype in pandas_obj.dtypes.items()
            }
        )


def compress_types(
        pandas_obj: Union[pandas.DataFrame, pandas.Series],
        int_type: type = Default.int_type,
        float_type: type = Default.float_type,
        missing_value: Any = Default.missing_value,
        cast_whole_nums_to_int: bool = Default.cast_whole_nums_to_int,
        literal_eval_columns: Optional[Set[str]] = None
) -> Union[pandas.DataFrame, pandas.Series]:
    f"""
    Decrease memory usage of pandas object by
        -reducing int and float types
        -compressing objects to categoricals. Note: NaN pandas can't use NaN as a category (NaNs can be in categorical
        series, but they don't show up in .cat.categories), so object Series have NaNs changed to missing_value (unless
        missing_value is set to NaN)
    Args:
        pandas_obj: Union[pandas.DataFrame, pandas.Series]
            pandas object with data to compress
        int_type: numpy field_type (Default = {Default.int_type})
            Minimum field_type for "begin" and "end". Other int columns are compressed into smallest viable integer
            field_type.
        float_type: numpy field_type (Default = {Default.float_type})
            Type for any float columns
        missing_value: str (Default = {Default.missing_value})
            Value to replace NaNs with in object columns
        cast_whole_nums_to_int: bool (Default={Default.cast_whole_nums_to_int})
            If True and the values are float but are all whole numbers, covert to ints
        literal_eval_columns: Optional[Set[str]] (default=None)
            Columns that need to be "evaluated" to convert from string to some more complex type. A typical example
            would be "alts", which is frequently a tuple. If None, then attempt to convert every object column.
    Returns:
        compressed_pandas_obj: Union[pandas.DataFrame, pandas.Series]
            pandas_obj with values compressed
    """
    if isinstance(pandas_obj, pandas.DataFrame):
        for col in pandas_obj.columns:
            pandas_obj[col] = compress_types(pandas_obj[col], int_type=int_type, float_type=float_type,
                                             missing_value=missing_value, literal_eval_columns=literal_eval_columns)
        return pandas_obj

    if pandas_obj.dtype == object:
        if literal_eval_columns is None or pandas_obj.name in literal_eval_columns:
            try:
                pandas_obj = pandas_obj.map(_literal_eval, na_action="ignore")
            except ValueError as value_error:
                common.add_exception_context(value_error, f"Error performing literal_eval on {pandas_obj.name}")
                raise
        return pandas_obj.fillna(missing_value).astype("category") if pandas_obj.nunique() < len(pandas_obj) // 2 \
            else pandas_obj.fillna(missing_value)
    else:
        numeric_dtype = FloatPropertyCollator.min_numeric_dtype(
            pandas_obj, prop_name=str(pandas_obj.name), cast_whole_nums_to_int=cast_whole_nums_to_int,
            int_type=int_type, float_type=float_type
        )
        # don't up-cast to larger type than is originally used
        return pandas_obj.astype(numeric_dtype) if numeric_dtype.type(0).nbytes <= pandas_obj.dtype.type(0).nbytes \
            else pandas_obj


def shift_origin(
        df: pandas.DataFrame, current_origin: int, desired_origin: int, copy_on_change: bool = False
) -> pandas.DataFrame:
    """
    Utility function to alter intervals dataframe between BED and VCF formats
    """
    # note, only shift "begin" because VCF and BED format have different origin but different open/closed structure, so
    # "end" doesn't change
    if current_origin == desired_origin or Keys.begin not in df.columns:
        return df
    if copy_on_change:
        df = df.copy()
    shift = desired_origin - current_origin
    df[Keys.begin] += shift
    return df


def _get_remapped_columns(
        df_columns: Union[pandas.Index, Sequence[str]],
        columns: Union[Sequence[str], Mapping[str, str]],
        first_columns: Sequence[str] = ()
) -> Dict[str, str]:
    """
    Helper function for _remap_columns, converts Sequence columns to a Mapping, and handles ordering issues.
    Args:
        df_columns: Union[pandas.Index, Sequence[str]]
            Columns of original, unaltered DataFrame
        columns: Union[Sequence[str], Mapping[str, str], None]
            if None:       return unchanged df
            if a Sequence: alter df to include only the requested columns in the requested order
            if a Mapping:  alter df to rename columns (from keys to values of mapping) and re-order them to the same
                           order as the mapping
        first_columns: Sequence[str] (Default=())
            For situations where some columns must be present, and go first in a specific order (typically BED files)
            but remapping is also necessary, the leading columns can be specified here.
    Returns:
        remapped_columns: Dict[str, str]
            Mapping from initial column names to output names, in the appropriate order for the final output DataFrame
    """
    if isinstance(columns, Mapping):
        invert_map = {val: key for key, val in columns.items()}
        # note this will throw an exception if anything is missing from first_columns
        remapped_columns = {
            **{invert_map[first_col]: first_col for first_col in first_columns},
            **{key: value for key, value in columns.items() if value not in first_columns and key in df_columns},
            **{col: col for col in df_columns if col not in columns}
        }
    elif columns is not None:
        remapped_columns = {
            **{col: col for col in first_columns},
            **{col: col for col in columns if col not in first_columns}
        }
    else:
        remapped_columns = {
            **{col: col for col in first_columns},
            **{col: col for col in df_columns if col not in first_columns}
        }
    return remapped_columns


def _remap_columns(
        df: pandas.DataFrame,
        columns: Union[Sequence[str], Mapping[str, str], None],
        copy_on_change=False,
        first_columns: Sequence[str] = ()
) -> pandas.DataFrame:
    """
    Rename and re-order columns of pandas Dataframe to match requested names and orders
    Args:
        df: pandas.DataFrame
            input DataFrame to alter
        columns: Union[Sequence[str], Mapping[str, str], None]
            if None:       return unchanged df
            if a Sequence: alter df to include only the requested columns in the requested order
            if a Mapping:  alter df to rename columns (from keys to values of mapping) and re-order them to the same
                           order as the mapping
        copy_on_change: bool (default=False)
            If no changes are necessary, always return original df
            if True, return a new DataFrame copy if changes are necessary
            if False, make any changes in place
        first_columns: Sequence[str] (Default=())
            For situations where some columns must be present, and go first in a specific order (typically BED files)
            but remapping is also necessary, the leading columns can be specified here.
    Returns:
        remapped_df: pandas.DataFrame
            DataFrame with columns renamed and re-ordered
    """
    if columns is not None:
        remapped_columns = _get_remapped_columns(df.columns, columns, first_columns)
        if list(remapped_columns.keys()) != list(remapped_columns.values()) or \
                list(remapped_columns.values()) != df.columns.to_list():
            # need to change dataframe columns
            if copy_on_change:
                df = df.loc[:, remapped_columns.keys()].rename(columns=remapped_columns, copy=False)
            else:
                df.drop(columns=list(set(df.columns).difference(remapped_columns.keys())), inplace=True)
                df.rename(columns=remapped_columns, inplace=True)
    return df


def _tuple_to_str(val: Any) -> Any:
    return str(val).replace(' ', '') if isinstance(val, tuple) else val


def pandas_to_tsv(
        data_file: Text,
        df: pandas.DataFrame,
        columns: Union[Sequence[str], Mapping[str, str], None] = None,
        first_columns: Sequence[str] = (),
        write_index: bool = False,
        write_header: bool = True,
        header_start: str = Default.header_start,
        genome_origin: int = Default.genome_origin,
        tsv_origin: int = Default.genome_origin,
        encoding: str = Default.encoding
):
    f"""
    Save pandas DataFrame into tab-delimited-gzipped file. Notable extra features:
        -Rename columns from internal schema to safe-file schema
        -Shift genomic coordinates from internal origin to save-file origin if different
    Args:
        data_file: Text
            Full path to save file
        df: pandas.DataFrame
            Table of data.
        columns: Mapping[str, str], Sequence[str], or None (Default = None)
            If a Mapping: when saving to disk, rename any column names that are keys in this mapping.
            If a Sequence: when saving to disk, set the column names to this sequence
        first_columns: Sequence[str] (default=())
            For situations where some columns must be present, and go first in a specific order (typically BED files)
            but remapping is also necessary, the leading columns can be specified here.
        write_index: bool (Default = False)
            If true, write the DataFrame index as the first output column, if False, omit the index
        write_header: bool (Default = True)
            If true, begin save file with header of column names
        header_start: str (Default={Default.header_start})
            Start header with this string.
        genome_origin: int (Default = {Default.genome_origin})
            Indexing origin for DataFrame.
        tsv_origin: int (Default = {Default.genome_origin})
            Indexing origin for output bed file. Shift output values if needed.
        encoding: str (Default = {Default.encoding}
            Encoding to use when writing strings.
    """
    if isinstance(df.columns, pandas.MultiIndex):
        raise ValueError("Unable to output multi-index columns as TSV. Manually flatten the column labels first.")
    df = _remap_columns(
        shift_origin(df, current_origin=genome_origin, desired_origin=tsv_origin, copy_on_change=True),
        columns=columns, first_columns=first_columns, copy_on_change=True
    ).applymap(_tuple_to_str)
    with pysam.BGZFile(data_file, "wb") as f_out:
        if write_header and header_start:
            f_out.write(header_start.encode(encoding))
        f_out.write(df.to_csv(sep='\t', index=write_index, header=write_header, na_rep=numpy.nan).encode(encoding))


def pandas_to_bed(
        data_file: Text,
        df: pandas.DataFrame,
        build_tabix_index: bool = True,
        columns: Union[Sequence[str], Mapping[str, str], None] = Default.save_bed_columns,
        write_header: bool = True,
        header_start: str = Default.header_start,
        genome_origin: int = Default.genome_origin,
        bed_origin: int = Default.bed_origin,
        encoding: str = Default.encoding
):
    f"""
    Save pandas DataFrame into .bed file. Notable extra features:
        -Rename columns from internal schema to safe-file schema
        -Shift genomic coordinates from internal origin to save-file origin if different
        -Build tabix index if requested.
    This is an ease-of-use function that just calls pandas_to_tsv() with appropriate defaults.
    Args:
        data_file: Text
            Full path to save file
        df: pandas.DataFrame
            Table of data.
        build_tabix_index: bool (Default = True)
            If true, create tabix index file.
        columns: Mapping[str, str], Sequence[str], or None (Default = {Default.save_bed_columns})
            If a Mapping: when saving to disk, rename any column names that are keys in this mapping.
            If a Sequence: when saving to disk, set the column names to this sequence
        write_header: bool (Default = True)
            If true, begin save file with header of column names
        header_start: str (Default={Default.header_start})
            Start header with this string.
        genome_origin: int (Default = {Default.genome_origin})
            Indexing origin for DataFrame.
        bed_origin: int (Default = {Default.bed_origin})
            Indexing origin for output bed file. Shift output values if needed.
        encoding: str (Default = {Default.encoding}
            Encoding to use when writing strings.
    """
    if isinstance(df.columns, pandas.MultiIndex):
        raise ValueError("Unable to output multi-index columns as BED. Manually flatten the column labels first.")
    # turn index into regular column so that mandatory bed ordering for first few columns can be adhered to
    df = df.reset_index(col_level=Keys.property)
    # now write out TSV file, with appropriate options
    pandas_to_tsv(
        data_file=data_file, df=df, columns=columns,
        first_columns=(BedKeys.chrom, BedKeys.chrom_start, BedKeys.chrom_end, BedKeys.name),
        write_header=write_header, header_start=header_start, write_index=False,
        genome_origin=genome_origin, tsv_origin=bed_origin, encoding=encoding
    )
    # build index if requested
    if build_tabix_index:
        pysam.tabix_index(data_file, preset="bed", force=True)


def tsv_to_pandas(
        data_file: Text,
        columns: Union[Sequence[Text], Mapping[Text, Text], None] = Default.load_bed_columns,
        int_type: type = Default.int_type,
        missing_value: str = Default.missing_value,
        float_type: type = Default.float_type,
        genome_origin: int = Default.genome_origin,
        tsv_origin: int = Default.genome_origin,
        sort_intervals: bool = False,
        header_start: str = '#',
        require_header: bool = True,
        encoding: str = Default.encoding,
        log_progress: bool = Default.log_progress,
        literal_eval_columns: Set[str] = frozenset(),
        *args, **kwargs
) -> pandas.DataFrame:
    f"""
    Load dataframe from generic tab-delimited-gzipped file. If specifically loading from a .bed file, use bed_to_pandas
    which calls this function with appropriate arguments.
       - Treat lines starting with {header_start} as comments and ignore.
       - Rename columns from file schema to desired schema
       - Shift genomic coordinates from file origin ({tsv_origin}) to desired origin ({genome_origin})
       - Attempt to compress data optimally
           - Use {float_type} for float columns
           - Use {int_type} for genomic coordinate columns
           - Use smallest possible int for all other int columns
           - Use categoricals for object columns provided the number of categories is < 1/2 the number of rows
    Args:
        data_file: Text
            Full path to file
        columns: Mapping[str, str], Sequence[str], or None (Default = {Default.load_bed_columns})
            If a Mapping: when loading from disk, rename any column names that are keys in this mapping.
            If a Sequence: when loading from disk, set the column names to this sequence
        int_type: numpy field_type (Default = {Default.int_type})
            Minimum field_type for "begin" and "end". Other int columns are compressed into smallest viable integer
            field_type.
        float_type: numpy field_type (Default = {Default.float_type})
            Type for any float columns
        missing_value: str (Default = {Default.missing_value})
            Value to replace NaNs with in object columns
        genome_origin: int (Default = {Default.genome_origin})
            Desired indexing origin for DataFrame. bed files use 0-indexing, so shift values in DataFrame if needed.
        tsv_origin: int (Default = {Default.genome_origin})
            Indexing origin for bed  file.
        sort_intervals: bool (Default=False)
            If set to True, sort the intervals into canonical order
        header_start: str (Default='#')
            Assume header starts with this string.
            If set to empty string, assume first line is header (if require_header is True)
        require_header: bool (Default=True)
            If set to True, throw error if a header (starting with header_start) is not present in the file.
            If False, use columns (if provided) to create a header) or else just integers
        encoding: str (Default={Default.encoding})
            Encoding to use for file
        log_progress: bool (Default={Default.log_progress})
            Display file name and requested properties on loading file
        literal_eval_columns: Iterable[str] (default=())
            Columns that need to be "evaluated" to convert from string to some more complex type. A typical example
            would be "alts", which is frequently a tuple.
        *args, **kwargs: passed to pandas.read_csv
    Returns:
        df: pandas.DataFrame
            Table of data.
    """
    if log_progress:
        if columns is None:
            print(f"Loading {data_file} ...", flush=True, file=sys.stderr, end="")
        else:
            target_columns = columns.values() if isinstance(columns, Mapping) else columns
            print(f"Reading {','.join(target_columns)} from {data_file} ...", flush=True, file=sys.stderr, end="")

    if not os.path.isfile(data_file):
        raise ValueError(f"{data_file} does not exist")
    # it is *VASTLY* faster to handle header manually, then load remaining file into buffer and call pandas.read_csv
    # on the buffer than it is to load line-by-line with python code.
    header_start = header_start.encode(encoding)
    buffer = io.StringIO()
    with pysam.BGZFile(data_file, "rb") as f_in:
        file_start = f_in.read(len(header_start)) if header_start else "".encode(encoding)
        if require_header:
            if header_start and file_start != header_start:
                # require a header but one wasn't found
                raise ValueError(f".bed ({data_file}) file does not start with {header_start.decode(encoding)} header")
        else:
            # don't require a header in the file, create one
            if columns is None:
                first_line = (file_start + f_in.readline()).decode(encoding) + '\n'
                # make header by counting number of columns in first line, and numbering columns
                header = '\t'.join(str(x) for x in range(len(first_line.split('\t')))) + '\n'
                buffer.write(header)  # insert header onto front of file
                buffer.write(first_line)  # put back first line, after the new header
            else:
                # make header from columns
                header = '\t'.join(columns if isinstance(columns, Sequence) else columns.keys()) + '\n'
                buffer.write(header)  # insert header onto front of file
                buffer.write(file_start.decode(encoding))  # put back file_start, after the new header

        # now that header stuff is taken care of, read rest of file into buffer
        buffer.write(f_in.read().decode(encoding))

    buffer.seek(0)
    # the specified optional arguments to read_csv dramatically increase speed and stability
    # noinspection PyTypeChecker
    df = pandas.read_csv(
        buffer, sep='\t', low_memory=True, engine='c', memory_map=False, *args, **kwargs
    )

    df = shift_origin(
        _remap_columns(df, columns=columns), current_origin=tsv_origin, desired_origin=genome_origin
    )
    df = compress_types(df, int_type=int_type, float_type=float_type, missing_value=missing_value,
                        literal_eval_columns=literal_eval_columns)
    if Keys.id in df.columns:
        df.set_index(Keys.id, inplace=True)

    if sort_intervals:
        sort_intervals_table(df, inplace=True)
    if log_progress:
        print(" done", flush=True, file=sys.stderr, end="")
    return df


def bed_to_pandas(
        data_file: Text,
        columns: Union[Sequence[Text], Mapping[Text, Text], None] = Default.load_bed_columns,
        int_type: type = Default.int_type,
        float_type: type = Default.float_type,
        missing_value: str = Default.missing_value,
        genome_origin: int = Default.genome_origin,
        bed_origin: int = Default.bed_origin,
        sort_intervals: bool = False,
        header_start: str = '#',
        require_header: bool = True,
        encoding: str = Default.encoding,
        log_progress: bool = Default.log_progress,
        literal_eval_columns: Set[str] = frozenset(),
        *args, **kwargs
) -> pandas.DataFrame:
    f"""
    Load dataframe from tab-delimited-gzipped file (usually a .bed file).
       - Treat lines starting with {header_start} as comments and ignore.
       - Rename columns from file schema to desired schema
       - Shift genomic coordinates from file origin ({bed_origin}) to desired origin ({genome_origin})
       - Attempt to compress data optimally
           - Use {float_type} for float columns
           - Use {int_type} for genomic coordinate columns
           - Use smallest possible int for all other int columns
           - Use categoricals for object columns provided the number of categories is < 1/2 the number of rows
    This is just an ease-of-use function that calls tsv_to_pandas with appropriate defaults.
    Args:
        data_file: Text
            Full path to file
        columns: Mapping[str, str], Sequence[str], or None (Default = {Default.load_bed_columns})
            If a Mapping: when loading from disk, rename any column names that are keys in this mapping.
            If a Sequence: when loading from disk, set the column names to this sequence
        int_type: numpy field_type (Default = {Default.int_type})
            Minimum field_type for "begin" and "end". Other int columns are compressed into smallest viable integer
            field_type.
        float_type: numpy field_type (Default = {Default.float_type})
            Type for any float columns
        missing_value: str (Default = {Default.missing_value})
            Value to replace NaNs with in object columns
        genome_origin: int (Default = {Default.genome_origin})
            Desired indexing origin for DataFrame. bed files use 0-indexing, so shift values in DataFrame if needed.
        bed_origin: int (Default = {Default.bed_origin})
            Indexing origin for bed  file.
        sort_intervals: bool (Default=False)
            If set to True, sort the intervals into canonical order
        header_start: str (Default='#')
            Assume header starts with this string.
            If set to empty string, assume first line is header (if require_header is True)
        require_header: bool (Default=True)
            If set to True, throw error if a header (starting with header_start) is not present in the file.
            If False, use columns (if provided) to create a header) or else just integers
        encoding: str (Default={Default.encoding})
            Encoding to use for file
        log_progress: bool (Default={Default.log_progress})
            Display file name and requested properties on loading file
        literal_eval_columns: Iterable[str] (default=())
            Columns that need to be "evaluated" to convert from string to some more complex type. A typical example
            would be "alts", which is frequently a tuple.
        *args, **kwargs: passed to pandas.read_csv
    Returns:
        df: pandas.DataFrame
            Table of data.
    """

    return tsv_to_pandas(
        data_file=data_file, columns=columns, int_type=int, float_type=float_type, missing_value=missing_value,
        genome_origin=genome_origin, tsv_origin=bed_origin, sort_intervals=sort_intervals, header_start=header_start,
        require_header=require_header, encoding=encoding, log_progress=log_progress,
        literal_eval_columns=literal_eval_columns, *args, **kwargs
    )


@common.static_vars(is_autosome_map={})
def is_autosome(contig: Text) -> bool:
    _is_autosome = is_autosome.is_autosome_map.get(contig, None)
    if _is_autosome is None:
        if contig.startswith("chr"):
            try:
                contig_num = int(contig[3:])
                _is_autosome = 1 <= contig_num <= 22
            except ValueError:  # pragma: no cover
                _is_autosome = False
        else:
            _is_autosome = False
        is_autosome.is_autosome_map[contig] = _is_autosome
    return _is_autosome


def get_is_autosome(variants: pandas.DataFrame) -> pandas.Series:
    """ Get boolean Series specifying if each variant is in the autosome (standard diploid/non-sex chromosomes) """
    return variants[
        (None, Keys.contig) if variants.columns.nlevels == 2 else Keys.contig
    ].map(is_autosome).astype(bool)


def contig_sort_key(contig: Text) -> str:
    """
    Used to put contigs (including alt contigs) into a sensible sorted order.
    Args:
        contig: str
            Name of contig
    Returns:
        sort_key: str
            Key for sorting contig
    """
    if contig == "chrX" or contig == "chrY":
        return contig
    try:
        contig_num = int(contig[3:])
        return contig if contig_num >= 10 else "chr0" + contig[-1]
    except ValueError:  # pragma: no cover
        return "zzz" + contig


def sort_intervals_table(
        intervals_df: pandas.DataFrame,
        inplace: bool = False,
        drop_index: bool = True,
) -> pandas.DataFrame:
    """
    Sort table of genomic intervals so that
    a) contigs are in cannonical order (primary first, then lexical) and all the
       intervals in a contig are sequential in the table.
    b) within a contig, intervals have non-decreasing "begin"
    c) within a contig, intervals with the same "begin" have non-decreasing "end"
    Ensure that final output DataFrame has "contig" as a categorical column.
    Args:
        intervals_df: pandas.DataFrame
            table of genomic intervals (with columns "contig", "begin", "end")
        inplace: bool (Default=False)
            if True, update intervals_df in place
            if False, return a new intervals_df table
        drop_index: bool (Default=True)
            if True, replace old index with an index from 0 to num_intervals - 1
            if False, keep old index (which will now be out of order)
    Returns:
        intervals: pandas.DataFrame
            sorted table of genomic intervals
    """
    if pandas.api.types.is_categorical_dtype(intervals_df[Keys.contig]):
        # if contig column is already categorical, ensure that the categories are in sorted order
        contigs = intervals_df[Keys.contig].cat.categories
        sorted_contigs = sorted(contigs, key=contig_sort_key)
        if not (intervals_df[Keys.contig].cat.ordered and contigs.equals(sorted_contigs)):
            # they're not ordered and in order, so swap them around and set to ordered
            intervals_df[Keys.contig] = intervals_df[Keys.contig].cat.reorder_categories(sorted_contigs, ordered=True)
    else:
        # get a list of all the contigs, sorted in correct order (e.g. chr11 comes *after* chr3)
        contigs = sorted(intervals_df[Keys.contig].unique(), key=contig_sort_key)
        # make contig column categorical, with category codes in desired order
        intervals_df[Keys.contig] = pandas.Categorical(intervals_df[Keys.contig].values, ordered=True,
                                                       categories=contigs)

    # now sorting by contig will actually sort by category code, which is in the desired order
    # sort intervals_df by contig, then within contig, by begin and end
    if inplace:
        intervals_df.sort_values([Keys.contig, Keys.begin, Keys.end], inplace=True)
    else:
        intervals_df = intervals_df.sort_values([Keys.contig, Keys.begin, Keys.end])
    if drop_index:
        intervals_df.reset_index(drop=True, inplace=True)
    return intervals_df


def vcat_with_categoricals(
        dataframes: Sequence[pandas.DataFrame],
        categorical_selection: str = "union",
        **kwargs
) -> pandas.DataFrame:
    """
    Concatenate DataFrames vertically while preserving categorical columns.
    NOTE: this function alters the categories in-place for input DataFrames
    Args:
        dataframes: Sequence[pandas.DataFrame]
            Variable number of interval tables to concatenate into one, column-wise
        categorical_selection: str (default = "union")
            If "union" make column categorical if any Dataframe with that column is categorical
            If "intersecion" make column categorical if all Dataframes with that column are categorical
            Otherwise throw ValueError
        **kwargs:
            keyword args to pass to pandas.concat. Don't attempt to override axis, this is only for vertical concat.
    Returns:
        vcat_table: pandas.DataFrame
            concatenated
    """
    if len(dataframes) <= 1:
        if len(dataframes) == 0:
            raise ValueError("Empty list of DataFrames supplied")
        return dataframes[0]
    if categorical_selection == "union":
        select_func = set.union
    elif categorical_selection == "intersection":
        select_func = set.intersection
    else:
        raise ValueError(f"categorical_selection must be 'union' or 'intersection', got '{categorical_selection}'")
    # Iterate on columns that are categorical in any dataframe
    for col in select_func(
        *[
            set(dataframe.select_dtypes(include="category").columns)
            for dataframe in dataframes
        ]
    ):
        dtypes = None
        try:
            dtypes = [dataframe[col].dtype for dataframe in dataframes if col in dataframe]
            if all(pandas.api.types.is_categorical_dtype(_dtype) for _dtype in dtypes) and \
                    all(_dtype == dtypes[0] for _dtype in dtypes):
                continue  # don't need to mess with this column, all dtypes are consistent
        except TypeError as err:
            context = f"{len(dataframes)} dataframes, {len(dtypes)} with col={col}, type(col)={type(col)}\n"
            if dtypes is not None:
                context += '\n'.join(f"\tdtype={str(_dtype)}" for _dtype in dtypes)

            common.add_exception_context(err, context)
            raise

        # Generate the union category across dfs for this column. convert to category if dataframe doesn't already have
        # column as categorical
        union_categories = set.union(*(
            set(dataframe[col].cat.categories if pandas.api.types.is_categorical_dtype(dataframe[col])
                else dataframe[col].unique())
            for dataframe in dataframes
        ))
        # Change to union category for all dataframes. We don't know anything about the ordering of categories, so
        # set ordered=False
        for dataframe in dataframes:
            dataframe[col] = pandas.Categorical(dataframe[col].values, categories=union_categories, ordered=False)
    return pandas.concat(dataframes, axis=0, **kwargs)


def vapor_to_pandas(vapor_file: Text) -> pandas.DataFrame:
    return bed_to_pandas(
        vapor_file,
        columns=(Keys.contig, Keys.begin, Keys.end, Keys.svtype, Keys.id, VaporKeys.qs, VaporKeys.gs, VaporKeys.gt,
                 VaporKeys.gq, VaporKeys.read_scores),
        require_header=False
    )


def _get_all_properties(
        header: pysam.VariantHeader,
        basic_vcf_columns: Mapping[str, str] = Default.basic_vcf_columns
) -> Dict[str, str]:
    f"""
    Get all the unique and useful properties in the VCF:
        -basic properties like CHROM, POS, END, FILTER, etc
        -INFO fields
        -per-GT fields specified in FORMATS
    Args:
        header: pysam.VariantHeader
            Header object for VCF
        basic_vcf_columns: Mapping[str, str] (Default={Default.basic_vcf_columns})
            Mapping of the basic fields, common to all VCFs
    Returns:
        all_properties: Dict[str, str]
            Mapping from desired VCF property names to the corresponding pandas DataFrame names
    """
    # note .info["END"] is reserved and accessed by .stop, don't get it
    table_columns = set(basic_vcf_columns.values())
    all_properties = {
        **basic_vcf_columns,
        **{k: k.lower() for k in header.info.keys()
           if k != VcfKeys.end and k not in basic_vcf_columns and k.lower() not in table_columns},
        **{k: k.lower() for k in header.formats.keys()}
    }
    return all_properties


def _get_allowed_properties(
        header: pysam.VariantHeader,
        basic_vcf_columns: Mapping[str, str] = Default.basic_vcf_columns
) -> Dict[str, str]:
    f"""
    Get all the unique and useful properties in the VCF, and add in a few properties that are less useful because they
    are duplicates, or not particularly useful, hence all allowed properties:
        -basic properties like CHROM, POS, END, FILTER, etc
        -INFO fields
        -per-GT fields specified in FORMATS
        -extra basic fields: rlen, rid, alleles
    Args:
        header: pysam.VariantHeader
            Header object for VCF
        basic_vcf_columns: Mapping[str, str] (Default={Default.basic_vcf_columns})
            Mapping of the basic fields, common to all VCFs
    Returns:
        allowed_properties: Dict[str, str]
            Mapping from desired VCF property names to the corresponding pandas DataFrame names
    """
    # return duplicates or weird things
    return {**_get_all_properties(header, basic_vcf_columns=basic_vcf_columns),
            **{VcfKeys.rlen: Keys.rlen, VcfKeys.rid: Keys.rid, VcfKeys.alleles: Keys.alleles}}


def get_vcf_properties(vcf: str) -> ValuesView[str]:
    with pysam.VariantFile(vcf, 'r', threads=common.num_physical_cpus) as f_in:
        return _get_allowed_properties(f_in.header).values()


def get_vcf_sample_ids(vcf: str) -> List[str]:
    with pysam.VariantFile(vcf, 'r') as f_in:
        return list(f_in.header.samples)


def get_vcf_variant_ids(vcf: str) -> pandas.Index:
    return vcf_to_pandas(vcf, wanted_properties=(Keys.id,)).index


def drop_trivial_columns_multi_index(variants: pandas.DataFrame) -> pandas.DataFrame:
    if len(variants.columns) > 0:
        for level in (Keys.sample_id, Keys.property):
            level_values = variants.columns.get_level_values(level)
            # noinspection PyUnresolvedReferences
            if level_values.isnull().all() or (level_values == level_values[0]).all():
                variants.columns = variants.columns.droplevel(level)  # drop this multi-index
    return variants


def vcf_to_pandas(
        vcf: str,
        samples: Optional[Collection[str]] = None,
        wanted_properties: Optional[Collection[str]] = None,
        genome_origin: int = Default.genome_origin,
        vcf_origin: int = Default.vcf_origin,
        int_type: type = Default.int_type,
        float_type: type = Default.float_type,
        missing_value: str = Default.missing_value,
        drop_missing_columns: bool = Default.drop_missing_columns,
        drop_trivial_multi_index: bool = Default.drop_trivial_multi_index,
        cast_whole_nums_to_int: bool = Default.cast_whole_nums_to_int,
        missing_samples_action: ErrorAction = Default.missing_samples_action,
        missing_properties_action: ErrorAction = Default.missing_properties_action,
        category_prop_getter_ordered: bool = Default.category_prop_getter_ordered,
        log_progress: bool = Default.log_progress,
        num_threads: int = common.num_logical_cpus
) -> pandas.DataFrame:
    f"""
    Load data from VCF into a pandas DataFrame.
      -Because there may be multiple samples, each with multiple properties, the DataFrame columns are multi-index by
       default. The first level indicates the sample the property refers to (None for variant properties like
       {Keys.begin}), and the second level gives the property name.
      -Property (DataFrame column) names will be lowercase. Standard properties will be renamed to the elements of the
       static Keys class, providing consistency between VCFs, BEDs, and other formats for genomic loci.
      -Column values are compressed like so:
        -position columns are set to {int_type}, other integers are compressed to smallest int type that will hold their
         value
        -float columns are set to {float_type}
        -string values (e.g. SVTYPE) are stored as CategoricalDTypes
    Args:
        vcf: str
            Full path to VCF
        samples: Optional[Sequence[str]] (default=None)
            If None, get data for all samples in VCF. Otherwise, get only the samples provided. If extra sample IDs are
            requested that are not present in the VCF, behavior is controlled by the value of missing_samples_action.
        wanted_properties: Optional[Collection[str] (default=None)
            If None, get all available properties in the VCF (determined by call to _get_all_properties()).
            Otherwise get the requested properties (which should be lower-case, and standard properties are elements of
            Keys). If extra property names are requested that are not present in the VCF, behavior is controlled by the
            value of missing_properties_action
        genome_origin: int (default={Default.genome_origin})
            Desired origin of genome (e.g. 0 or 1) for the table. If needed, {Keys.begin} will be shifted to match
            the desired origin.
        vcf_origin: int (default={Default.vcf_origin})
            Origin of genome for data in the VCF file.
        int_type: type (default={Default.int_type})
            Type for position values ({Keys.begin}, {Keys.end})
        float_type: (default={Default.float_type})
            Type for any floating values.
        missing_value: str (default={Default.missing_value})
            Value to use for missing non-numerical (e.g. string/categorical) properties.
        drop_missing_columns: bool (default={Default.drop_missing_columns})
            If True and an entire column is missing data, drop the column. Otherwise keep it.
        drop_trivial_multi_index: bool (default={Default.drop_trivial_multi_index})
            If True and one of the levels of the multi-index columns is trivial (e.g. there is no format data)
            the drop the trivial level. Otherwise keep it.
        cast_whole_nums_to_int: bool (default={Default.cast_whole_nums_to_int})
            If True then check if all the elements of a float column are whole numbers, and if so, cast to the minimal
            int that can store them (to save memory). If False, don't check and don't convert.
        missing_samples_action: ErrorAction (default={Default.missing_samples_action})
            If there are requested samples that are not present in the VCF, then based on this value:
                -If Ignore, then do nothing
                -If Warn, then raise a warning showing which samples are missing
                -If RaiseException, then raise a ValueError
        missing_properties_action: bool (default={Default.missing_properties_action})
            If properties are requested that are not present in the VCF, then based on this value:
                -If Ignore, then do nothing
                -If Warn, then raise a warning showing which properties are missing
                -If RaiseException, then raise a ValueError
        category_prop_getter_ordered: bool (default={Default.category_prop_getter_ordered})
            If True, create categorical properties as ordered CategoricalDTypes. Don't mess with this unless you're 100%
            sure you know what you're doing.
        log_progress: bool (default={Default.log_progress})
            If True, display a few lines to indicate how loading the VCF is proceeding.
        num_threads: int (default={common.num_logical_cpus})
            Use this many threads for decompression while reading the VCF.
    Returns:
        variants: pandas.DataFrame
            Table with variant data.
    """
    if log_progress:
        if wanted_properties is None:
            print(f"Loading {vcf} ...", flush=True, file=sys.stderr, end="")
        else:
            print(f"Reading {','.join(wanted_properties)} from {vcf} ...", flush=True, file=sys.stderr, end="")
    with pysam.VariantFile(vcf, 'r', threads=num_threads) as f_in:
        if samples is None:
            samples = tuple(f_in.header.samples)
        else:
            samples = set(samples)
            bad_samples = samples.difference(f_in.header.samples)
            if bad_samples:
                missing_samples_action.handle_error(f"Requested samples are not present in VCF: {bad_samples}")
            samples = tuple(sample for sample in f_in.header.samples if sample in samples)
            f_in.subset_samples(samples)
        if wanted_properties is None:
            properties = _get_all_properties(f_in.header)
        else:
            # we always want variant id, but it's easy to forget
            wanted_properties = {Keys.id}.union(wanted_properties)
            properties = _get_allowed_properties(f_in.header)
            bad_properties = wanted_properties.difference(properties.values())
            if bad_properties:
                missing_properties_action.handle_error(
                    f"Requested properties are not present in VCF: {bad_properties}\n"
                    f"Allowed_properties: {properties.values}"
                )
            properties = {key: value for key, value in properties.items() if value in wanted_properties}

        genotype_properties = {
            vcf_prop_name: table_prop_name for vcf_prop_name, table_prop_name in properties.items()
            if vcf_prop_name in f_in.header.formats
        }
        variant_properties = {
            vcf_prop_name: table_prop_name for vcf_prop_name, table_prop_name in properties.items()
            if vcf_prop_name not in genotype_properties
        }

        variant_extractors = tuple(
            _get_variant_extractor(vcf_prop_name, header=f_in.header, missing_value=missing_value)
            for vcf_prop_name in variant_properties.keys()
        )
        if samples:
            genotype_extractors = tuple(
                _get_genotype_extractor(vcf_prop_name, header=f_in.header)
                for vcf_prop_name in genotype_properties.keys()
            )

            def _row_extractor(variant_record: pysam.VariantRecord) -> Iterator[EncodedVcfField]:
                for variant_extractor in variant_extractors:
                    yield variant_extractor(variant_record)
                for genotype_record in variant_record.samples.itervalues():
                    for genotype_extractor in genotype_extractors:
                        yield genotype_extractor(genotype_record)
        else:
            def _row_extractor(variant_record: pysam.VariantRecord) -> Iterator[EncodedVcfField]:
                for variant_extractor in variant_extractors:
                    yield variant_extractor(variant_record)

        encoded_rows_iter = (
            tuple(_row_extractor(variant_record))
            for variant_record in f_in.fetch()
        )

        property_collators = tuple(
            VcfPropertyCollator.get_collator(
                vcf_prop_name=vcf_prop_name, column_name=(None, table_prop_name), header=f_in.header,
                missing_value=missing_value, int_type=int_type, float_type=float_type,
                cast_whole_nums_to_int=cast_whole_nums_to_int, ordered=category_prop_getter_ordered
            )
            for vcf_prop_name, table_prop_name in variant_properties.items()
        ) + tuple(
            VcfPropertyCollator.get_collator(
                vcf_prop_name=vcf_prop_name, column_name=(sample_id, table_prop_name), header=f_in.header,
                missing_value=missing_value, int_type=int_type, float_type=float_type,
                cast_whole_nums_to_int=cast_whole_nums_to_int, ordered=category_prop_getter_ordered
            )
            for sample_id in samples
            for vcf_prop_name, table_prop_name in genotype_properties.items()
        )

        variants = pandas.DataFrame(
            {property_collator.column_name: property_collator.get_pandas_series(property_values)
             for property_collator, property_values in zip(property_collators, zip(*encoded_rows_iter))}
        ).rename_axis(columns=(Keys.sample_id, Keys.property))\
            .sort_index(axis=1, level=Keys.sample_id)

    if variants[(None, Keys.id)].nunique() == len(variants[(None, Keys.id)]):
        # all IDs are unique, use ID as the index
        variants.set_index((None, Keys.id), inplace=True)
        variants.rename_axis(index=Keys.id, inplace=True)

    if drop_missing_columns:  # drop non-sample columns where all data is missing
        for column in set(variants.columns).difference(samples):
            if is_missing(variants[column], missing_value=missing_value).all():
                variants.drop(columns=column, inplace=True)

    if drop_trivial_multi_index:
        drop_trivial_columns_multi_index(variants)

    variants = shift_origin(variants, desired_origin=genome_origin, current_origin=vcf_origin, copy_on_change=False)

    # add custom field listing the samples
    variants.attrs[Keys.sample_ids] = samples
    if log_progress:
        print(" done", flush=True, file=sys.stderr, end='\n')
    return variants


def get_sample_ids_from_variants(variants: pandas.DataFrame) -> List[str]:
    """
    Get sample IDs from variants table
    Args:
        variants: pandas.DataFrame
            Table of variants, as loaded by vcf_to_pandas
    Returns:
        sample_ids: List[str]
            List of all samples in the variants table.
    """
    sample_ids = variants.columns.get_level_values(Keys.sample_id)
    return sample_ids[~sample_ids.isnull()].unique().to_list()


def is_missing(series: pandas.Series, missing_value: Any) -> pandas.Series:
    return series.isna() if pandas.api.types.is_numeric_dtype(series) else series.isna() | (series == missing_value)


def has_missing(series: pandas.Series, missing_value: Any) -> bool:
    return series.hasnans if pandas.api.types.is_numeric_dtype(series) else \
        series.hasnans or (series == missing_value).any()


def list_info_properties(variants: pandas.DataFrame) -> Sequence[str]:
    info_columns = variants.loc[:, (None, slice(None))].columns.droplevel(Keys.sample_id)
    return info_columns


def list_format_properties(variants: pandas.DataFrame) -> Set[str]:
    return set(variants.columns.droplevel(Keys.sample_id)).difference(list_info_properties(variants))


def get_format_property(variants: pandas.DataFrame, property_name: str) -> pandas.DataFrame:
    """ """
    try:
        prop_df = variants.loc[:, (slice(None), property_name)]
        prop_df.columns = prop_df.columns.droplevel(Keys.property)
        return prop_df
    except KeyError as key_error:
        common.add_exception_context(
            key_error,
            f"{property_name} missing from columns. Available properties: "
            f"{','.join(list_format_properties(variants))}"
        )
        raise


def get_info_property(variants: pandas.DataFrame, property_name: Union[str, Iterable[str]]) -> Union[pandas.Series, pandas.DataFrame]:
    try:
        if isinstance(property_name, str):
            return variants.loc[(None, property_name)].rename(property_name)
        else:
            prop_df = variants.loc[:, (None, list(property_name))]
            prop_df.columns = prop_df.columns.droplevel(Keys.sample_id)
            return prop_df
    except KeyError as key_error:
        common.add_exception_context(
            key_error,
            f"{property_name} missing from columns. Available properties: "
            f"{','.join(list_info_properties(variants))}"
        )
        raise


def combine_call_properties(
        preferred_properties: pandas.DataFrame,
        alternate_properties: pandas.DataFrame,
        use_alternate_properties: pandas.DataFrame
) -> pandas.DataFrame:
    """
    Combine two call properties, by selecting preferred_properties where they are non-null, and otherwise using
    alternate_properties (thus if both are non-null, preferred_properties will be selected).
    Each of the input DataFrames must have the same indices and columns
    Args:
        preferred_properties: pandas.DataFrame
            DataFrame of integer (possibly but not necessarily nullable) properties.
        alternate_properties: pandas.DataFrame
            DataFrame of integer (possibly but not necessarily nullable) properties.
        use_alternate_properties: pandas.DataFrame
            DataFrame of non-null bools specifying when to use each property
    Returns:
        combined_properties: pandas.DataFrame
            DataFrame of combined preferred and alternate properties
    """
    # the fastest way to combine is using DataFrame.where, but it can't handle promoting dtypes, so we'll have to do
    # that first
    dtypes = {dt for df in (preferred_properties, alternate_properties) for dt in df.dtypes}
    iinfos = [numpy.iinfo(dt.numpy_dtype if hasattr(dt, "numpy_dtype") else dt) for dt in dtypes]
    min_type = IntPropertyCollator.range_to_min_int_dtype(
        min_value=min(iinfo.min for iinfo in iinfos), max_value=max(iinfo.max for iinfo in iinfos), nullable=True
    )
    # convert to the minimum type that can handle all the values, then return DataFrame that is preferred_properties
    # when that is not bad, and otherwise alternate_properties
    return alternate_properties.astype(min_type).where(use_alternate_properties, preferred_properties)


def get_copy_number(variants: pandas.DataFrame, use_cn: bool = Default.use_cn) -> Optional[pandas.DataFrame]:
    """ Return table of copy number if it is available, otherwise None """
    available_props = set(variants.columns.get_level_values(Keys.property))
    if use_cn and Keys.cn in available_props:
        cn = get_format_property(variants, Keys.cn)
        if Keys.rd_cn in available_props:
            # CN *and* RD_CN are available. Sometimes one is in the header but not listed, so return whichever is
            # non-null, and if they're both non-null, prefer CN
            return combine_call_properties(
                preferred_properties=cn,
                alternate_properties=get_format_property(variants, Keys.rd_cn),
                use_alternate_properties=cn.isna()
            )
        else:
            return cn
    elif Keys.rd_cn in available_props:
        return get_format_property(variants, Keys.rd_cn)
    else:
        return None


def get_copy_number_and_quality(
        variants: pandas.DataFrame, use_cn: bool = Default.use_cn
) -> Union[Tuple[pandas.DataFrame, pandas.DataFrame], Tuple[None, None]]:
    f"""
    Return tables with copy number and associated phred-scaled quality if they are both available, otherwise None, None
    There are two possible sources of copy-number data:
        -{Keys.rd_cn} with quality {Keys.rd_gq}
        -{Keys.cn} with quality {Keys.cnq}
    By the end of CleanVcf, {Keys.cn}/{Keys.cnq} are identical to {Keys.rd_cn}/{Keys.rd_gq} when they are present, so
    these data are only examined if use_cn is set to True. In which case, the output table will choose whichever is
    available, and if both are, it will take whichever is non-null, preferring {Keys.cn} if both are non-null. It will
    always use the above-listed quality metric. If the correct quality is not available, that will cause an Exception.
    Args:
        variants: pandas.DataFrame
            Table of raw properties loaded from VCF
        use_cn: bool (default={Default.use_cn})
            If False, ignore any {Keys.cn} data. If True, reconcile it as described above.
    Returns:
            copy_number: pandas.DataFrame
                Table of copy-number with rows corresponding to variants and columns corresponding to sample IDs
            copy_number_quality:
                Corresponding table of copy-number qualities.
    """
    available_props = set(variants.columns.get_level_values(Keys.property))
    if use_cn and Keys.cn in available_props:
        if Keys.cnq not in available_props:
            raise ValueError(f"{Keys.cn} in variants but not associated quality {Keys.cnq}")
        if Keys.rd_cn in available_props:
            if Keys.rd_gq not in available_props:
                raise ValueError(f"{Keys.rd_cn} in variants but not associated quality {Keys.rd_gq}")
            cn = get_format_property(variants, Keys.cn)
            cn_is_null = cn.isna()
            cn_calls = combine_call_properties(
                preferred_properties=cn,
                alternate_properties=get_format_property(variants, Keys.rd_cn),
                use_alternate_properties=cn_is_null
            )
            cn_qualities = combine_call_properties(
                preferred_properties=get_format_property(variants, Keys.cnq),
                alternate_properties=get_format_property(variants, Keys.rd_gq),
                use_alternate_properties=cn_is_null
            )
            return cn_calls, cn_qualities
        else:
            return get_format_property(variants, Keys.cn), get_format_property(variants, Keys.cnq)
    elif Keys.rd_cn in available_props:
        if Keys.rd_gq not in available_props:
            raise ValueError(f"{Keys.rd_cn} in variants but not associated quality {Keys.rd_gq}")
        return get_format_property(variants, Keys.rd_cn), get_format_property(variants, Keys.rd_gq)
    else:
        return None, None


@common.static_vars(genotype_to_carrier_status_map={})
def genotype_to_carrier_status(genotype: Genotype) -> Union[bool, NAType]:
    if genotype in genotype_to_carrier_status.genotype_to_carrier_status_map:
        return genotype_to_carrier_status.genotype_to_carrier_status_map[genotype]
    else:
        if genotype is None:
            carrier_status = pandas.NA
        else:
            carrier_status = any(a is not None and a > 0 for a in genotype)
            if carrier_status == 0 and any(a is None for a in genotype):
                carrier_status = pandas.NA
        genotype_to_carrier_status.genotype_to_carrier_status_map[genotype] = carrier_status
        return carrier_status


def get_carrier_status(
        variants: pandas.DataFrame,
        use_copy_number: bool = Default.use_copy_number,
        use_cn: bool = Default.use_cn
) -> pandas.DataFrame:
    f"""
    return carrier status as num_variants x num_samples table, with each entry being a BooleanDType:
      - True if the sample is called non-REF for that variant
      - False if it is called REF
      - pandas.NA if it is not called
    Args:
        variants: pandas.DataFrame
            Multi-index table of raw properties loaded from VCF
        use_copy_number: bool (default={Default.use_copy_number})
            If use_copy_number is True and the genotype is not called, examine copy number info to determine carrier
            status. Otherwise return pandas.NA
        use_cn: bool (default={Default.use_cn})
            If False, ignore CN data, which duplicates RD_CN by the end of CleanVcf. If True, check CN *and* RD_CN*.
            This value has no effect if use_copy_number is False.
    Returns:
        carrier_status: pandas.DataFrame
            Normal (non-multi-index) table of carrier status
    """
    # get the carrier status as implied by the genotype
    gt_carrier_status = get_format_property(variants, Keys.gt)\
        .applymap(genotype_to_carrier_status)\
        .astype(pandas.BooleanDtype())
    if not use_copy_number:
        return gt_carrier_status
    # genotypes are often left no-call or REF for multi-allelic CNVs, so check if we can get copy number
    copy_number = get_copy_number(variants, use_cn=use_cn)
    if copy_number is None:
        return gt_carrier_status
    else:
        # Compute copy number carrier status by assuming copy number 2 is REF. Currently our pipeline sets copy number
        # to 2 for REF even on X and Y.
        # Use copy number status only where gt_carrier_status is null
        return (copy_number != 2).astype(pandas.BooleanDtype()).where(gt_carrier_status.isna(), gt_carrier_status)


def get_or_estimate_allele_frequency(
        variants: pandas.DataFrame, use_copy_number: bool = Default.use_copy_number, use_cn: bool = Default.use_cn
) -> pandas.Series:
    f"""
    If {Keys.allele_frequency} is present in the variants, just return that. Otherwise estimate it:
        -Estimate AF from cohort by counting called non-REF alleles and dividing by number of called alleles
        -If use_copy_number is true *AND* there are no called alleles (or they were all called REF but there is copy
         number info), then use copy number to estimate allele frequency. This is even dicier because we won't have
         genotype, only carrier status. Estimate allele frequency assuming Hardy_Weinberg Equilibrium.
    Args:
        variants: pandas.DataFrame
            Multi-index table of variant properties, genotypes, etc
        use_copy_number: boolean (default={Default.use_copy_number})
            If True, use copy number info (from {Keys.cn} and {Keys.rd_cn} fields) if it's helpful
            If False, only use genotypes
        use_cn: bool (default={Default.use_cn})
            If False, ignore CN data, which duplicates RD_CN by the end of CleanVcf. If True, check CN *and* RD_CN*.
            This value has no effect if use_copy_number is False.
    Returns:
        allele_frequency: pandas.Series
            Estimated allele frequency for each variant
    """
    available_props = set(variants.columns.get_level_values(Keys.property))
    if Keys.allele_frequency in available_props:
        return get_info_property(variants, Keys.allele_frequency)
    # get the allele_counts as implied by the genotype
    gt_allele_counts = get_allele_counts(variants, use_copy_number=False).sum(axis=1, skipna=True)
    gt_num_called_alleles = get_num_called_alleles(variants).sum(axis=1, skipna=True)
    gt_allele_frequency = gt_allele_counts / gt_num_called_alleles
    if not use_copy_number:
        return gt_allele_frequency
    # genotypes are often left no-call or REF for multi-allelic CNVs, so check if we can get copy number
    copy_number = get_copy_number(variants, use_cn=use_cn)
    if copy_number is None:
        return gt_allele_frequency
    else:
        # Can't really get allele count from copy number. Estimate by assuming HWE and diploid genome
        copy_number_carrier_freq = (copy_number != 2).mean(axis=1, skipna=True)
        copy_number_allele_frequency = 1.0 - (1.0 - copy_number_carrier_freq) ** 0.5
        # prefer GT based allele frequency, but where copy number estimate is better than nothing, use that
        return copy_number_allele_frequency.where(gt_allele_frequency.isna(), gt_allele_frequency)


@common.static_vars(genotype_to_allele_count_map={})
def genotype_to_allele_count(genotype: Genotype) -> Union[int, NAType]:
    if genotype in genotype_to_allele_count.genotype_to_allele_count_map:
        return genotype_to_allele_count.genotype_to_allele_count_map[genotype]
    else:
        if genotype is None:
            allele_count = pandas.NA
        else:
            allele_count = sum(1 for a in genotype if a is not None and a > 0)
            if allele_count == 0 and any(a is None for a in genotype):
                allele_count = pandas.NA
        genotype_to_allele_count.genotype_to_allele_count_map[genotype] = allele_count
        return allele_count


def get_allele_counts(
        variants: pandas.DataFrame, use_copy_number: bool = Default.use_copy_number, use_cn: bool = Default.use_cn
) -> pandas.DataFrame:
    f"""
    return allele_counts as num_variants x num_samples table, with each entry being a Uint8DType:
      - Number of non-REF counts if the sample GT is called for that variant
      - pandas.NA if it is not called
    Args:
        variants: pandas.DataFrame
            Multi-index table of raw properties loaded from VCF
        use_copy_number: bool (default={Default.use_copy_number})
            If use_copy_number is True and the genotype is not called, examine copy number info to determine carrier
            status. If the CN is called and 2, set to ref (AC=0), otherwise assume HET (AC=1). If not called,
            leave the AC as pandas.NA
        use_cn: bool (default={Default.use_cn})
            If False, ignore CN data, which duplicates RD_CN by the end of CleanVcf. If True, check CN *and* RD_CN*.
            This value has no effect if use_copy_number is False.
    Returns:
        allele_counts:
            Normal (non-multi-index) table of allele counts
    """
    gt_allele_count = get_format_property(variants, Keys.gt).applymap(genotype_to_allele_count).astype(pandas.UInt8Dtype())
    if not use_copy_number:
        return gt_allele_count
    copy_number = get_copy_number(variants, use_cn=use_cn)
    if copy_number is None:
        return gt_allele_count
    # Use copy number where GT is no-call, otherwise use GT allele counts
    return combine_call_properties(
        preferred_properties=gt_allele_count,
        alternate_properties=(copy_number != 2).astype(pandas.UInt8Dtype()),
        use_alternate_properties=gt_allele_count.isna()
    )


def get_allele_counts_and_quality(
        variants: pandas.DataFrame,
        use_copy_number: bool = Default.use_copy_number,
        use_cn: bool = Default.use_cn
) -> (pandas.DataFrame, pandas.DataFrame):
    f"""
    Form allele_counts as num_variants x num_samples table, with each entry being a Uint8DType:
      - Number of non-REF counts if the sample GT is called for that variant
      - pandas.NA if it is not called
    Also compute corresponding qualities of the calls.
    Args:
        variants: pandas.DataFrame
            Multi-index table of raw properties loaded from VCF
        use_copy_number: bool (default={Default.use_copy_number})
            If use_copy_number is True and the genotype is not called, examine copy number info to determine carrier
            status. If the copy umber is called and 2, set to ref (AC=0), otherwise assume HET (AC=1). If not called,
            leave the AC as pandas.NA. If the copy number is used, use the corresponding copy number quality rather than
            GQ.
        use_cn: bool (default={Default.use_cn})
            If False, ignore CN data, which duplicates RD_CN by the end of CleanVcf. If True, check CN *and* RD_CN*.
            This value has no effect if use_copy_number is False.
    Returns:
        allele_counts:
            Normal (non-multi-index) table of allele counts
        allele_counts_quality:
            Corresponding table of qualities of the allele counts
    """
    # get genotype-only allele count (we'll combine with copy number along-side quality)
    gt_allele_count = get_allele_counts(variants, use_copy_number=False)
    gt_quality = get_format_property(variants, Keys.gq)
    if not use_copy_number:
        return gt_allele_count, gt_quality
    copy_number, copy_number_quality = get_copy_number_and_quality(variants, use_cn=use_cn)
    if copy_number is None:
        return gt_allele_count, gt_quality
    # Use copy number where GT is no-call. Sometimes for copy-number calls, GT is set to HOMREF, so if GT is HOMREF and
    # copy number is not no-call, use copy number
    gt_allele_count_is_na = gt_allele_count.isna()
    allele_count = combine_call_properties(
        preferred_properties=gt_allele_count,
        alternate_properties=(copy_number != 2).astype(pandas.UInt8Dtype()),
        use_alternate_properties=gt_allele_count_is_na
    )
    allele_count_quality = combine_call_properties(
        preferred_properties=gt_quality,
        alternate_properties=copy_number_quality,
        use_alternate_properties=gt_allele_count_is_na
    )
    return allele_count, allele_count_quality


@common.static_vars(genotype_to_called_allele_count_map={})
def genotype_to_called_allele_count(genotype: Genotype) -> int:
    if genotype in genotype_to_called_allele_count.genotype_to_called_allele_count_map:
        return genotype_to_called_allele_count.genotype_to_called_allele_count_map[genotype]
    else:
        if genotype is None:
            called_allele_count = 0
        else:
            called_allele_count = sum(1 for a in genotype if a is not None and a > 0)
        genotype_to_called_allele_count.genotype_to_called_allele_count_map[genotype] = called_allele_count
        return called_allele_count


def get_called_allele_count(variants: pandas.DataFrame) -> pandas.DataFrame:
    return get_format_property(variants, Keys.gt).applymap(genotype_to_called_allele_count).astype(numpy.uint8)


@common.static_vars(genotype_to_num_called_alleles_map={})
def genotype_to_num_called_alleles(genotype: Genotype) -> int:
    if genotype in genotype_to_num_called_alleles.genotype_to_num_called_alleles_map:
        return genotype_to_num_called_alleles.genotype_to_num_called_alleles_map[genotype]
    else:
        if genotype is None:
            num_called_alleles = 0
        else:
            num_called_alleles = sum(1 for a in genotype if a is not None)
        genotype_to_num_called_alleles.genotype_to_num_called_alleles_map[genotype] = num_called_alleles
        return num_called_alleles


def get_num_called_alleles(variants: pandas.DataFrame) -> pandas.DataFrame:
    return get_format_property(variants, Keys.gt).applymap(genotype_to_num_called_alleles).astype(numpy.uint8)


def sort_multi_index_dataframe_columns(df: pandas.DataFrame):
    df.sort_index(
        axis=1, inplace=True, level=[0, 1],
        key=lambda ind: pandas.Index(["" if isinstance(name, float) and numpy.isnan(name) else name for name in ind])
    )


def assign_dataframe_column(
        df: pandas.DataFrame,
        values: Union[pandas.DataFrame, pandas.Series],
        name: str
) -> pandas.DataFrame:
    """ helper function to assign a column values to a DataFrame, potentially with multi-index columns """
    if isinstance(values, pandas.Series):
        # adding a single column (e.g. adding an info field to a variants DataFrame)
        if isinstance(df.columns, pandas.MultiIndex):
            values.name = (None, name)
            is_new_prop = name in df.columns.get_level_values(level=Keys.property)
            if is_new_prop:
                df = df.join(values)
            else:
                df.loc[:, (None, name)] = values
            # need to sort to put the table back into well-organized form
            sort_multi_index_dataframe_columns(df)
        else:
            # Adding series to non-multi-index DataFrame, the usual deal:
            df[name] = values
    else:
        # adding a multi-index property (e.g. adding a genotype field to a variants DataFrame)
        if not isinstance(df.columns, pandas.MultiIndex):
            # adding a multi-index property to a regular DataFrame: need to promote df to multi-index
            df.columns = pandas.MultiIndex.from_product([[None], df.columns], names=[Keys.sample_id, Keys.property])

        values.columns = pandas.MultiIndex.from_product([values.columns, [name]], names=[Keys.sample_id, Keys.property])
        is_new_prop = name not in df.columns.get_level_values(level=Keys.property)
        if is_new_prop:
            df = df.join(values)
        else:
            df.loc[:, (slice(None), name)] = values
        # need to sort to put the table back into well-organized form
        sort_multi_index_dataframe_columns(df)

    return df
