import collections.abc
import numpy
from typing import Union, Iterable, Tuple, List, Set, Dict, Text, Iterator, Optional, Sequence, Collection
from enum import Enum

from sv_utils import common


class TrioIndices(collections.abc.Iterable, collections.abc.Hashable):
    """
    Simple class for holding sample indices pointing to members of a trio
    """
    __slots__ = ("mother", "father", "proband")

    def __init__(
            self,
            mother: int,
            father: int,
            proband: int
    ):
        self.mother = mother
        self.father = father
        self.proband = proband

    def __hash__(self):
        return hash(tuple(self))

    def __iter__(self):
        yield self.mother
        yield self.father
        yield self.proband

    def __repr__(self):
        return f"({self.mother},{self.father},{self.proband})"


class Sex(Enum):
    Unknown = '0'
    Male = '1'
    Female = '2'

    def __str__(self):
        return self.value


class Phenotype(Enum):
    Missing = '0'
    Unaffected = '1'
    Affected = '2'

    @classmethod
    def _missing_(cls, value):
        if value == "-9":
            return cls.Missing
        else:
            raise ValueError("Invalid Phenotype: %s" % str(value))

    def __str__(self):
        return self.value


class ParticipantId(str):
    UNKNOWN: "ParticipantId" = '0'

    @property
    def is_unknown(self):
        return self == self.UNKNOWN


class FamilyId(str):
    UNKNOWN: "FamilyId" = '0'

    @property
    def is_unknown(self):
        return self == self.UNKNOWN


class PedigreeLine:
    """
    Basic class to read in / store info from one line of a pedigree file. Can alternatively initialize by specifying
    the values of individual fields as strings.
    Checks during object creation:
        Family and proband IDs are asserted to be known
    """
    __slots__ = ("family_id", "proband_id", "father_id", "mother_id", "sex", "phenotype", "provenance")
    unique_sep: str = '\t'

    family_id: FamilyId
    proband_id: ParticipantId
    father_id: ParticipantId
    mother_id: ParticipantId
    sex: Sex
    phenotype: Phenotype
    provenance: Text

    def __init__(
            self,
            family_id: Union[FamilyId, Text],
            proband_id: Union[ParticipantId, Text],
            father_id: Union[ParticipantId, Text],
            mother_id: Union[ParticipantId, Text],
            sex: Union[Sex, Text],
            phenotype: Union[Phenotype, Text],
            provenance: Text
    ):
        self.family_id = FamilyId(family_id)
        self.proband_id = ParticipantId(proband_id)
        self.father_id = ParticipantId(father_id)
        self.mother_id = ParticipantId(mother_id)

        # set sex and phenotype to appropriate enumerations
        self.sex = Sex(sex)
        self.phenotype = Phenotype(phenotype)
        self.provenance = provenance

        if self.family_id.is_unknown:
            raise ValueError(f"While processing {provenance}: Unknown family ID")
        if self.proband_id.is_unknown:
            raise ValueError(f"While processing {provenance}: Unknown proband ID for family {family_id}")

    @staticmethod
    def from_ped_file_line(pedigree_line: Text, file_name: Text, line_number: int):
        return PedigreeLine(*pedigree_line.split(maxsplit=7)[:6], provenance=f"{file_name}: line {line_number}")

    def __str__(self):
        # re-create the original line, minus the provenance
        return '\t'.join(str(getattr(self, _slot)) for _slot in self.__slots__[:-1])

    @property
    def any_unknown(self) -> bool:
        return self.father_id.is_unknown or self.mother_id.is_unknown or self.proband_id.is_unknown

    @property
    def parents_id(self) -> FamilyId:
        return FamilyId.UNKNOWN if self.father_id.is_unknown or self.mother_id.is_unknown \
            else FamilyId(self.unique_sep.join((self.father_id, self.mother_id)))

    @property
    def trio_ids(self) -> Tuple[ParticipantId, ParticipantId, ParticipantId]:
        """ return participant IDs of trio in order: proband, father, mother """
        return self.proband_id, self.father_id, self.mother_id

    def isolate_proband_from_parents(self) -> "PedigreeLine":
        """
        Return info in this pedigree file line altered so that the proband is isolated from the family:
            set father_id and mother_id to UNKNOWN
        Useful for taking subsets of data to attempt to estimate power.
        Returns:
            new_pedigree_line_info: PedigreeLine
                Pedigree line info corresponding to isolated proband
        """
        # copy self:
        new_pedigree_line_info = PedigreeLine(*tuple(getattr(self, _slot) for _slot in self.__slots__))
        # update:
        new_pedigree_line_info.father_id = ParticipantId.UNKNOWN
        new_pedigree_line_info.mother_id = ParticipantId.UNKNOWN
        new_pedigree_line_info.provenance += " isolated-proband"
        return new_pedigree_line_info


class PedigreeFileInfo:
    """
    Class to read / store / write data in pedigree file.
    Checks during object creation:
        Family IDs are asserted to be consistent with parent IDs
        Each participant is asserted to be a proband in at most one trio
    Notes:
        IDs are participant IDs, *not* sample IDs (e.g. it's possible for one participant to have multiple samples)
    """
    __slots__ = ("_pedigree_lines", "parents_line_indices", "participant_indices", "proband_indices",
                 "_num_indexed_lines")
    _pedigree_lines: List[PedigreeLine]
    parents_line_indices: Dict[FamilyId, Tuple[int, ...]]  # map from family id to relevant indices in pedigree_lines
    # map from unique participant id to relevant indices in pedigree_lines:
    participant_indices: Dict[ParticipantId, Tuple[int, ...]]
    # map from unique participant id to index in _pedigree lines where that participant is a proband:
    proband_indices: Dict[ParticipantId, int]
    _num_indexed_lines: int  # number of pedigree lines that have been indexed

    def __init__(self, pedigree_lines: Iterable[PedigreeLine]):
        self.pedigree_lines = pedigree_lines  # note: this will trigger setting of indices

    @staticmethod
    def _read_pedigree_file(pedigree_file: str) -> Iterator[PedigreeLine]:
        """ Private method that reads pedigree files and yields PedigreeLines, used internally by the class """
        with open(pedigree_file, 'r') as f_in:
            for line_number, pedigree_file_line in enumerate(f_in):
                if pedigree_file_line.startswith('#'):
                    continue
                try:
                    pedigree_line = PedigreeLine.from_ped_file_line(pedigree_file_line, pedigree_file, line_number)
                except Exception as err:
                    message = "Error reading line %d of pedigree file: %s" \
                              % (line_number + 1, pedigree_file)
                    common.add_exception_context(err, message)
                    raise
                yield pedigree_line

    @staticmethod
    def load(pedigree_files: Union[str, Iterable[str]]) -> "PedigreeFileInfo":
        """ Load one or more pedigree files from disk, and return PedigreeFileInfo object """
        if isinstance(pedigree_files, str):
            pedigree_files = (pedigree_files,)
        return PedigreeFileInfo([
            pedigree_line
            for pedigree_file in pedigree_files
            for pedigree_line in PedigreeFileInfo._read_pedigree_file(pedigree_file)
        ])

    @property
    def pedigree_lines(self) -> List[PedigreeLine]:
        """ Getter for pedigree_lines """
        return self._pedigree_lines

    @pedigree_lines.setter
    def pedigree_lines(self, pedigree_lines: Iterable[PedigreeLine]):
        """
        Setter for pedigree_lines, it ensures the internally-managed indices are up to date when pedigree_lines are
        changed
        """
        # when pedigree_lines are set, update indices
        if not isinstance(pedigree_lines, Iterable):
            raise ValueError(
                f"pedigree_lines should only be set to an Iterable[PedigreeLine] (was set to {type(pedigree_lines)})"
            )
        self.parents_line_indices = {}
        self.participant_indices = {}
        self.proband_indices = {}
        self._num_indexed_lines = 0
        self._pedigree_lines = list(pedigree_lines)
        for pedigree_line in self._pedigree_lines:
            if not isinstance(pedigree_line, PedigreeLine):
                raise ValueError(
                    f"pedigree_lines should only be set to an Iterable[PedigreeLine] (was set to Iterable["
                    f"{type(pedigree_line)}])"
                )
        self._update_indices()

    def extend(self, new_pedigree_lines: List[PedigreeLine]):
        """ Method to add new pedigree_lines to existing PedigreeFileInfo object """
        self._pedigree_lines.extend(new_pedigree_lines)
        self._update_indices()

    def _update_indices(self):
        """ Private method that manages various internal indices when pedigree_lines changes """
        for line_number, pedigree_line in enumerate(self.pedigree_lines[self._num_indexed_lines:],
                                                    start=self._num_indexed_lines):
            # keep track of line numbers that each set of parents_line_indices, participant, and proband are in
            self.parents_line_indices[pedigree_line.parents_id] = \
                self.parents_line_indices.get(pedigree_line.parents_id, ()) + (line_number,)
            for participant_id in pedigree_line.trio_ids:
                if participant_id.is_unknown:
                    continue
                self.participant_indices[participant_id] = \
                    self.participant_indices.get(participant_id, ()) + (line_number,)
            if pedigree_line.proband_id in self.proband_indices:
                old_index = self.proband_indices[pedigree_line.proband_id]
                old_line = self._pedigree_lines[old_index]
                raise ValueError(
                    f"{pedigree_line.proband_id} is a proband in multiple places: {old_line.provenance} line "
                    f"{old_index}, and {pedigree_line.provenance} line {line_number}"
                )
            self.proband_indices[pedigree_line.proband_id] = line_number
        self._num_indexed_lines = len(self.pedigree_lines)

    @property
    def num_families(self) -> int:
        """
        Return number of unique families: in this case, "family" refers to a unique pair of parents with one or more
        offspring, not an extended family, e.g. grandparents
        """
        return len(self.parents_ids)

    def subset_participants(self, participant_ids: Collection[ParticipantId]) -> "PedigreeFileInfo":
        """
        Return subset of these data by including only the supplied participant IDs. Note that some participant IDs
        """
        allowed_participant_ids = set(participant_ids).union({ParticipantId.UNKNOWN})
        return PedigreeFileInfo(
            [pedigree_line for pedigree_line in self.pedigree_lines
             if all(participant_id in allowed_participant_ids for participant_id in pedigree_line.trio_ids)]
        )

    def subset_families(self, parents_ids: Collection[FamilyId]) -> "PedigreeFileInfo":
        """
        Return subset of these data by including only the lines corresponding to the supplied family IDs
        """
        return PedigreeFileInfo(
            [pedigree_line
             for parents_id in parents_ids
             for pedigree_line in common.deref_iterable_with_sorted_indices(self.pedigree_lines,
                                                                            self.parents_line_indices[parents_id])
             ]
        )

    def subset_num_families(
            self,
            num_families: int,
            random_state: Union[None, int, numpy.random.RandomState] = None
    ) -> "PedigreeFileInfo":
        """
        Return subset of these data by randomly selecting the requested number of families, and only including the
        relevant pedigree lines. This is mainly useful for power analysis, to help determine how results scale with the
        number of families in the data set.
        Args:
            num_families: int
                Number of families to keep. Will raise ValueError if more than the number of families available.
            random_state: numpy.random.RandomState, int or None
                If None, will create a random seed, otherwise, use the supplied seed.
        Returns:
            subset_pedigree_file_info: PedigreeFileInfo
        """
        if num_families >= self.num_families:
            if num_families == self.num_families:
                return self
            else:
                raise ValueError(
                    f"Requested num_families ({num_families}) greater than number in data set ({self.num_families}"
                )
        if not isinstance(random_state, numpy.random.RandomState):
            random_state = numpy.random.RandomState(random_state)
        selected_families = random_state.choice(sorted(self.parents_ids), size=num_families, replace=False)
        return self.subset_families(selected_families)

    @property
    def num_trios(self) -> int:
        """
        Return number of trios where every participant ID is known
        """
        return sum(1 for pedigree_line in self.pedigree_lines if not pedigree_line.any_unknown)

    @property
    def participant_ids(self) -> Set[ParticipantId]:
        """ Return the set of all (known) participant IDs """
        return set(self.participant_indices.keys())

    @property
    def parents_ids(self) -> Set[FamilyId]:
        """ Return the set of all family IDs where there is at least one valid trio with no unknown participant IDs """
        return {
            parents_id
            for parents_id, parents_inds in self.parents_line_indices.items()
            if not all(self.pedigree_lines[parents_ind].any_unknown for parents_ind in parents_inds)
        }

    @property
    def num_participants(self) -> int:
        """ Return the number of known participant IDs """
        return len(self.participant_indices)

    def get_family_sizes(self) -> Dict[int, Tuple[FamilyId, ...]]:
        """
        Return a dictionary with a mapping from family size to tuple of FamilyIds listing the families with this size.
        This function is useful for down-sampling a data set, maintaining families of a given size (quads, tuples, etc).
        Doesn't consider data from pedigree lines with any unknown participants
        """
        sizes_dict: Dict[int, Tuple[FamilyId, ...]] = {}
        for parents_id, family_participant_ids in self.iter_family_participant_ids():
            family_size = len(family_participant_ids)
            sizes_dict[family_size] = sizes_dict.get(family_size, ()) + (parents_id,)
        return sizes_dict

    def iter_family_participant_ids(self) -> Iterator[Tuple[str, Set[str]]]:
        """
        Iterate over each nuclear family, and yield the parents_id, and set of all the participant IDs for each family.
        Useful for down-sampling / filtering families in data-sets
        """
        for parents_id, parents_inds in self.parents_line_indices.items():
            # Note that anyone with both parents unknown will be lumped into a massive "family". That line will be
            # filtered out by the call to _get_family_participant_ids though
            family_participant_ids = self._get_family_participant_ids(parents_id, parents_inds)
            if family_participant_ids:
                yield parents_id, family_participant_ids

    def _get_family_participant_ids(
            self,
            parents_id: FamilyId,
            parents_inds: Optional[Iterable[int]] = None
    ) -> Set[ParticipantId]:
        """
        Private method: for a given parents_id, return the set of unique participants in that family.
        Args:
            parents_id: str
                Unique nuclear family identifier (father_id + "\t" + mother_id)
            parents_inds: Optional[Iterable[int]] (default=None)
                If provided, must be sorted list of ints listing the indices to pedigree_lines that concern this family
                If not provided, it's looked up via self.parents_line_indices
        """
        if parents_inds is None:
            parents_inds = self.parents_line_indices[parents_id]
        return {
            participant_id
            for parents_line in common.deref_iterable_with_sorted_indices(self.pedigree_lines, parents_inds)
            for participant_id in parents_line.trio_ids
            if not parents_line.any_unknown
        }

    def get_trio_indices(self, participant_ids: Sequence[ParticipantId]) -> Tuple[TrioIndices, ...]:
        """
        Given an ordered list of participant_ids, return Tuple of TrioIndices corresponding to the indices of each
        member of a trio in the data set. Useful for functions that need to loop over trio properties that are part of
        a matrix.
        Args:
            participant_ids: Sequence[ParticipantId]
                The participant ids to get TrioIndices from
        """
        check_ids = set(participant_ids)
        return tuple(
            TrioIndices(mother=participant_ids.index(pedigree_line.mother_id),
                        father=participant_ids.index(pedigree_line.father_id),
                        proband=participant_ids.index(pedigree_line.proband_id))
            for pedigree_line in self.pedigree_lines
            if (not pedigree_line.any_unknown) and
            all(participant_id in check_ids for participant_id in pedigree_line.trio_ids)
        )
