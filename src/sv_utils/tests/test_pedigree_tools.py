import os
from typing import Optional, Sequence, Tuple, Set

from sv_utils import pedigree_tools


class Default:
    resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    ped_file = os.path.join(resources_dir, "1KGP_2504_and_698_with_GIAB.ped")
    random_seed = 0


class SubsetData:
    # manually curated subset data for test_subset_participants
    p_ids = ("HG04204", "HG03679", "HG04215", "HG03642", "HG04140", "HG04141", "HG04142", "HG04148", "HG04147",
             "HG04146")
    num_trios = 4
    num_families = 3
    trio_inds = (
        pedigree_tools.TrioIndices(mother=3, father=1, proband=0),
        pedigree_tools.TrioIndices(mother=3, father=1, proband=2),
        pedigree_tools.TrioIndices(mother=5, father=4, proband=6),
        pedigree_tools.TrioIndices(mother=8, father=9, proband=7)
    )


def _trio_inds_to_tuples(trio_inds: Tuple[pedigree_tools.TrioIndices, ...]) -> Set[Tuple[int, int, int]]:
    trio_tuples = set(tuple(trio_ind) for trio_ind in trio_inds)
    # better not be any repeats
    assert len(trio_tuples) == len(trio_inds)
    return trio_tuples


def _assert_trio_inds_equal(trio_inds_a: Tuple[pedigree_tools.TrioIndices, ...],
                            trio_inds_b: Tuple[pedigree_tools.TrioIndices, ...]):
    # we don't care about the order the trios are listed in, only that the same trios are iterated over
    assert not _trio_inds_to_tuples(trio_inds_a).symmetric_difference(_trio_inds_to_tuples(trio_inds_b))


def test_subset_participants(ped_file: str = Default.ped_file,
                             subset_pids: Sequence[pedigree_tools.ParticipantId] = SubsetData.p_ids,
                             num_subset_trios: int = SubsetData.num_trios,
                             num_subset_families: int = SubsetData.num_families,
                             trio_inds: Tuple[pedigree_tools.TrioIndices, ...] = SubsetData.trio_inds):
    pedigree_file_info = pedigree_tools.PedigreeFileInfo.load(ped_file)
    num_trios = pedigree_file_info.num_trios
    num_families = pedigree_file_info.num_families
    num_participants = pedigree_file_info.num_participants
    original_participants = pedigree_file_info.participant_ids
    original_families = pedigree_file_info.parents_ids

    # get the TrioIndices corresponding to these participants
    _assert_trio_inds_equal(pedigree_file_info.get_trio_indices(subset_pids), trio_inds)

    # check that the restrictions produce the expected results
    subset_pedigree_file_info = pedigree_file_info.subset_participants(subset_pids)
    assert subset_pedigree_file_info.participant_ids == set(subset_pids)
    assert subset_pedigree_file_info.num_trios == num_subset_trios
    assert subset_pedigree_file_info.num_families == num_subset_families

    # check it's a true subset, we haven't somehow corrupted participant IDs or family IDs
    assert all(participant_id in original_participants for participant_id in subset_pedigree_file_info.participant_ids)
    assert all(family_id in original_families for family_id in subset_pedigree_file_info.parents_ids)
    # check we didn't change the original PedigreeFileInfo
    assert pedigree_file_info.num_trios == num_trios
    assert pedigree_file_info.num_families == num_families
    assert pedigree_file_info.num_participants == num_participants
    assert pedigree_file_info.participant_ids == original_participants
    assert pedigree_file_info.parents_ids == original_families


def test_subset_num_families(ped_file: str = Default.ped_file, random_seed: Optional[int] = Default.random_seed):
    pedigree_file_info = pedigree_tools.PedigreeFileInfo.load(ped_file)
    num_trios = pedigree_file_info.num_trios
    num_families = pedigree_file_info.num_families
    num_participants = pedigree_file_info.num_participants
    original_participants = pedigree_file_info.participant_ids
    original_families = pedigree_file_info.parents_ids
    assert len(original_participants) == num_participants
    assert len(original_families) == num_families

    num_subset_families = num_families // 2
    subset_pedigree_file_info = pedigree_file_info.subset_num_families(num_families=num_subset_families,
                                                                       random_state=random_seed)
    assert isinstance(subset_pedigree_file_info, pedigree_tools.PedigreeFileInfo)
    assert subset_pedigree_file_info.num_families == num_subset_families
    # not every family is a trio (we have some quads). Can test the trios are within possible bounds
    assert num_subset_families <= subset_pedigree_file_info.num_trios <= 2 * num_subset_families
    assert 3 * num_subset_families <= subset_pedigree_file_info.num_participants <= 4 * num_subset_families
    # check it's a true subset, we haven't somehow corrupted participant IDs or family IDs
    assert all(participant_id in original_participants for participant_id in subset_pedigree_file_info.participant_ids)
    assert all(family_id in original_families for family_id in subset_pedigree_file_info.parents_ids)
    # check we didn't change the original PedigreeFileInfo
    assert pedigree_file_info.num_trios == num_trios
    assert pedigree_file_info.num_families == num_families
    assert pedigree_file_info.num_participants == num_participants
    assert pedigree_file_info.participant_ids == original_participants
    assert pedigree_file_info.parents_ids == original_families
