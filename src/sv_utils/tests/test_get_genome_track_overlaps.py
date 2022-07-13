import os
import glob
import pytest
from typing import Iterable

from sv_utils import get_genome_track_overlaps, genomics_io


class Default:
    resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    small_vcfs_dir = os.path.join(resources_dir, "small_vcfs")
    small_vcfs = tuple(sorted(glob.glob(os.path.join(small_vcfs_dir, "*.vcf.gz"))))
    small_vcf = small_vcfs[0]
    tracks_dir = os.path.join(resources_dir, "tracks")
    genome_tracks = tuple(glob.glob(os.path.join(tracks_dir, "hg38-*.bed.gz")))


@pytest.mark.integration_test
def test_get_genome_track_overlaps(
        vcf: str = Default.small_vcf, genome_tracks: Iterable[str] = Default.genome_tracks
):
    overlaps = get_genome_track_overlaps.get_genome_track_overlaps(vcf=vcf, genome_track_files=genome_tracks)
    # basic sanity checks
    #     -no negative overlaps, some positive overlaps
    assert (overlaps >= 0).all(axis=None)
    assert (overlaps > 0).any(axis=None)
    #     -all genome tracks are present in the columns
    for track in genome_tracks:
        track_name = os.path.basename(track).split('.', maxsplit=1)[0]
        assert track_name in overlaps.columns
    #     -the overlap index has all the variants (it is the same as the index of the vcf)
    variant_ids = genomics_io.get_vcf_variant_ids(vcf)
    assert overlaps.index.equals(variant_ids)
