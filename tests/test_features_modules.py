import pytest

from fragmentomics import Fragment, dict_ref_lookup
from fragmentomics.lengths import binned_short_long_ratio, length_features
from fragmentomics.motifs import end_motif_counts, motif_features
from fragmentomics.wps import Region, call_nucleosomes, wps_track


def test_fragment_geometry():
    f = Fragment("chr1", 100, 267)
    assert f.length == 167
    assert f.midpoint == pytest.approx(183.5)


def test_length_features_basic():
    frags = [Fragment("chrT", 0, n) for n in (120, 130, 140, 160, 170, 200)]
    feats = length_features(frags)
    assert feats.n_fragments == 6
    # short=100-150 -> {120,130,140}=3 ; long=151-220 -> {160,170,200}=3
    assert feats.short_long_ratio == pytest.approx(1.0)
    assert 100 <= feats.median <= 200


def test_length_features_empty():
    feats = length_features([])
    assert feats.n_fragments == 0
    assert feats.short_long_ratio != feats.short_long_ratio  # nan


def test_binned_ratio_keys():
    frags = [Fragment("chrT", 1_000, 1_130), Fragment("chrT", 2_000, 2_180)]
    out = binned_short_long_ratio(frags, bin_size=1_000_000)
    assert all(isinstance(k, tuple) and len(k) == 2 for k in out)


def test_end_motif_counts_both_strands():
    # reference where we control the ends precisely
    ref = dict_ref_lookup({"chrT": "ACGT" + "N" * 10 + "TTTT" + "A" * 80})
    # fragment from 0..? ; watson end at 0 -> "ACGT"
    f = Fragment("chrT", 0, 18)  # end=18 -> ref[14:18] = "TTTT", revcomp -> "AAAA"
    counts = end_motif_counts([f], ref, k=4)
    assert counts["ACGT"] == 1
    assert counts["AAAA"] == 1


def test_motif_features_mds_range(biased_sample):
    frags, ref_dict, _ = biased_sample
    ref = dict_ref_lookup(ref_dict)
    feats = motif_features(frags, ref, k=4)
    assert feats.n_ends > 0
    assert 0.0 <= feats.mds <= 1.0
    assert len(feats.top_motifs) == 10


def test_wps_spanning_positive_center():
    # one long fragment fully spanning a central window -> positive WPS in middle
    region = Region("chrT", 0, 400)
    frags = [Fragment("chrT", 100, 280)]  # length 180, within [120,180]
    track = wps_track(frags, region, window=120)
    assert track.shape[0] == region.length
    # near the fragment centre the window is spanned (WPS = +1)
    assert track[190] == pytest.approx(1.0)
    # far outside the fragment WPS is 0
    assert track[10] == pytest.approx(0.0)


def test_call_nucleosomes_recovers_positions(biased_sample):
    frags, _, true_dyads = biased_sample
    region = Region("chrT", 0, 20_000)
    track = wps_track(frags, region, window=120)
    called = call_nucleosomes(track, region, min_distance=150)
    assert len(called) > 5
    # each called peak should be reasonably close to some true dyad
    for c in called:
        nearest = min(true_dyads, key=lambda d: abs(d - c))
        assert abs(c - nearest) < 60
