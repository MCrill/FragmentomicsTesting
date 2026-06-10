"""Tests for the NuPEM coupling metric.

The key behavioural contract: when end-motif context depends on nucleosome phase
(``phase_bias=True``), the coupling score must be substantially higher than when it
does not (``phase_bias=False``). This is the property that makes NuPEM meaningful.
"""

from fragmentomics import Fragment, Region, dict_ref_lookup
from fragmentomics.features import compute_sample_features
from fragmentomics.nupem import nupem
from fragmentomics.wps import call_nucleosomes, wps_track


def _call_dyads(frags, chrom="chrT", length=20_000):
    region = Region(chrom, 0, length)
    track = wps_track(frags, region, window=120)
    return {chrom: call_nucleosomes(track, region, min_distance=150)}


def test_nupem_detects_phase_coupling(biased_sample, unbiased_sample):
    bf, bref, _ = biased_sample
    uf, uref, _ = unbiased_sample

    biased = nupem(bf, _call_dyads(bf), dict_ref_lookup(bref), k=4)
    unbiased = nupem(uf, _call_dyads(uf), dict_ref_lookup(uref), k=4)

    assert biased.n_core_ends > 0 and biased.n_linker_ends > 0
    # injected coupling should be clearly detectable above the unbiased baseline
    assert biased.coupling_score > unbiased.coupling_score
    assert biased.coupling_score > 0.05


def test_nupem_recovers_injected_motifs(biased_sample):
    frags, ref_dict, _ = biased_sample
    res = nupem(frags, _call_dyads(frags), dict_ref_lookup(ref_dict), k=4)
    core_motifs = {m for m, _ in res.core_top_enriched}
    linker_motifs = {m for m, _ in res.linker_top_enriched}
    # we stamped CGCG near dyads and TTTT in linkers
    assert "CGCG" in core_motifs
    assert "TTTT" in linker_motifs


def test_nupem_empty_dyads_is_safe():
    frags = [Fragment("chrT", 100, 267)]
    ref = dict_ref_lookup({"chrT": "ACGT" * 100})
    res = nupem(frags, {"chrT": []}, ref)
    assert res.coupling_score == 0.0
    assert res.n_core_ends == 0


def test_compute_sample_features_end_to_end(biased_sample):
    frags, ref_dict, _ = biased_sample
    ref = dict_ref_lookup(ref_dict)
    regions = [Region("chrT", 0, 20_000)]
    feats = compute_sample_features("SAMPLE_01", frags, ref, regions, k=4)

    d = feats.to_dict()
    assert d["sample_id"] == "SAMPLE_01"
    assert d["n_fragments"] == len(frags)
    assert d["n_nucleosomes"] > 0
    assert 0.0 <= d["motifs"]["mds"] <= 1.0
    assert d["nupem"]["coupling_score"] >= 0.0
    assert "lengths" in d and "motifs" in d and "nupem" in d
