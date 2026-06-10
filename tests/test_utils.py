
import pytest

from fragmentomics.utils import (
    jensen_shannon_divergence,
    log2_enrichment,
    normalised_entropy,
    reverse_complement,
    shannon_entropy,
)


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GATTACA") == "TGTAATC"
    assert reverse_complement("ACGTN") == "NACGT"


def test_shannon_entropy_uniform_vs_peaked():
    uniform = shannon_entropy([1, 1, 1, 1])
    peaked = shannon_entropy([97, 1, 1, 1])
    assert uniform == pytest.approx(2.0)  # log2(4)
    assert peaked < uniform


def test_normalised_entropy_bounds():
    assert normalised_entropy([1, 1, 1, 1], 4) == pytest.approx(1.0)
    assert normalised_entropy([1, 0, 0, 0], 4) == pytest.approx(0.0)
    assert 0.0 < normalised_entropy([5, 3, 1, 1], 4) < 1.0


def test_jsd_identical_is_zero_and_symmetric():
    p = {"AA": 3, "CC": 1}
    q = {"AA": 6, "CC": 2}  # same proportions
    assert jensen_shannon_divergence(p, q) == pytest.approx(0.0, abs=1e-9)
    a = {"AA": 9, "TT": 1}
    b = {"AA": 1, "TT": 9}
    assert jensen_shannon_divergence(a, b) == pytest.approx(
        jensen_shannon_divergence(b, a)
    )
    assert jensen_shannon_divergence(a, b) > 0


def test_jsd_handles_empty():
    assert jensen_shannon_divergence({}, {"AA": 1}) == 0.0


def test_log2_enrichment_direction():
    fg = {"CGCG": 100, "AAAA": 1}
    bg = {"CGCG": 1, "AAAA": 100}
    enr = log2_enrichment(fg, bg)
    assert enr["CGCG"] > 0
    assert enr["AAAA"] < 0
