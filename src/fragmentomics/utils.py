"""Small, dependency-light numerical and sequence helpers.

These are kept free of any I/O so they are trivially unit-testable and reusable
across the feature modules.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Mapping

import numpy as np

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA string (N-safe)."""
    return seq.translate(_COMPLEMENT)[::-1]


def shannon_entropy(counts: Iterable[float], base: float = 2.0) -> float:
    """Shannon entropy of a (possibly unnormalised) count vector.

    Zero-count categories contribute nothing. Returns 0.0 for an empty/degenerate
    distribution.
    """
    arr = np.asarray(list(counts), dtype=float)
    total = arr.sum()
    if total <= 0:
        return 0.0
    p = arr[arr > 0] / total
    return float(-(p * (np.log(p) / np.log(base))).sum())


def normalised_entropy(counts: Iterable[float], n_categories: int, base: float = 2.0) -> float:
    """Entropy normalised to [0, 1] by the maximum possible entropy.

    This is the basis of the Motif Diversity Score (MDS) of Jiang et al. 2020,
    where ``n_categories == 4 ** k`` for k-mer end motifs.
    """
    if n_categories <= 1:
        return 0.0
    h = shannon_entropy(counts, base=base)
    return h / (math.log(n_categories) / math.log(base))


def jensen_shannon_divergence(
    p: Mapping[str, float] | Iterable[float],
    q: Mapping[str, float] | Iterable[float],
    base: float = 2.0,
) -> float:
    """Jensen-Shannon divergence between two distributions.

    Accepts either aligned count/probability vectors or dict-like mappings keyed
    by category (e.g. k-mer string). Mappings are aligned on the union of keys so
    that callers do not have to pre-align motif vocabularies. Result is in [0, 1]
    when ``base == 2`` (square root would give the bounded metric; we return the
    divergence itself, which is the quantity of interest for coupling strength).
    """
    if isinstance(p, Mapping) or isinstance(q, Mapping):
        keys = sorted(set(p) | set(q))  # type: ignore[arg-type]
        pv = np.array([float(p.get(k, 0.0)) for k in keys])  # type: ignore[union-attr]
        qv = np.array([float(q.get(k, 0.0)) for k in keys])  # type: ignore[union-attr]
    else:
        pv = np.asarray(list(p), dtype=float)
        qv = np.asarray(list(q), dtype=float)

    if pv.sum() <= 0 or qv.sum() <= 0:
        return 0.0
    pv = pv / pv.sum()
    qv = qv / qv.sum()
    m = 0.5 * (pv + qv)

    def _kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0
        return float((a[mask] * (np.log(a[mask] / b[mask]) / np.log(base))).sum())

    return 0.5 * _kl(pv, m) + 0.5 * _kl(qv, m)


def log2_enrichment(
    foreground: Mapping[str, float],
    background: Mapping[str, float],
    pseudocount: float = 1.0,
) -> dict[str, float]:
    """Per-category log2 fold enrichment of foreground over background.

    Both inputs are normalised internally; a pseudocount guards against zeros.
    """
    keys = sorted(set(foreground) | set(background))
    fg = np.array([foreground.get(k, 0.0) + pseudocount for k in keys])
    bg = np.array([background.get(k, 0.0) + pseudocount for k in keys])
    fg = fg / fg.sum()
    bg = bg / bg.sum()
    return {k: float(np.log2(f / b)) for k, f, b in zip(keys, fg, bg, strict=True)}
