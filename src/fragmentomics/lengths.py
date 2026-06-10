"""Fragment-length features.

Two complementary views:

* the full length *histogram* and summary statistics, and
* a DELFI-style *short:long ratio* (Cristiano et al. 2019), which is the workhorse
  of genome-wide fragmentation profiling and is sensitive to the mononucleosomal
  vs sub-nucleosomal balance that shifts in tumour-derived cfDNA.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

import numpy as np

from .fragments import Fragment


@dataclass
class LengthFeatures:
    n_fragments: int
    mean: float
    median: float
    mode: int
    short_long_ratio: float
    frac_sub_nucleosomal: float  # fraction < 100 bp
    frac_mononucleosomal: float  # fraction in [100, 220]
    histogram: dict[int, int] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = self.__dict__.copy()
        # keep the JSON compact: histogram is large, emit separately if needed
        d["histogram"] = {str(k): v for k, v in self.histogram.items()}
        return d


def length_features(
    fragments: Iterable[Fragment],
    *,
    short: tuple[int, int] = (100, 150),
    long: tuple[int, int] = (151, 220),
    hist_range: tuple[int, int] = (1, 500),
) -> LengthFeatures:
    lengths = np.fromiter((f.length for f in fragments), dtype=np.int32)
    if lengths.size == 0:
        return LengthFeatures(0, 0.0, 0.0, 0, float("nan"), 0.0, 0.0, {})

    lo, hi = hist_range
    clipped = lengths[(lengths >= lo) & (lengths <= hi)]
    bincounts = np.bincount(clipped, minlength=hi + 1)
    mode = int(np.argmax(bincounts))

    n_short = int(((lengths >= short[0]) & (lengths <= short[1])).sum())
    n_long = int(((lengths >= long[0]) & (lengths <= long[1])).sum())
    slr = n_short / n_long if n_long else float("nan")

    frac_sub = float((lengths < 100).mean())
    frac_mono = float(((lengths >= 100) & (lengths <= 220)).mean())

    hist = {int(i): int(c) for i, c in enumerate(bincounts) if c > 0}
    return LengthFeatures(
        n_fragments=int(lengths.size),
        mean=float(lengths.mean()),
        median=float(np.median(lengths)),
        mode=mode,
        short_long_ratio=float(slr),
        frac_sub_nucleosomal=frac_sub,
        frac_mononucleosomal=frac_mono,
        histogram=hist,
    )


def binned_short_long_ratio(
    fragments: Iterable[Fragment],
    *,
    bin_size: int = 5_000_000,
    short: tuple[int, int] = (100, 150),
    long: tuple[int, int] = (151, 220),
) -> dict[tuple[str, int], float]:
    """DELFI-style genome-wide profile: short:long ratio per fixed-size bin.

    Returns a mapping ``{(chrom, bin_index): ratio}``. Downstream, these vectors
    are typically z-scored per sample and used as ML features for classification.
    """
    short_counts: dict[tuple[str, int], int] = {}
    long_counts: dict[tuple[str, int], int] = {}
    for f in fragments:
        key = (f.chrom, int(f.midpoint) // bin_size)
        if short[0] <= f.length <= short[1]:
            short_counts[key] = short_counts.get(key, 0) + 1
        elif long[0] <= f.length <= long[1]:
            long_counts[key] = long_counts.get(key, 0) + 1
    keys = set(short_counts) | set(long_counts)
    out = {}
    for k in sorted(keys):
        n_long = long_counts.get(k, 0)
        out[k] = short_counts.get(k, 0) / n_long if n_long else float("nan")
    return out
