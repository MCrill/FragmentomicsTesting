"""Windowed Protection Score (WPS) and nucleosome dyad calling.

WPS (Snyder et al. 2016) summarises, at each genomic position, the balance between
fragments that *protect* (fully span) a window and fragments whose *endpoint* falls
inside it. Peaks in the smoothed L-WPS track mark the centres (dyads) of
well-positioned nucleosomes.

The default long-fraction parameters (120 bp window, 120-180 bp fragments) follow
the original L-WPS definition.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np

from .fragments import Fragment


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def wps_track(
    fragments: Iterable[Fragment],
    region: Region,
    *,
    window: int = 120,
    min_len: int = 120,
    max_len: int = 180,
) -> np.ndarray:
    """Per-base L-WPS over ``region`` (length == region.length).

    For each position ``p`` and window ``[p - w/2, p + w/2)``:
        WPS(p) = (#fragments fully spanning the window)
                 - (#fragments with an endpoint inside the window)

    Implemented with cumulative-sum arrays so the cost is O(region + fragments)
    rather than O(region * fragments).
    """
    half = window // 2
    L = region.length
    # span_delta accumulates +1 over positions where a fragment fully spans the
    # window; endpoint_delta accumulates +1 where an endpoint lies in the window.
    span = np.zeros(L + 1, dtype=np.int32)
    endp = np.zeros(L + 1, dtype=np.int32)

    for f in fragments:
        if f.chrom != region.chrom or not (min_len <= f.length <= max_len):
            continue
        s = f.start - region.start
        e = f.end - region.start  # exclusive

        # A window centred at p is fully spanned when s <= p-half and e-1 >= p+half-1
        # => p in [s + half, e - half). Clamp to region.
        span_lo = max(s + half, 0)
        span_hi = min(e - half, L)
        if span_hi > span_lo:
            span[span_lo] += 1
            span[span_hi] -= 1

        # An endpoint at position x contributes to windows centred at p where
        # p - half <= x < p + half  =>  p in (x - half, x + half].
        for x in (s, e - 1):
            ep_lo = max(x - half + 1, 0)
            ep_hi = min(x + half + 1, L)
            if ep_hi > ep_lo:
                endp[ep_lo] += 1
                endp[ep_hi] -= 1

    span_cov = np.cumsum(span[:-1])
    endp_cov = np.cumsum(endp[:-1])
    return (span_cov - endp_cov).astype(np.float64)


def _smooth(track: np.ndarray, window: int = 21) -> np.ndarray:
    """Light moving-average smoothing (odd window) to stabilise peak calling."""
    if window <= 1 or track.size < window:
        return track
    if window % 2 == 0:
        window += 1
    kernel = np.ones(window) / window
    return np.convolve(track, kernel, mode="same")


def call_nucleosomes(
    track: np.ndarray,
    region: Region,
    *,
    min_distance: int = 150,
    smooth_window: int = 21,
    height_quantile: float = 0.5,
) -> list[int]:
    """Return absolute genomic positions of inferred nucleosome dyads.

    Peaks are local maxima of the smoothed track, separated by at least
    ``min_distance`` bp (~ nucleosome repeat length) and above the
    ``height_quantile`` of the track. Uses scipy if available, else a NumPy
    fallback, so the module has no hard scipy requirement at import time.
    """
    sm = _smooth(track, smooth_window)
    if sm.size == 0:
        return []
    height = float(np.quantile(sm, height_quantile))

    try:
        from scipy.signal import find_peaks

        idx, _ = find_peaks(sm, distance=min_distance, height=height)
        peaks = list(idx)
    except Exception:  # pragma: no cover - fallback path
        peaks = _find_peaks_numpy(sm, min_distance, height)

    return [region.start + int(i) for i in peaks]


def _find_peaks_numpy(sm: np.ndarray, min_distance: int, height: float) -> list[int]:
    candidates = np.where((sm[1:-1] > sm[:-2]) & (sm[1:-1] >= sm[2:]) & (sm[1:-1] >= height))[0] + 1
    selected: list[int] = []
    for c in sorted(candidates, key=lambda i: sm[i], reverse=True):
        if all(abs(c - s) >= min_distance for s in selected):
            selected.append(int(c))
    return sorted(selected)
