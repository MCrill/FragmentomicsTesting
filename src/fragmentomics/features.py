"""Top-level orchestration: turn a set of fragments + reference into a single
per-sample feature record combining all modules (length, motif, WPS-derived
nucleosome map, and the NuPEM coupling score).
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from .fragments import Fragment, RefLookup
from .lengths import LengthFeatures, length_features
from .motifs import MotifFeatures, motif_features
from .nupem import NuPEMResult, nupem
from .wps import Region, call_nucleosomes, wps_track


@dataclass
class SampleFeatures:
    sample_id: str
    n_fragments: int
    lengths: LengthFeatures
    motifs: MotifFeatures
    nupem: NuPEMResult
    n_nucleosomes: int
    regions_analysed: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "n_fragments": self.n_fragments,
            "n_nucleosomes": self.n_nucleosomes,
            "regions_analysed": self.regions_analysed,
            "lengths": self.lengths.to_dict(),
            "motifs": self.motifs.to_dict(),
            "nupem": self.nupem.to_dict(),
        }


def compute_sample_features(
    sample_id: str,
    fragments: Sequence[Fragment],
    ref: RefLookup,
    regions: Sequence[Region],
    *,
    k: int = 4,
    wps_window: int = 120,
) -> SampleFeatures:
    """Compute the full feature panel for one sample.

    Nucleosome dyads are called per analysis region from the L-WPS track, then
    pooled to drive the NuPEM phase assignment. ``regions`` should be a curated
    set of mappable, blacklist-free windows (e.g. around housekeeping-gene TSSs or
    tiled open-chromatin regions) - passing the whole genome is intentionally not
    supported here to keep memory bounded.
    """
    fragments = list(fragments)

    dyads_by_chrom: dict[str, list[int]] = {}
    for region in regions:
        track = wps_track(fragments, region, window=wps_window)
        dyads = call_nucleosomes(track, region)
        dyads_by_chrom.setdefault(region.chrom, []).extend(dyads)
    for c in dyads_by_chrom:
        dyads_by_chrom[c] = sorted(set(dyads_by_chrom[c]))

    lengths = length_features(fragments)
    motifs = motif_features(fragments, ref, k=k)
    nupem_res = nupem(fragments, dyads_by_chrom, ref, k=k)

    return SampleFeatures(
        sample_id=sample_id,
        n_fragments=len(fragments),
        lengths=lengths,
        motifs=motifs,
        nupem=nupem_res,
        n_nucleosomes=sum(len(v) for v in dyads_by_chrom.values()),
        regions_analysed=[f"{r.chrom}:{r.start}-{r.end}" for r in regions],
    )
