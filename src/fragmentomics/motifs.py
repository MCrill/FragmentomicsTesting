"""Fragment end-motif features.

For each fragment we read the k-mer (default 4-mer) immediately 5' of each strand's
cleavage site from the *reference* genome:

* Watson 5' end: ``ref[start : start + k]``
* Crick 5' end:  ``reverse_complement(ref[end - k : end])``

End-motif frequencies reflect the sequence preference of the nucleases that
generate cfDNA (DFFB, DNASE1L3, DNASE1). The Motif Diversity Score (MDS; Jiang
et al. 2020) is the Shannon entropy of the 4^k motif distribution normalised to
[0, 1]; a drop in MDS (more uniform cleavage) is a recurrent cancer signal.
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from itertools import product

from .fragments import Fragment, RefLookup
from .utils import normalised_entropy, reverse_complement

_VALID = set("ACGT")


def all_kmers(k: int) -> list[str]:
    return ["".join(p) for p in product("ACGT", repeat=k)]


def end_motif_counts(
    fragments: Iterable[Fragment],
    ref: RefLookup,
    *,
    k: int = 4,
) -> Counter:
    """Count both-strand 5' end k-mers across all fragments.

    Motifs containing any non-ACGT base (e.g. spanning an N-masked region) are
    skipped so they do not distort the diversity score.
    """
    counts: Counter = Counter()
    for f in fragments:
        watson = ref(f.chrom, f.start, f.start + k)
        crick = reverse_complement(ref(f.chrom, f.end - k, f.end))
        for motif in (watson, crick):
            if len(motif) == k and set(motif) <= _VALID:
                counts[motif] += 1
    return counts


@dataclass
class MotifFeatures:
    k: int
    n_ends: int
    mds: float
    frequencies: dict[str, float]
    top_motifs: list[tuple[str, float]]

    def to_dict(self) -> dict:
        return {
            "k": self.k,
            "n_ends": self.n_ends,
            "mds": self.mds,
            "top_motifs": self.top_motifs,
        }


def motif_features(
    fragments: Iterable[Fragment],
    ref: RefLookup,
    *,
    k: int = 4,
    top_n: int = 10,
) -> MotifFeatures:
    counts = end_motif_counts(fragments, ref, k=k)
    total = sum(counts.values())
    n_categories = 4 ** k
    mds = normalised_entropy(
        [counts.get(m, 0) for m in all_kmers(k)], n_categories=n_categories
    )
    freqs = {m: counts.get(m, 0) / total for m in all_kmers(k)} if total else {}
    top = sorted(freqs.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
    return MotifFeatures(
        k=k,
        n_ends=total,
        mds=mds,
        frequencies=freqs,
        top_motifs=[(m, round(v, 5)) for m, v in top],
    )
