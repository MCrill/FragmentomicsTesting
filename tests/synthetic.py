"""Generate synthetic cfDNA fragments with *known* structure for testing.

Biologically faithful in the ways that matter for these tests:

* Mononucleosomal fragments (~167 bp) are centred *on* nucleosome dyads, so they
  protect the dyad and their ends fall in the flanking linkers - this is what lets
  L-WPS recover dyad positions.
* A sub-nucleosomal population (~50 bp), also centred on dyads, places ends close
  to the dyad - this populates the NuPEM "dyad-core" compartment.

With ``phase_bias=True`` the reference is edited so dyad-proximal ends carry the
motif ``CGCG`` and linker ends carry ``TTTT``, creating a detectable coupling
signal; with ``phase_bias=False`` no such structure exists and coupling ~ 0.
"""
from __future__ import annotations

import random

from fragmentomics.fragments import Fragment


def make_reference(chrom: str = "chrT", length: int = 20_000, seed: int = 7) -> dict[str, str]:
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    return {chrom: seq}


def _stamp(seq_list: list[str], pos: int, motif: str) -> None:
    for i, base in enumerate(motif):
        if 0 <= pos + i < len(seq_list):
            seq_list[pos + i] = base


def simulate_fragments(
    reference: dict[str, str],
    *,
    chrom: str = "chrT",
    n_fragments: int = 4000,
    repeat: int = 190,
    mono_len_mean: int = 167,
    mono_len_sd: int = 12,
    short_len_mean: int = 50,
    short_len_sd: int = 6,
    short_fraction: float = 0.35,
    phase_bias: bool = True,
    seed: int = 11,
):
    """Return (fragments, possibly-edited reference, true dyad positions)."""
    rng = random.Random(seed)
    seq = list(reference[chrom])
    L = len(seq)
    dyads = list(range(repeat, L - repeat, repeat))
    fragments: list[Fragment] = []

    def stamp_end(pos: int) -> None:
        nearest = min(dyads, key=lambda x: abs(x - pos))
        d = abs(pos - nearest)
        if d <= 35:
            _stamp(seq, pos, "CGCG")
        elif 60 <= d <= 95:
            _stamp(seq, pos, "TTTT")

    for _ in range(n_fragments):
        d = rng.choice(dyads)
        center = d + rng.randint(-12, 12)
        if rng.random() < short_fraction:
            length = max(35, int(rng.gauss(short_len_mean, short_len_sd)))
        else:
            length = max(120, int(rng.gauss(mono_len_mean, mono_len_sd)))
        start = center - length // 2
        end = start + length
        if start < 5 or end > L - 5:
            continue
        fragments.append(Fragment(chrom, start, end))
        if phase_bias:
            stamp_end(start)
            stamp_end(end - 4)

    return fragments, {chrom: "".join(seq)}, dyads
