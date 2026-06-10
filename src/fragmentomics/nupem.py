"""NuPEM - Nucleosome-Phased End-Motif coupling.

Motivation
----------
Two cfDNA fragmentomic signals are usually computed and modelled independently:

1. **End motifs** capture the *sequence* preference of the nucleases that cut DNA
   (DNASE1L3, DFFB, DNASE1).
2. **Nucleosome positioning** (from WPS) captures the *structural* accessibility
   of chromatin - where the cutting can physically happen.

But these are mechanistically linked: a nuclease can only act where chromatin
lets it, so the *sequence* context of cut sites should depend on *where, relative
to the nucleosome, the cut occurred*. Prior work has noted "differentially phased
fragment-end signals", yet end-motif composition is typically pooled over all
fragments, discarding this phase dependence.

NuPEM makes that dependence an explicit, single-number feature:

* assign every fragment end a **phase** = signed distance to the nearest inferred
  dyad, folded onto the nucleosome repeat;
* partition ends into a **dyad-core** compartment and a **linker** compartment;
* compute the end-motif distribution within each compartment;
* the **NuPEM coupling score** is the Jensen-Shannon divergence between the two
  distributions - high coupling means the nuclease's sequence preference changes
  sharply with chromatin context.

This is presented as a *novel composite* of established primitives, not a claim of
a first-in-literature observation. Validate on labelled cohorts before drawing
biological conclusions.
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .fragments import Fragment, RefLookup
from .motifs import _VALID, all_kmers
from .utils import jensen_shannon_divergence, log2_enrichment, reverse_complement


@dataclass
class NuPEMResult:
    coupling_score: float                 # Jensen-Shannon divergence (bits), >= 0
    n_core_ends: int
    n_linker_ends: int
    core_top_enriched: list[tuple[str, float]] = field(default_factory=list)
    linker_top_enriched: list[tuple[str, float]] = field(default_factory=list)
    phase_profile: dict[int, float] = field(default_factory=dict)  # phase-bin -> MDS

    def to_dict(self) -> dict:
        return {
            "coupling_score": self.coupling_score,
            "n_core_ends": self.n_core_ends,
            "n_linker_ends": self.n_linker_ends,
            "core_top_enriched": self.core_top_enriched,
            "linker_top_enriched": self.linker_top_enriched,
            "phase_profile_mds": {str(k): v for k, v in self.phase_profile.items()},
        }


def _signed_phase(pos: float, dyads: np.ndarray, repeat: int) -> float | None:
    """Distance from ``pos`` to nearest dyad, folded into [-repeat/2, repeat/2]."""
    if dyads.size == 0:
        return None
    i = int(np.searchsorted(dyads, pos))
    candidates = []
    if i < dyads.size:
        candidates.append(dyads[i])
    if i > 0:
        candidates.append(dyads[i - 1])
    nearest = min(candidates, key=lambda d: abs(d - pos))
    d = pos - nearest
    # fold onto the repeat so positions one nucleosome away align in phase
    folded = ((d + repeat / 2) % repeat) - repeat / 2
    return folded


def _motif_at(ref: RefLookup, chrom: str, pos: int, strand: str, k: int) -> str | None:
    """5' end k-mer at an end position. strand '+' is a Watson 5' end at ``pos``;
    strand '-' is a Crick 5' end whose 5'->3' motif is the reverse complement."""
    if strand == "+":
        motif = ref(chrom, pos, pos + k)
    else:
        motif = ref(chrom, pos - k + 1, pos + 1)
        motif = reverse_complement(motif)
    if len(motif) == k and set(motif) <= _VALID:
        return motif
    return None


def nupem(
    fragments: Sequence[Fragment],
    dyads_by_chrom: dict[str, Sequence[int]],
    ref: RefLookup,
    *,
    k: int = 4,
    repeat: int = 190,
    core_halfwidth: int = 35,
    linker_window: tuple[int, int] = (60, 95),
    n_phase_bins: int = 12,
    top_n: int = 8,
) -> NuPEMResult:
    """Compute the NuPEM coupling score and supporting profiles.

    Parameters
    ----------
    dyads_by_chrom
        Inferred nucleosome dyad positions per chromosome (from
        :func:`fragmentomics.wps.call_nucleosomes`), assumed sortable.
    repeat
        Nucleosome repeat length used to fold phase (~185-200 bp in human).
    core_halfwidth
        Ends within +/- this distance of a dyad are "dyad-core".
    linker_window
        Absolute folded-distance band counted as "linker".
    """
    sorted_dyads = {c: np.asarray(sorted(d), dtype=float) for c, d in dyads_by_chrom.items()}

    core: Counter = Counter()
    linker: Counter = Counter()
    # phase-binned motif counts for the diversity profile
    bin_counts: dict[int, Counter] = {b: Counter() for b in range(n_phase_bins)}
    bin_width = repeat / n_phase_bins

    n_core = n_linker = 0
    for f in fragments:
        dyads = sorted_dyads.get(f.chrom)
        if dyads is None or dyads.size == 0:
            continue
        ends = ((f.start, "+"), (f.end - 1, "-"))
        for pos, strand in ends:
            phase = _signed_phase(pos, dyads, repeat)
            if phase is None:
                continue
            motif = _motif_at(ref, f.chrom, pos, strand, k)
            if motif is None:
                continue

            b = int((phase + repeat / 2) // bin_width)
            b = min(max(b, 0), n_phase_bins - 1)
            bin_counts[b][motif] += 1

            ad = abs(phase)
            if ad <= core_halfwidth:
                core[motif] += 1
                n_core += 1
            elif linker_window[0] <= ad <= linker_window[1]:
                linker[motif] += 1
                n_linker += 1

    coupling = jensen_shannon_divergence(dict(core), dict(linker)) if (core and linker) else 0.0

    core_enr = log2_enrichment(dict(core), dict(linker)) if (core and linker) else {}
    linker_enr = log2_enrichment(dict(linker), dict(core)) if (core and linker) else {}
    core_top = sorted(core_enr.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
    linker_top = sorted(linker_enr.items(), key=lambda kv: kv[1], reverse=True)[:top_n]

    from .utils import normalised_entropy

    phase_profile = {
        b: normalised_entropy([c.get(m, 0) for m in all_kmers(k)], n_categories=4 ** k)
        for b, c in bin_counts.items()
        if sum(c.values()) > 0
    }

    return NuPEMResult(
        coupling_score=round(float(coupling), 6),
        n_core_ends=n_core,
        n_linker_ends=n_linker,
        core_top_enriched=[(m, round(v, 4)) for m, v in core_top],
        linker_top_enriched=[(m, round(v, 4)) for m, v in linker_top],
        phase_profile={b: round(v, 5) for b, v in phase_profile.items()},
    )
