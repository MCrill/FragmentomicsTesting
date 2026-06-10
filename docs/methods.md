# Methods

This document defines each feature precisely so results are reproducible and the
maths is auditable. Coordinates are 0-based, half-open (BED convention) throughout.

## 1. Fragment definition

A *fragment* is the inferred original cfDNA molecule. From a coordinate-sorted,
duplicate-marked, properly-paired BAM we take the leftmost read of each pair with a
positive template length (TLEN):

```
fragment = (chrom, reference_start, reference_start + TLEN)
length   = TLEN
```

Using only the positive-TLEN mate counts each molecule exactly once.

## 2. Fragment-length features

* Summary stats: mean, median, mode of the length histogram.
* **Short:long ratio** (DELFI; Cristiano et al., *Nature* 2019): with short = 100-150 bp
  and long = 151-220 bp, the per-bin ratio across fixed genomic bins forms a
  genome-wide fragmentation profile. Tumour-derived cfDNA is shifted toward shorter
  fragments, raising this ratio in affected bins.

## 3. End motifs and the Motif Diversity Score (MDS)

For each fragment we read the k-mer (default 4) immediately 5' of each strand's
cleavage site **from the reference**:

* Watson 5' end: `ref[start : start+k]`
* Crick 5' end: `revcomp(ref[end-k : end])`

The **MDS** (Jiang et al., *PNAS* 2020) is the Shannon entropy of the 4^k motif
distribution normalised to [0, 1]:

```
MDS = H(motif_freqs) / log2(4^k)
```

A drop in MDS reflects more uniform cleavage (e.g. loss of DNASE1L3 substrate
preference) and is a recurrent signal in cancer plasma.

## 4. Windowed Protection Score (WPS) and nucleosome calling

L-WPS (Snyder et al., *Cell* 2016). For a window of width *w* (default 120 bp) and
the long fragment fraction (120-180 bp), at each position *p*:

```
WPS(p) = (# fragments fully spanning the window centred at p)
       - (# fragments with an endpoint inside that window)
```

Implemented with cumulative-sum arrays (O(region + fragments)). The smoothed track
is peak-called (min spacing ~150 bp) to yield nucleosome **dyad** positions.

## 5. NuPEM - Nucleosome-Phased End-Motif coupling (novel composite)

**Idea.** End motifs encode *nuclease sequence preference*; nucleosome positioning
encodes *structural accessibility*. These are mechanistically linked but are almost
always modelled independently, pooling end motifs over all fragments and discarding
the dependence on where, relative to the nucleosome, a cut occurred. NuPEM makes
that dependence a single interpretable number.

**Procedure.**

1. Call nucleosome dyads from L-WPS (§4).
2. For every fragment end, compute its **phase** = signed distance to the nearest
   dyad, folded onto the nucleosome repeat length (default 190 bp).
3. Partition ends into two compartments:
   * **dyad-core**: `|phase| <= 35 bp`
   * **linker**: `60 <= |phase| <= 95 bp`
4. Compute the end-motif distribution within each compartment.
5. **NuPEM coupling score** = Jensen-Shannon divergence between the dyad-core and
   linker motif distributions.

Also reported: per-compartment log2-enriched motifs, and a phase-binned MDS profile
(motif diversity as a function of position relative to the dyad).

**Interpretation.** Coupling ~ 0 means cleavage sequence preference is independent
of chromatin context; high coupling means the nuclease "sees" a different sequence
landscape at dyads vs linkers. The hypothesis worth testing on real cohorts is that
*the strength and shape of this coupling differs between healthy and disease cfDNA*,
because disease shifts both the nuclease repertoire and chromatin organisation.

**Honesty note.** NuPEM is a novel *composition* of well-established primitives
(end motifs + WPS nucleosome maps). Phase-dependent fragment-end signals have been
described qualitatively in the literature; the explicit JSD-based coupling score and
the phase-binned MDS profile here are, to our knowledge, not a standard named
method, but no claim of biological priority is made. Establish discriminative value
on labelled cohorts with proper cross-validation and batch controls before drawing
conclusions.

## References

* Snyder MW et al. Cell-free DNA comprises an in vivo nucleosome footprint. *Cell* 2016.
* Cristiano S et al. Genome-wide cell-free DNA fragmentation in patients with cancer (DELFI). *Nature* 2019.
* Jiang P et al. Plasma DNA end-motif profiling (MDS). *PNAS* / *Cancer Discovery* 2020.
* Han DSC et al. The biology of cell-free DNA fragmentation (DNASE1, DNASE1L3, DFFB). *AJHG* 2020.
* Lo YMD, Han DSC, Jiang P, Chiu RWK. Epigenetics, fragmentomics, and topology of cfDNA. *Science* 2021.

*Verify exact citations/years against the primary sources before publication.*
