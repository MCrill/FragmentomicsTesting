"""Fragment data model and I/O.

A *fragment* is the inferred original cfDNA molecule: the genomic interval from
the leftmost mapped base of a properly-paired read pair to the rightmost. We work
in 0-based, half-open coordinates throughout (BED convention).

The BAM reader is the only component that depends on ``pysam``; the rest of the
library operates on plain :class:`Fragment` objects so it can be tested and reused
without alignment files.
"""
from __future__ import annotations

import gzip
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

# A reference lookup returns the uppercase reference sequence for a 0-based,
# half-open interval [start, end) on ``chrom``. Production code backs this with a
# pysam.FastaFile; tests back it with an in-memory dict.
RefLookup = Callable[[str, int, int], str]


@dataclass(frozen=True, slots=True)
class Fragment:
    chrom: str
    start: int  # 0-based, inclusive (Watson 5' end)
    end: int    # 0-based, exclusive (Crick 5' end)

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


def read_fragments_from_bam(
    bam_path: str | Path,
    *,
    min_mapq: int = 30,
    min_length: int = 50,
    max_length: int = 1000,
    region: str | None = None,
    exclude_flags: int = 0xF04,  # unmapped, secondary, qcfail, dup, supplementary
) -> Iterator[Fragment]:
    """Yield :class:`Fragment` objects from a coordinate-sorted, indexed BAM.

    Only the forward read of each properly-paired pair with a positive template
    length (TLEN) is used, which yields each fragment exactly once. ``region`` is
    an htslib region string (e.g. ``"chr1:1-1000000"``).
    """
    import pysam  # local import keeps pysam optional for pure-Python users

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(region=region):
            if read.flag & exclude_flags:
                continue
            if not read.is_proper_pair:
                continue
            if read.mapping_quality < min_mapq:
                continue
            tlen = read.template_length
            if tlen <= 0:  # take only the leftmost mate to avoid double counting
                continue
            length = tlen
            if length < min_length or length > max_length:
                continue
            start = read.reference_start
            yield Fragment(read.reference_name, start, start + length)


def read_fragments_from_bedpe(path: str | Path) -> Iterator[Fragment]:
    """Read fragments from a (optionally bgzipped) 3+ column BED of fragments."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            yield Fragment(chrom, int(start), int(end))


def write_bed(fragments: Iterable[Fragment], path: str | Path) -> int:
    """Write fragments as BED3+length. Returns the number of records written."""
    n = 0
    with open(path, "w") as fh:
        for f in fragments:
            fh.write(f"{f.chrom}\t{f.start}\t{f.end}\t{f.length}\n")
            n += 1
    return n


def dict_ref_lookup(reference: dict[str, str]) -> RefLookup:
    """Build a :data:`RefLookup` from an in-memory ``{chrom: sequence}`` dict.

    Out-of-bounds requests return Ns so motif extraction degrades gracefully at
    contig edges rather than raising.
    """
    def _lookup(chrom: str, start: int, end: int) -> str:
        seq = reference.get(chrom, "")
        if start < 0:
            pad = "N" * (-start)
            return (pad + seq[0:end]).upper()[: end - start]
        out = seq[start:end].upper()
        if len(out) < end - start:
            out = out + "N" * (end - start - len(out))
        return out

    return _lookup


def fasta_ref_lookup(fasta_path: str | Path) -> RefLookup:
    """Build a :data:`RefLookup` backed by an indexed FASTA via pysam."""
    import pysam

    fasta = pysam.FastaFile(str(fasta_path))

    def _lookup(chrom: str, start: int, end: int) -> str:
        s = max(start, 0)
        try:
            seq = fasta.fetch(chrom, s, end).upper()
        except (KeyError, ValueError):
            return "N" * (end - start)
        if start < 0:
            seq = "N" * (-start) + seq
        if len(seq) < end - start:
            seq = seq + "N" * (end - start - len(seq))
        return seq

    return _lookup
