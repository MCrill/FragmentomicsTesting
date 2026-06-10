"""Command-line interface for the fragmentomics toolkit.

Subcommands::

    fragmentomics extract   BAM  -> fragment BED
    fragmentomics features  fragments + reference + regions -> JSON feature record

The CLI is deliberately thin: it wires I/O to the pure-Python feature functions so
the same logic is exercised by the unit tests.
"""
from __future__ import annotations

import json
from pathlib import Path

import click

from . import __version__
from .features import compute_sample_features
from .fragments import (
    Fragment,
    fasta_ref_lookup,
    read_fragments_from_bam,
    read_fragments_from_bedpe,
    write_bed,
)
from .wps import Region


@click.group()
@click.version_option(__version__)
def main() -> None:
    """cfDNA fragmentomics feature extraction (length, end motifs, WPS, NuPEM)."""


@main.command()
@click.option("--bam", required=True, type=click.Path(exists=True), help="Indexed, sorted BAM.")
@click.option("--out", required=True, type=click.Path(), help="Output fragment BED.")
@click.option("--min-mapq", default=30, show_default=True)
@click.option("--min-length", default=50, show_default=True)
@click.option("--max-length", default=1000, show_default=True)
@click.option("--region", default=None, help="htslib region string, e.g. chr1:1-1000000")
def extract(bam: str, out: str, min_mapq: int, min_length: int, max_length: int, region: str | None) -> None:
    """Extract cfDNA fragments from a BAM into a BED file."""
    frags = read_fragments_from_bam(
        bam, min_mapq=min_mapq, min_length=min_length, max_length=max_length, region=region
    )
    n = write_bed(frags, out)
    click.echo(f"Wrote {n} fragments to {out}", err=True)


def _parse_regions(path: str) -> list[Region]:
    regions: list[Region] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, *_ = line.split("\t")
            regions.append(Region(chrom, int(start), int(end)))
    return regions


@main.command()
@click.option("--fragments", "frag_path", required=True, type=click.Path(exists=True))
@click.option("--reference", required=True, type=click.Path(exists=True), help="Indexed FASTA (.fai).")
@click.option("--regions", required=True, type=click.Path(exists=True), help="BED of analysis regions.")
@click.option("--sample-id", required=True)
@click.option("--k", default=4, show_default=True, help="End-motif k-mer size.")
@click.option("--out", default="-", help="Output JSON path ('-' for stdout).")
def features(frag_path: str, reference: str, regions: str, sample_id: str, k: int, out: str) -> None:
    """Compute the full fragmentomic feature panel for one sample."""
    frags: list[Fragment] = list(read_fragments_from_bedpe(frag_path))
    ref = fasta_ref_lookup(reference)
    region_list = _parse_regions(regions)
    result = compute_sample_features(sample_id, frags, ref, region_list, k=k)
    payload = json.dumps(result.to_dict(), indent=2)
    if out == "-":
        click.echo(payload)
    else:
        Path(out).write_text(payload)
        click.echo(f"Wrote features for {sample_id} to {out}", err=True)


if __name__ == "__main__":  # pragma: no cover
    main()
