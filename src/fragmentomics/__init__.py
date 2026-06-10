"""cfDNA fragmentomics toolkit.

Public API re-exports the most commonly used building blocks so users can
``from fragmentomics import Fragment, compute_sample_features`` directly.
"""
from __future__ import annotations

__version__ = "0.1.0"

from .features import SampleFeatures, compute_sample_features
from .fragments import (
    Fragment,
    dict_ref_lookup,
    fasta_ref_lookup,
    read_fragments_from_bam,
    read_fragments_from_bedpe,
)
from .lengths import LengthFeatures, binned_short_long_ratio, length_features
from .motifs import MotifFeatures, end_motif_counts, motif_features
from .nupem import NuPEMResult, nupem
from .wps import Region, call_nucleosomes, wps_track

__all__ = [
    "__version__",
    "Fragment",
    "Region",
    "SampleFeatures",
    "LengthFeatures",
    "MotifFeatures",
    "NuPEMResult",
    "compute_sample_features",
    "length_features",
    "binned_short_long_ratio",
    "motif_features",
    "end_motif_counts",
    "wps_track",
    "call_nucleosomes",
    "nupem",
    "read_fragments_from_bam",
    "read_fragments_from_bedpe",
    "dict_ref_lookup",
    "fasta_ref_lookup",
]
