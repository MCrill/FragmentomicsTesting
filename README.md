# cfDNA Fragmentomics + NuPEM

A reproducible, end-to-end **cell-free DNA fragmentomics** pipeline that goes from
raw NGS FASTQs to a panel of fragmentomic features — and introduces a novel composite
metric, **NuPEM** (*Nucleosome-Phased End-Motif coupling*), which links nuclease
sequence preference to chromatin structure in a single interpretable number.

[![CI](https://github.com/morgan/cfdna-fragmentomics-nupem/actions/workflows/ci.yml/badge.svg)](https://github.com/morgan/cfdna-fragmentomics-nupem/actions)
![python](https://img.shields.io/badge/python-3.10%2B-blue)
![nextflow](https://img.shields.io/badge/nextflow-DSL2-23aa62)
![license](https://img.shields.io/badge/license-MIT-green)

---

## Why this exists

cfDNA is fragmented non-randomly: the nucleases that cut it (DNASE1L3, DFFB,
DNASE1) leave characteristic **end motifs**, and the **nucleosomes** they cut around
leave characteristic protection footprints. These two signals are biologically
coupled — a nuclease can only cut where chromatin permits — yet they are almost
always measured *separately*, with end motifs pooled across all fragments.

This project computes the standard fragmentomic feature set **and** asks a question
that the pooled view throws away:

> Does the sequence context of cfDNA cut sites depend on *where, relative to the
> nucleosome,* the cut happened — and does that coupling carry signal?

That question is operationalised as the **NuPEM coupling score**. See
[`docs/methods.md`](docs/methods.md) for the precise definition. NuPEM is a novel
*composition* of established primitives, not a claim of biological priority — it is
designed to be a clean, testable feature you can validate on labelled cohorts.

## The feature panel

| Feature | What it captures | Reference |
|---|---|---|
| Fragment-length distribution & **short:long ratio** | Sub- vs mono-nucleosomal balance; tumour size shift | DELFI, Cristiano 2019 |
| **End motifs** + **Motif Diversity Score (MDS)** | Nuclease sequence preference / cleavage uniformity | Jiang 2020 |
| **L-WPS** + nucleosome dyad calls | Chromatin protection footprint | Snyder 2016 |
| **NuPEM coupling** (novel) | Dependence of end-motif composition on nucleosome phase | this work |

## Pipeline

```
FASTQ (PE)
   │  fastp            adapter/quality trim + QC
   ▼
   │  bwa-mem2         align to GRCh38
   ▼
   │  samtools         sort, markdup
   ▼
   │  samtools+bedtools  filter: proper pairs, MAPQ≥30, autosomes, −ENCODE blacklist
   ▼
   │  fragmentomics extract   BAM → bgzipped+tabix fragment BED
   ▼
   │  fragmentomics features  length • motifs/MDS • WPS/nucleosomes • NuPEM → JSON
   ▼
   └─ (optional) fragmentomics.modeling   feature matrix → cross-validated classifier
```

Orchestrated in **Nextflow DSL2** with `docker` / `singularity` / `conda` / `slurm`
profiles, per-process resource labels, and MultiQC aggregation.

## Quickstart

### Install the analysis library

```bash
pip install -e ".[all,dev]"     # numpy, click, pysam, scipy, pytest, ruff
pytest -q                       # run the test suite (uses synthetic data, no genome needed)
```

### Run on a sample's BAM (library only)

```bash
# 1. extract fragments
fragmentomics extract --bam sample.bam --out sample.fragments.bed
sort -k1,1 -k2,2n sample.fragments.bed | bgzip > sample.fragments.bed.gz && tabix -p bed sample.fragments.bed.gz

# 2. compute the full feature panel (incl. NuPEM)
fragmentomics features \
    --fragments sample.fragments.bed.gz \
    --reference GRCh38.fa \
    --regions assets/regions.example.bed \
    --sample-id sample --out sample.features.json
```

### Run the whole pipeline from FASTQ

```bash
nextflow run main.nf -profile docker \
    --samplesheet samplesheet.csv \
    --reference GRCh38.fa \
    --bwa_index 'index/GRCh38.fa.*' \
    --regions_bed assets/regions.example.bed \
    --blacklist ENCODE_blacklist.bed \
    --outdir results
```

### Train a classifier on the features

```bash
python -m fragmentomics.modeling \
    --features 'results/features/*.features.json' \
    --labels labels.csv \
    --out model_report.json
```

## Using NuPEM programmatically

```python
from fragmentomics import (
    Region, read_fragments_from_bedpe, fasta_ref_lookup, call_nucleosomes, wps_track, nupem,
)

frags  = list(read_fragments_from_bedpe("sample.fragments.bed.gz"))
ref    = fasta_ref_lookup("GRCh38.fa")
region = Region("chr12", 6_533_927, 6_543_927)

dyads  = {"chr12": call_nucleosomes(wps_track(frags, region), region)}
result = nupem(frags, dyads, ref, k=4)

print(result.coupling_score)        # Jensen–Shannon divergence, core vs linker
print(result.core_top_enriched)     # motifs enriched at dyad-proximal cut sites
```

## Validated behaviour

The metric is unit-tested against synthetic cfDNA with *known* nucleosome positions
and an injected phase→motif bias. On that ground truth:

* L-WPS recovers 104/104 true dyad positions;
* NuPEM coupling rises from **0.04** (no phase structure) to **~1.0** (strong
  phase-dependent motif bias) and recovers the injected core/linker motifs.

See `tests/test_nupem.py`.

## Repository layout

```
src/fragmentomics/      # pure-Python feature library (the core)
  ├─ fragments.py       #   data model + BAM/BEDPE/FASTA I/O
  ├─ lengths.py         #   length distribution + DELFI short:long ratio
  ├─ motifs.py          #   end motifs + MDS
  ├─ wps.py             #   L-WPS + nucleosome calling
  ├─ nupem.py           #   ★ the novel coupling metric
  ├─ features.py        #   per-sample orchestrator
  ├─ modeling.py        #   feature matrix + classifier
  └─ cli.py             #   `fragmentomics` CLI
main.nf, modules/       # Nextflow DSL2 pipeline
tests/                  # pytest suite + synthetic data generator
docs/methods.md         # precise feature definitions + references
Dockerfile, environment.yml, .github/workflows/ci.yml
```

## Design choices worth noting

* **I/O is isolated.** Only `fragments.py` touches `pysam`; every feature function
  operates on plain `Fragment` objects and a `RefLookup` callable, so the science is
  testable without alignment files and reusable outside the pipeline.
* **WPS is vectorised** with cumulative sums — O(region + fragments), not O(region ×
  fragments).
* **Reproducibility is first-class** — pinned conda env, container, CI across Python
  3.10–3.12, and Nextflow execution reports.

## Caveats & honest scope

This is a research/portfolio pipeline, **not** a clinical assay. The NuPEM score is a
hypothesis-generating feature: demonstrate discriminative value on labelled cohorts
with proper cross-validation, batch-effect control, and GC correction before drawing
any biological or clinical conclusion. Citations in `docs/methods.md` should be
verified against primary sources.

## License

MIT — see [`LICENSE`](LICENSE).
