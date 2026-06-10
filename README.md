# cfDNA Fragmentomics + Tool

Reproducible end-to-end cell-free DNA fragmentomics pipeline that goes from raw NGS FASTQs to a panel of fragmentomic features. The new feature, NuPEM (Nucleosome-Phased End-Motif coupling), links nuclease sequence preference to chromatin structure in a single interpretable number.

## Why this exists

cfDNA is fragmented non-randomly: the nucleases that cut it (DNASE1L3, DFFB, DNASE1) leave characteristic end motifs, and the nucleosomes they cut around leave characteristic protection footprints. These two signals are biologically coupled in that a nuclease can only cut where chromatin permits. Yet, they are almost always measured separately, with end motifs pooled across all fragments. 

This project computes the standard fragmentomic feature set and asks a question that the pooled view throws away - 

Does the sequence context of cfDNA cut sites depend on: where relative to the nucleosome the cut happened, and does that coupling carry a signal?

This question can be answered via the NuPEM coupling score. 

## The Feature Panel

