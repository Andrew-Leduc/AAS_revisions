# AAS Revisions

Revision analyses for **"Alternate RNA decoding results in stable and abundant proteins in mammals"**
Tsour, Machne, Leduc et al. — [bioRxiv 2024](https://doi.org/10.1101/2024.08.26.609665)

---

## Overview

Non-standard translation produces peptides carrying single amino acid substitutions (**SAAPs** — Substituted Amino Acid Peptides) alongside the canonical reference peptide (**BP** — Base Peptide). The **RAAS** (Rate of Amino Acid Substitution) is quantified as the MS1 intensity of a SAAP divided by that of its BP. A higher RAAS indicates more frequent substitution at that site.

The original study identified >60,000 fragmentation spectra representing 8,801 unique substitution sites across 1,782 genes from >1,000 human samples. Key findings include: hundreds of SAAPs are more abundant than their canonical BPs, substitution rates are shaped by protein stability and codon usage, and patterns show tissue- and cancer-specificity conserved between human and mouse.

These analyses extend the study by examining how protein degradation pathways affect SAAP accumulation in iPSC-derived neurons and mouse brain.

The analysis covers two datasets:

| Dataset | Notebook | Searched data |
|---------|----------|---------------|
| iNeuron (iPSC-derived neurons) | `Analysis/iNeuron.ipynb` | Google Drive — see below |
| Mouse Brain | `Analysis/MouseBrain.ipynb` | TBD |

---

## iNeuron Experiment

Cells were treated with inhibitors of the two main protein degradation pathways to probe how degradation affects SAAP accumulation:

| Condition | Inhibitor |
|-----------|-----------|
| DMSO | Vehicle control |
| Proteosome | Proteasome inhibitor |
| Lysosome | Lysosomal/autophagy inhibitor |
| Both | Proteasome + Lysosomal inhibitor |

Two differentiation timepoints (7-day, 14-day), 4 replicates per condition × timepoint = **32 raw files total**.

### Raw data path (Google Drive)
```
My Drive/MS/Users/aleduc/AAS_rev/report.parquet
```
Full local path:
```
/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/My Drive/MS/Users/aleduc/AAS_rev/report.parquet
```
DIA-NN output in `.parquet` format (~2.4M rows, 71 columns).

---

## Repository Structure

```
AAS_revisions/
├── Analysis/
│   ├── iNeuron.ipynb       # iNeuron RAAS analysis
│   └── MouseBrain.ipynb    # Mouse brain RAAS analysis (in progress)
└── meta_files/
    ├── All_SAAP.fasta                        # SAAP search library (8,768 entries)
    ├── output_MTP_cl_fp.fasta                # Mouse brain SAAP library
    ├── Supplemental_Data_2.SAAP_proteins.xlsx # SAAP → BP peptide mapping (10,180 pairs)
    ├── neuron_meta.csv                       # Run-level metadata (condition, timepoint)
    ├── high_quality_SAAPs.xlsx               # Curated high-confidence SAAP list
    └── SILAC.xlsx                            # SILAC-derived degradation rates per peptide
```

---

## Analysis Pipeline (iNeuron)

1. **Load** DIA-NN report and SAAP→BP mapping
2. **Filter** to SAAP and BP peptides; sum MS1/MS2 intensities per run
3. **Compute RAAS ratios** (SAAP MS1 / BP MS1, SAAP MS2 / BP MS2) per run
4. **QC plots**: unique SAAPs per condition, unique total peptides per condition, MS1 vs MS2 ratio scatter
5. **Filter** to high-confidence SAAPs: peptides where SAAP degrades ≥2× more slowly than BP (from SILAC data) — stabilised SAAPs are expected to accumulate
6. **RAAS per condition**: beeswarm plot of mean log10(RAAS) per SAAP per condition with median overlay
7. **Normalised RAAS**: per-SAAP mean subtracted across conditions to highlight condition-driven shifts
