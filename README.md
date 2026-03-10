# AAS Revisions

Analysis of the **Rate of Amino Acid Substitution (RAAS)** — quantifying how frequently Substituted Amino Acid Peptides (SAAPs) occur relative to their cognate Base Peptides (BPs) using DIA-MS data.

## Overview

Mistranslation during protein synthesis produces peptides with single amino acid substitutions (SAAPs). The **RAAS** is computed as the MS1 intensity of a SAAP divided by the intensity of its BP (the same tryptic peptide carrying the reference amino acid). A higher RAAS indicates a higher rate of substitution at that position.

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
