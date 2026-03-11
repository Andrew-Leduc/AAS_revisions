# CLAUDE.md — Guide for AI-assisted analysis

This file tells Claude (or any AI assistant) how this repository is structured, what the analyses do, and how to work with the code safely.

---

## Project summary

Revision analyses for **"Alternate RNA decoding results in stable and abundant proteins in mammals"**
Tsour, Machne, Leduc et al. — [bioRxiv 2024](https://doi.org/10.1101/2024.08.26.609665)

We quantify **RAAS** (Rate of Amino Acid Substitution) — how frequently a given amino acid position is mistranslated — across two datasets:

| Dataset | Notebook | Quantification |
|---------|----------|----------------|
| iPSC-derived neurons (iNeuron) | `Analysis/iNeuron.ipynb` | DIA-NN label-free (MS1 intensity) |
| Mouse brain young vs old | `Analysis/MouseBrain.ipynb` | FragPipe TMT10 |

---

## Key terminology

| Term | Meaning |
|------|---------|
| **SAAP** | Substituted Amino Acid Peptide — a peptide carrying a single amino acid substitution relative to the canonical sequence |
| **BP** | Base Peptide — the canonical reference peptide that a SAAP is derived from |
| **RAAS** | Rate of Amino Acid Substitution = abundance(SAAP) / abundance(BP). Higher RAAS = more frequent mistranslation at that site |

---

## Repository layout

```
AAS_revisions/
├── Analysis/
│   ├── iNeuron.ipynb          # iNeuron RAAS analysis (DIA-NN / label-free)
│   └── MouseBrain.ipynb       # Mouse brain RAAS analysis (TMT10)
├── meta_files/
│   ├── All_SAAP.fasta                         # iNeuron SAAP search library (8,768 entries)
│   ├── output_MTP_cl_fp.fasta                 # Mouse brain SAAP search library (MTP_ + new_peptide entries)
│   ├── Supplemental_Data_2.SAAP_proteins.xlsx # iNeuron SAAP → BP mapping (10,180 pairs)
│   ├── neuron_meta.csv                        # iNeuron run metadata (Raw file, Condition, Day)
│   ├── high_quality_SAAPs.xlsx                # Curated high-confidence SAAP list (5,999 entries)
│   └── SILAC.xlsx                             # SILAC degradation rates per SAAP and BP
├── README.md
└── CLAUDE.md                                  # this file
```

All paths in notebooks use `pathlib.Path` with `META_DIR = Path('..') / 'meta_files'` so they work regardless of where the repo is cloned. Raw mass spec data lives on Google Drive (paths below).

---

## Raw data paths (Google Drive)

### iNeuron
```
My Drive/MS/Users/aleduc/AAS_rev/report.parquet
```
DIA-NN output, ~2.4M rows, 71 columns. Key columns: `Stripped.Sequence`, `Protein.Group`, `Run`, `Ms1.Area`, `Ms2.Area`.

### Mouse brain
```
My Drive/MS/Users/aleduc/AAS_rev/psm.tsv
```
FragPipe PSM output, ~242k rows. Key columns: `Peptide`, `Protein`, `Intensity` (precursor MS1), `Intensity mouse_brain_1_{channel}` (TMT reporters for channels 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N).

> Note: `PSM_all.tsv` also exists in the same folder but has zeroed-out TMT intensities due to a quantification issue — always use `psm.tsv`.

---

## iNeuron analysis (`Analysis/iNeuron.ipynb`)

### Experiment design
iPSC-derived neurons treated with degradation pathway inhibitors:

| Condition | Treatment |
|-----------|-----------|
| DMSO | Vehicle control |
| Proteasome | Proteasome inhibitor |
| Lysosome | Lysosomal/autophagy inhibitor |
| Both | Proteasome + Lysosomal inhibitor |

Two differentiation timepoints (7-day, 14-day), 4 replicates each = 32 raw files.

### RAAS computation (label-free)
```
RAAS = MS1_area(SAAP) / MS1_area(BP)    per run
```
Averaged within condition after log10 transform.

### Filters applied
- SAAP must be in `high_quality_SAAPs.xlsx`  ← currently active
- SAAP degradation rate ≤ BP degradation rate / 2 (from `SILAC.xlsx`) — SAAPs that degrade at least 2× slower than their BP, meaning they are stabilised relative to the canonical peptide

### Key outputs
1. Unique SAAPs per raw file / condition (dot plot)
2. Unique total peptides per condition
3. log10(RAAS) beeswarm per condition with median lines
4. Normalised RAAS boxplot (per-SAAP mean subtracted)
5. Scatter: RAAS fold change (Both vs DMSO) vs degradation fold change, coloured by log2(BP deg rate)

---

## Mouse brain analysis (`Analysis/MouseBrain.ipynb`)

### Experiment design
Mouse brain tissue, TMT10 multiplexed, young vs old comparison.

| Age group | TMT channels |
|-----------|-------------|
| Young | 126, 127N, 127C, 128N, 128C |
| Old | 129N, 129C, 130N, 130C, 131N |

> **TODO**: confirm channel→age mapping. The `channel_age` dict in the paths cell is a placeholder — update once the experiment annotation is verified.

### SAAP identification
SAAPs are identified by their `Protein` column in the PSM file:
- `MTP_` prefix entries (e.g. `sp|EK008270|MTP_EK008270-ek`)
- `new_peptide` entries (e.g. `sp|ST000001|new_peptide_1-st`)

Both types come from the custom SAAP fasta `output_MTP_cl_fp.fasta` used in the FragPipe search. There are **82,954 unique peptides** total (82,368 canonical + 586 SAAPs).

### SAAP → BP mapping
BPs are not encoded in the fasta headers. Instead, for each SAAP we find its BP by **1-AA Hamming matching**: search canonical PSM peptides of the same length for one that differs by exactly one amino acid. This recovers the large majority of SAAP→BP pairs directly from the experiment.

### RAAS computation (TMT)
TMT reporter intensities are relative within a peptide and cannot be directly compared between SAAP and BP. The correct formula is:

```
RAAS_channel = (TMT_SAAP_channel × Precursor_SAAP) / (TMT_BP_channel × Precursor_BP)
```

- `TMT_channel` = reporter ion intensity for that sample channel
- `Precursor` = MS1 precursor intensity (column `Intensity`) — scales between different peptides

This gives an absolute per-channel RAAS that is comparable across age groups.

### Key outputs
1. log2(RAAS) swarmplot: Young vs Old with median lines
2. Normalised log2(RAAS) boxplot (per-SAAP mean subtracted across channels)
3. Volcano plot: Welch t-test Old vs Young per SAAP, x = log2 FC, y = −log10(p-value)

---

## Important implementation notes for Claude

- **Do not use `fillna(0)` on TMT columns** — zeros are treated as missing and excluded from RAAS calculation. Keep them as NaN or filter with `> 0`.
- **Do not reintroduce `PSM_all.tsv`** — use `psm.tsv` for mouse brain data.
- **plt.savefig() calls are commented out** throughout both notebooks — do not uncomment unless asked.
- **Paths use `pathlib.Path`** — do not hardcode string paths.
- The iNeuron `ratios` dataframe is computed per-run (not per-condition); averaging across conditions happens downstream.
- The mouse brain `raas_df` has one row per SAAP × channel observation. Pivot with `pivot_table` before running statistics.
- When adding new plot cells, follow the existing style: `sns.despine()`, `plt.tight_layout()`, commented-out `savefig`, log2 scale for RAAS.
