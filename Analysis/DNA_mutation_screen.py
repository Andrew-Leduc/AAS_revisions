"""
DNA_mutation_screen.py

For each SAAP in Supplemental_Data_7, check whether the same amino acid change
at the same protein position in the same gene is also observed as a somatic/germline
missense mutation in any sample in all_missense_mutations.tsv.

A match means the SAAP *could* be explained by a DNA mutation rather than
genuine mistranslation — useful for flagging or excluding such cases.

Match criteria (all three must agree):
  - Ensembl gene ID
  - Protein position
  - Amino acid change (ref → alt)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
MISSENSE_PATH = Path('/Users/andrewleduc/Desktop/all_missense_mutations.tsv')
SAAP_PATH     = Path('/Users/andrewleduc/Downloads/Supplemental_Data_7.SAAP_coordinates.xlsx')
META_DIR      = Path('/Users/andrewleduc/Desktop/Github/AAS_revisions/meta_files')
SD2_PATH      = META_DIR / 'Supplemental_Data_2.SAAP_proteins.xlsx'
HQ_PATH       = META_DIR / 'high_quality_SAAPs.xlsx'
OUT_PATH      = Path('/Users/andrewleduc/Desktop/Github/AAS_revisions/Analysis/SAAP_mutation_matches.tsv')

# ── 1. Load annotation files ─────────────────────────────────────────────────
print('Loading annotation files...')
sd2 = pd.read_excel(SD2_PATH, usecols=['SAAP', 'Mean precursor RAAS'])
mean_raas = sd2.groupby('SAAP')['Mean precursor RAAS'].mean().rename('mean_RAAS')

hq_saaps = set(pd.read_excel(HQ_PATH, usecols=['SAAP'])['SAAP'].dropna())
print(f'  SD2: {len(mean_raas)} unique SAAPs with RAAS')
print(f'  High-quality SAAP list: {len(hq_saaps)} entries')

# ── 2. Load SAAP coordinates ──────────────────────────────────────────────────
print('Loading SAAP coordinates...')
saap = pd.read_excel(SAAP_PATH, usecols=['SAAP', 'BP', 'fromto', 'protein.position', 'gene', 'chr', 'coor'])
saap = saap.dropna(subset=['gene', 'protein.position', 'fromto'])

# Parse "P:N" → ref_aa='P', alt_aa='N'
saap[['ref_aa', 'alt_aa']] = saap['fromto'].str.split(':', expand=True)
saap['protein.position'] = pd.to_numeric(saap['protein.position'], errors='coerce').dropna().astype(int)
saap = saap.dropna(subset=['protein.position'])
saap['protein.position'] = saap['protein.position'].astype(int)

print(f'  {len(saap)} SAAP entries, {saap["gene"].nunique()} unique genes')

# Build lookup set: (gene, position, ref, alt)  — O(1) per missense row
saap_keys = set(
    zip(saap['gene'], saap['protein.position'], saap['ref_aa'], saap['alt_aa'])
)
saap_genes = set(saap['gene'])
print(f'  {len(saap_keys)} unique (gene, position, change) keys to screen')

# ── 2. Stream missense file in chunks ─────────────────────────────────────────
# Only load columns needed for matching to keep memory low
USECOLS   = ['sample_id', 'Gene', 'Protein_position', 'Amino_acids', 'HGVSp', 'VAF', 'gnomADe_AF']
CHUNKSIZE = 500_000

print(f'\nStreaming {MISSENSE_PATH.name} in chunks of {CHUNKSIZE:,}...')

matches = []
total_rows = 0

reader = pd.read_csv(
    MISSENSE_PATH,
    sep='\t',
    usecols=USECOLS,
    dtype={'Gene': str, 'Protein_position': str, 'Amino_acids': str},
    chunksize=CHUNKSIZE,
    low_memory=False,
)

for i, chunk in enumerate(reader):
    total_rows += len(chunk)

    # Fast pre-filter: keep only rows whose gene is in our SAAP gene set
    chunk = chunk[chunk['Gene'].isin(saap_genes)].copy()
    if chunk.empty:
        if (i + 1) % 5 == 0:
            print(f'  chunk {i+1}: {total_rows:,} rows processed...')
        continue

    # Parse protein position (may be "162" or "162-163"; take first number)
    chunk['pos_int'] = pd.to_numeric(
        chunk['Protein_position'].str.extract(r'(\d+)', expand=False),
        errors='coerce'
    )
    chunk = chunk.dropna(subset=['pos_int'])
    chunk['pos_int'] = chunk['pos_int'].astype(int)

    # Parse amino acids: "T/A" → ref='T', alt='A'
    aa_split = chunk['Amino_acids'].str.split('/', expand=True)
    chunk['mut_ref'] = aa_split[0]
    chunk['mut_alt'] = aa_split[1] if 1 in aa_split.columns else np.nan
    chunk = chunk.dropna(subset=['mut_ref', 'mut_alt'])

    # Check each row against SAAP lookup set
    hit_mask = chunk.apply(
        lambda r: (r['Gene'], r['pos_int'], r['mut_ref'], r['mut_alt']) in saap_keys,
        axis=1
    )
    hits = chunk[hit_mask]
    if not hits.empty:
        matches.append(hits)

    if (i + 1) % 5 == 0:
        print(f'  chunk {i+1}: {total_rows:,} rows processed, {sum(len(m) for m in matches):,} matches so far...')

print(f'\nFinished. {total_rows:,} missense rows screened.')

# ── 3. Compile results ────────────────────────────────────────────────────────
if matches:
    mut_hits = pd.concat(matches, ignore_index=True)
    # Join back to SAAP metadata
    saap_meta = saap[['SAAP', 'BP', 'fromto', 'protein.position', 'gene', 'chr', 'coor']].copy()
    saap_meta = saap_meta.rename(columns={'protein.position': 'pos_int', 'gene': 'Gene'})
    saap_meta['pos_int'] = saap_meta['pos_int'].astype(int)

    # Parse fromto to get ref/alt for join
    saap_meta[['mut_ref', 'mut_alt']] = saap_meta['fromto'].str.split(':', expand=True)

    result = mut_hits.merge(
        saap_meta[['SAAP', 'BP', 'fromto', 'Gene', 'pos_int', 'mut_ref', 'mut_alt', 'chr', 'coor']],
        on=['Gene', 'pos_int', 'mut_ref', 'mut_alt'],
        how='left'
    )

    # Collapse to one row per SAAP: count patients, average VAF
    summary = (result.groupby(['SAAP', 'BP', 'fromto', 'Gene', 'chr', 'coor'])
                     .agg(n_patients=('sample_id', 'nunique'),
                          mean_VAF=('VAF', 'mean'),
                          gnomADe_AF=('gnomADe_AF', 'mean'))
                     .reset_index()
                     .sort_values('n_patients', ascending=False))

    # Annotate with mean RAAS and high-quality flag
    summary = summary.join(mean_raas, on='SAAP')
    summary['high_quality'] = summary['SAAP'].isin(hq_saaps)

    print(f'\n── Results ──────────────────────────────────────────────────')
    print(f'  SAAPs with at least one matching DNA mutation : {len(summary)}')
    print(f'  Unique patient samples (across all SAAPs)     : {result["sample_id"].nunique()}')
    print(f'  Of matched SAAPs, high-quality                : {summary["high_quality"].sum()}')
    print(f'  Of matched SAAPs, with RAAS data              : {summary["mean_RAAS"].notna().sum()}')

    summary.to_csv(OUT_PATH, sep='\t', index=False)
    print(f'\nResults written to: {OUT_PATH}')

    print('\nTop mutation-explained SAAPs (by number of patients):')
    print(summary.head(20).to_string(index=False))

    # ── Correlations ──────────────────────────────────────────────────────────
    from scipy.stats import pearsonr

    plot_df = summary.dropna(subset=['mean_RAAS']).copy()

    r_pat, p_pat = pearsonr(plot_df['n_patients'], plot_df['mean_RAAS'])
    gnomad_df = plot_df.dropna(subset=['gnomADe_AF'])
    r_gnomad, p_gnomad = pearsonr(gnomad_df['gnomADe_AF'], gnomad_df['mean_RAAS'])

    print(f'\nPearson r (RAAS vs n_patients) : r={r_pat:.3f}, p={p_pat:.3e}')
    print(f'Pearson r (RAAS vs gnomADe_AF) : r={r_gnomad:.3f}, p={p_gnomad:.3e}')

    # ── Plot: mean RAAS vs n_patients, coloured by high_quality ──────────────
    palette = {True: '#E64B35', False: '#8C8C8C'}
    labels  = {True: 'High quality', False: 'Other'}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, x_col, xlabel, r, p in [
        (axes[0], 'n_patients', 'Number of patients with matching DNA mutation', r_pat,    p_pat),
        (axes[1], 'gnomADe_AF', 'gnomAD allele frequency',                       r_gnomad, p_gnomad),
    ]:
        df_ax = plot_df.dropna(subset=[x_col])
        for hq, grp in df_ax.groupby('high_quality'):
            ax.scatter(grp[x_col], grp['mean_RAAS'],
                       c=palette[hq], label=labels[hq],
                       s=40, alpha=0.8, edgecolors='none')
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel('Mean precursor RAAS (log10)', fontsize=11)
        ax.set_title(f'r = {r:.3f}, p = {p:.2e}', fontsize=10)
        ax.legend(fontsize=9)
        sns.despine(ax=ax)

    fig.suptitle('Mutation-explainable SAAPs: RAAS correlations', fontsize=13)
    plt.tight_layout()
    plot_path = OUT_PATH.parent / 'SAAP_mutation_RAAS_vs_patients.jpeg'
    plt.savefig(plot_path, bbox_inches='tight', dpi=150)
    print(f'Plot saved to: {plot_path}')
    plt.show()

else:
    print('\nNo SAAPs matched any missense mutation.')
