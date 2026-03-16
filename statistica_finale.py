# -*- coding: utf-8 -*-
# =============================================================================
# PIPELINE UNIFICATA — Cardiomiopatie Varianti Missense
# Versione: 5.0  (ClinVar fix + Gradiente v2 + Robustezza + Domini + Figure
#                  + Predittività gene + IEV + Identikit + VUS candidate)
#
# CHANGELOG v5.0:
#   - Aggiunta: Sezione 12 — Predittività per gene (AUC-ROC P+LP vs B+LB)
#   - Aggiunta: Sezione 13 — Banda IEV (AF 1e-4 — 1e-2), composizione ClinVar
#   - Aggiunta: Sezione 14 — Identikit patogenico gene-specifico + VUS candidate
#   - Master metrics aggiornato con nuove variabili
#   - Figure rinumerate a Sezione 16
#
# CHANGELOG v4.1 (post-audit):
#   - Fix: warnings('once') invece di ('ignore')
#   - Fix: conf_by_cl .agg() sintassi corretta (stringhe, non tuple)
#   - Fix: set_xticks() prima di set_xticklabels() in fig1
#   - Fix: df_full salvato come .csv.gz (non .xlsx lento con 42k righe)
#   - Fix: bbox sulle annotazioni volcano fig5
#   - Miglioramento: bootstrap n=100 per stabilità ARI
#   - Aggiunta: test Shapiro-Wilk (giustifica Kruskal-Wallis)
#   - Aggiunta: test monotonia Spearman cluster ~ score_composite/AM/REVEL
#   - Aggiunta: analisi VUS nello spazio score (evidenza zona grigia)
#
# INCLUDE:
#   01_statistica  — ClinVar fix, BIC, GMM, propagazione, metriche qualità
#   02_validazione — gradiente ClinVar, OR Fisher, VUS risk, AF validation
#   03_sensitivity — bootstrap (n=100), seed (n=10), CV (n=20)
#   04_domain      — enrichment Fisher + FDR, conflicting x dominio
#   05_gene        — AUC-ROC per gene, IEV band, identikit + VUS candidate
#   gradiente v2   — Spearman ordinale, Kruskal + post-hoc Bonferroni,
#                    AF per cluster con IEV, heatmap P(path) AF×REVEL
#   Figure 1–6
#
# OTTIMIZZAZIONI MEMORIA:
#   - Caricamento colonne minime + dtype float32 per score
#   - ClinVar fix vettorizzato (nessun .apply(axis=1))
#   - del + gc.collect() dopo ogni sezione pesante
#   - Bootstrap/CV senza conservare i modelli
#   - Figure su campione (≤6000 punti), plt.close() immediato
#   - FDR Benjamini-Hochberg inline (no statsmodels)
# =============================================================================

import gc
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from scipy.stats import (fisher_exact, spearmanr, kruskal,
                          mannwhitneyu)
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (silhouette_score, davies_bouldin_score,
                              adjusted_rand_score)
from sklearn.metrics import pairwise_distances_argmin
from sklearn.utils import resample
from sklearn.model_selection import train_test_split

warnings.filterwarnings('once')   # 'once' per non mascherare errori reali
sys.stdout.reconfigure(encoding='utf-8')
sns.set_theme(style='whitegrid', font_scale=1.0)

# =============================================================================
# CONFIGURAZIONE
# =============================================================================

INPUT_FILE  = 'cardiomyopathy_variants_with_interpro_local.xlsx'  # fallback .csv
OUTDIR      = 'PIPELINE_RESULTS'
os.makedirs(OUTDIR, exist_ok=True)

# Colonne obbligatorie — tutto il resto viene caricato ma non duplicato
AM_COL    = 'am_pathogenicity'
REVEL_COL = 'REVEL_max'
AF_COL    = 'AF_grpmax'
LOG_AF    = 'log10_AF_grpmax'
POS_COL   = 'aa_pos'
LEN_COL   = 'Protein_Length'
GENE_COL  = 'GeneSymbol'

ORDINE_CLASSI = [
    'Pathogenic', 'Pathogenic/Likely_path', 'Likely_pathogenic',
    'VUS', 'Conflicting', 'Likely_benign', 'Benign/Likely_benign', 'Benign'
]
PALETTE_CLASS = {
    'Pathogenic':             '#d62728',
    'Pathogenic/Likely_path': '#ff7f0e',
    'Likely_pathogenic':      '#e377c2',
    'VUS':                    '#7f7f7f',
    'Conflicting':            '#bcbd22',
    'Likely_benign':          '#17becf',
    'Benign/Likely_benign':   '#aec7e8',
    'Benign':                 '#1f77b4',
}
PALETTE_CLUSTER = {0: '#1f77b4', 1: '#ff7f0e', 2: '#2ca02c', 3: '#d62728'}

# =============================================================================
# UTIL: Benjamini-Hochberg FDR (evita statsmodels)
# =============================================================================

def bh_fdr(pvalues):
    """Corregge p-values con metodo Benjamini-Hochberg."""
    pv  = np.asarray(pvalues, dtype=float)
    n   = len(pv)
    if n == 0:
        return pv
    order  = np.argsort(pv)
    adj    = np.empty(n, dtype=float)
    # Calcola adj = p * n / rank, poi enforcea monotonia dal fondo
    sorted_pv = pv[order]
    adj_sorted = np.minimum(1.0, sorted_pv * n / np.arange(1, n + 1))
    # Monotonia: dal penultimo al primo, adj[i] = min(adj[i], adj[i+1])
    for i in range(n - 2, -1, -1):
        adj_sorted[i] = min(adj_sorted[i], adj_sorted[i + 1])
    # Rimetti nell'ordine originale
    adj[order] = adj_sorted
    return adj


def _bh_safe(df, col):
    """Aggiunge colonna padj a DataFrame con colonna 'col' di p-values."""
    valid = df[col].notna()
    if valid.sum() == 0:
        df['padj'] = np.nan
        return df
    padj = np.ones(len(df))
    padj[valid.values] = bh_fdr(df.loc[valid, col].values)
    df['padj'] = padj
    return df

# =============================================================================
# SEZIONE 1 — CARICAMENTO + RINOMINA + OTTIMIZZAZIONE DTYPE
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 1 — CARICAMENTO DATASET')
print('='*70)

try:
    df = pd.read_excel(INPUT_FILE)
    print(f'  Caricato Excel: {INPUT_FILE}')
except FileNotFoundError:
    csv_path = INPUT_FILE.replace('.xlsx', '.csv')
    df = pd.read_csv(csv_path, sep=None, engine='python', encoding='latin1')
    df.columns = df.columns.str.strip()
    print(f'  Caricato CSV fallback: {csv_path}')

if '#AlleleID' in df.columns:
    df = df.rename(columns={'#AlleleID': 'AlleleID'})

# Rinomina colonne AF/pos se necessario
rename_map = {
    'AF_grpmax':           AF_COL,
    'log10_AF_grpmax':     LOG_AF,
    'aa_pos_norm_protein': 'aa_pos_norm',
}
for old, new in rename_map.items():
    if old in df.columns and new not in df.columns:
        df = df.rename(columns={old: new})

# Conversione numerica + float32 per score (risparmio ~50% RAM)
for col in [AM_COL, REVEL_COL, AF_COL, LOG_AF]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').astype('float32')

for col in [POS_COL, LEN_COL]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

# Colonna log10AF (usa preesistente o calcola)
if LOG_AF not in df.columns and AF_COL in df.columns:
    df[LOG_AF] = np.log10(df[AF_COL].replace(0, np.nan)).astype('float32')

# Score composito
df['score_composite'] = np.sqrt(df[AM_COL] * df[REVEL_COL]).astype('float32')

# Posizione normalizzata
if 'aa_pos_norm' not in df.columns and POS_COL in df.columns and LEN_COL in df.columns:
    df['aa_pos_norm'] = (df[POS_COL] / df[LEN_COL]).astype('float32')

print(f'  Varianti totali: {len(df):,}')
print(f'  Colonne: {len(df.columns)}')
print(f'  RAM stimata: ~{df.memory_usage(deep=True).sum() / 1e6:.0f} MB')

# =============================================================================
# SEZIONE 2 — CLINVAR FIX (vettorizzato, nessun .apply)
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 2 — CLINVAR FIX (categorie mutually exclusive)')
print('='*70)

# Recupera flag booleani esistenti
def _flag(col):
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].fillna(0).astype(bool)

is_conf  = _flag('is_conflicting')
is_path  = _flag('is_pathogenic')
is_lp    = _flag('is_likely_pathogenic')
is_vus   = _flag('is_vus')
is_lb    = _flag('is_likely_benign')
is_b     = _flag('is_benign')

# Diagnostico doppio conteggio
n_both = int((is_path & is_conf).sum())
print(f'  is_pathogenic AND is_conflicting (doppio conteggio): {n_both:,}')

# Assegna categoria esclusiva — is_conflicting ha PRIORITÀ MASSIMA
cat = pd.Series('unclassified', index=df.index)
cat[is_b   & ~is_conf] = 'benign'
cat[is_lb  & ~is_conf] = 'likely_benign'
cat[is_vus & ~is_conf] = 'vus'
cat[is_lp  & ~is_conf] = 'likely_pathogenic'
cat[is_path & ~is_conf] = 'pathogenic'
cat[is_conf]            = 'conflicting'        # sovrascrive tutto
df['clinvar_category']  = cat.astype('category')

# Flag derivati utili
df['is_pathogenic_clean'] = (is_path & ~is_conf).astype('int8')
df['is_path_group']       = ((is_path | is_lp) & ~is_conf).astype('int8')
df['is_benign_group']     = (is_b | is_lb).astype('int8')
df['is_conflicting_any']  = is_conf.astype('int8')

# Etichetta leggibile per le figure (Capitalized, con '/')
_lbl_map = {
    'pathogenic':        'Pathogenic',
    'likely_pathogenic': 'Likely_pathogenic',
    'vus':               'VUS',
    'conflicting':       'Conflicting',
    'likely_benign':     'Likely_benign',
    'benign':            'Benign',
    'unclassified':      'Other',
}
df['class_clean'] = df['clinvar_category'].map(_lbl_map).fillna('Other')

print('\n  Distribuzione corretta (mutually exclusive):')
for cat_name, n in df['clinvar_category'].value_counts().items():
    print(f'    {cat_name:<30} {n:>6,}  ({100*n/len(df):.1f}%)')

del is_conf, is_path, is_lp, is_vus, is_lb, is_b, cat
gc.collect()

# Salva distribuzione corretta
_dist = df['clinvar_category'].value_counts().reset_index()
_dist.columns = ['category', 'n']
_dist['pct'] = (_dist['n'] / len(df) * 100).round(1)
_dist.to_csv(f'{OUTDIR}/clinvar_distribution_corrected.csv', index=False)

# =============================================================================
# SEZIONE 3 — GMM: BIC, CLUSTERING, LABEL PROPAGATION
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 3 — GMM CLUSTERING (BIC + StandardScaler)')
print('='*70)

# Subset con AM + REVEL disponibili
mask_gmm     = df[AM_COL].notna() & df[REVEL_COL].notna()
df_gmm       = df[mask_gmm].copy()
df_nocluster = df[~mask_gmm].copy()
print(f'  Subset GMM (AM+REVEL): {len(df_gmm):,}')
print(f'  Da propagare:          {len(df_nocluster):,}')

X_raw    = df_gmm[[AM_COL, REVEL_COL]].values.astype('float32')
scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X_raw)

# Selezione k via BIC con elbow 10%
print('\n  BIC model selection k=1–8:')
bic_scores, aic_scores = [], []
for k in range(1, 9):
    _gmm = GaussianMixture(n_components=k, random_state=42, n_init=10)
    _gmm.fit(X_scaled)
    bic_scores.append(_gmm.bic(X_scaled))
    aic_scores.append(_gmm.aic(X_scaled))
    print(f'    k={k}  BIC={bic_scores[-1]:.1f}  AIC={aic_scores[-1]:.1f}')
    del _gmm

pd.DataFrame({'k': range(1, 9), 'BIC': bic_scores, 'AIC': aic_scores}
             ).to_csv(f'{OUTDIR}/bic_model_selection.csv', index=False)

bic_deltas = [bic_scores[i-1] - bic_scores[i] for i in range(1, len(bic_scores))]
max_delta  = max(bic_deltas)
elbow_k    = next((i+2 for i, d in enumerate(bic_deltas) if d < 0.10*max_delta),
                   int(np.argmin(bic_scores)) + 1)
best_k_abs = int(np.argmin(bic_scores)) + 1
best_k     = elbow_k if best_k_abs > 6 else best_k_abs
print(f'\n  BIC min assoluto: k={best_k_abs}  |  Elbow (10%): k={elbow_k}  |  USATO: k={best_k}')

# Fit finale
gmm_final               = GaussianMixture(n_components=best_k, random_state=42, n_init=10)
df_gmm['cluster']       = gmm_final.fit_predict(X_scaled)

# Riordina cluster per AM crescente (0=benigno, max=HR)
c_means   = df_gmm.groupby('cluster')[AM_COL].mean().sort_values()
c_map     = {old: new for new, old in enumerate(c_means.index)}
df_gmm['cluster']    = df_gmm['cluster'].map(c_map).astype('int8')
HR_CLUSTER           = int(df_gmm['cluster'].max())
df_gmm['is_hr']      = (df_gmm['cluster'] == HR_CLUSTER).astype('int8')
df_gmm['cluster_propagated'] = np.int8(0)

print('\n  Cluster (AM crescente):')
for c in range(best_k):
    sub = df_gmm[df_gmm['cluster'] == c]
    tag = '  ← HIGH RISK' if c == HR_CLUSTER else ''
    print(f'    C{c}: n={len(sub):>6,}  AM={sub[AM_COL].mean():.3f}'
          f'  REVEL={sub[REVEL_COL].mean():.3f}{tag}')

# Metriche qualità
sil = silhouette_score(X_scaled, df_gmm['cluster'])
db  = davies_bouldin_score(X_scaled, df_gmm['cluster'])
rho_feat, p_feat = spearmanr(df_gmm[AM_COL], df_gmm[REVEL_COL])
print(f'\n  Silhouette={sil:.4f}  Davies-Bouldin={db:.4f}  Spearman(AM,REVEL)={rho_feat:.4f} p={p_feat:.2e}')

clusters_original = df_gmm['cluster'].values.copy()   # conserva per robustezza

# Label propagation su nearest centroid (aa_pos_norm)
if len(df_nocluster) > 0:
    centroids  = np.array([df_gmm[df_gmm['cluster'] == c]['aa_pos_norm'].mean()
                            for c in range(best_k)]).reshape(-1, 1)
    X_prop     = df_nocluster['aa_pos_norm'].fillna(0.5).values.reshape(-1, 1)
    prop_labels = pairwise_distances_argmin(X_prop, centroids)
    df_nocluster = df_nocluster.copy()
    df_nocluster['cluster']            = prop_labels.astype('int8')
    df_nocluster['is_hr']              = (df_nocluster['cluster'] == HR_CLUSTER).astype('int8')
    df_nocluster['cluster_propagated'] = np.int8(1)
    print(f'\n  Propagazione completata ({len(df_nocluster):,} varianti)')

# Merge dataset completo
common_cols = sorted(set(df_gmm.columns) & set(df_nocluster.columns))
df_full     = pd.concat([df_gmm[common_cols], df_nocluster[common_cols]], ignore_index=True)
print(f'\n  Dataset completo: {len(df_full):,}  (GMM={len(df_gmm):,}  prop={len(df_nocluster):,})')

del df_nocluster, X_raw
gc.collect()

# Salva dataset completo (CSV gzip — più veloce di xlsx con 42k righe)
df_full.to_csv(f'{OUTDIR}/clustered_dataset.csv.gz', index=False, compression='gzip')

# Distribuzione ClinVar per cluster
cc_rows = []
for c in range(best_k):
    sub  = df_gmm[df_gmm['cluster'] == c]
    cats = sub['clinvar_category'].value_counts().to_dict()
    row  = {'cluster': c, 'n_total': len(sub)}
    for cat_name in ['pathogenic','likely_pathogenic','vus',
                     'conflicting','likely_benign','benign','unclassified']:
        row[f'n_{cat_name}'] = cats.get(cat_name, 0)
    row['n_path_group']  = row['n_pathogenic'] + row['n_likely_pathogenic']
    row['n_benign_group'] = row['n_likely_benign'] + row['n_benign']
    cc_rows.append(row)
cc_df = pd.DataFrame(cc_rows)
cc_df.to_csv(f'{OUTDIR}/clinvar_by_cluster.csv', index=False)

# =============================================================================
# SEZIONE 4 — ANALISI 1: GRADIENTE ORDINALE (Spearman su rank ClinVar)
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 4 — GRADIENTE ORDINALE (Spearman rank ClinVar)')
print('='*70)

rank_map = {
    'Benign': 0, 'Benign/Likely_benign': 1, 'Likely_benign': 2,
    'VUS': 3, 'Conflicting': 4,
    'Likely_pathogenic': 5, 'Pathogenic/Likely_path': 6, 'Pathogenic': 7
}
df_gmm['rank'] = df_gmm['class_clean'].map(rank_map)
valid_rank = df_gmm.dropna(subset=['rank', AM_COL, REVEL_COL, 'score_composite'])

print(f'\n  {"Classe":<30} {"n":>6}  {"AM med":>8}  {"REVEL med":>10}')
print(f'  {"-"*58}')
for cls in ORDINE_CLASSI:
    sub = valid_rank[valid_rank['class_clean'] == cls]
    if len(sub) < 10:
        continue
    print(f'  {cls:<30} {len(sub):>6}  {sub[AM_COL].median():>8.3f}  {sub[REVEL_COL].median():>10.3f}')

rho_am,   p_am   = spearmanr(valid_rank['rank'], valid_rank[AM_COL])
rho_rv,   p_rv   = spearmanr(valid_rank['rank'], valid_rank[REVEL_COL])
rho_comp, p_comp = spearmanr(valid_rank['rank'], valid_rank['score_composite'])
print(f'\n  Spearman (n={len(valid_rank):,}):')
print(f'    AM              rho={rho_am:.4f}  p={p_am:.2e}')
print(f'    REVEL           rho={rho_rv:.4f}  p={p_rv:.2e}')
print(f'    Score composito rho={rho_comp:.4f}  p={p_comp:.2e}')

del valid_rank
gc.collect()

# =============================================================================
# SEZIONE 5 — ANALISI 2: PROFILI CLUSTER + GRADIENTE CLINVAR
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 5 — PROFILI CLUSTER + GRADIENTE CLINVAR')
print('='*70)

cluster_ids = list(range(best_k))

# Profili
prof_rows = []
for c in cluster_ids:
    sub = df_gmm[df_gmm['cluster'] == c]
    row = {
        'cluster':    c,
        'n':          len(sub),
        'pct_tot':    100 * len(sub) / len(df_gmm),
        'am_med':     float(sub[AM_COL].median()),
        'revel_med':  float(sub[REVEL_COL].median()),
        'pct_path':   100 * sub['is_pathogenic_clean'].mean(),
        'pct_path_group': 100 * sub['is_path_group'].mean(),
        'pct_benign': 100 * sub['is_benign_group'].mean(),
        'pct_vus':    100 * (sub['clinvar_category'] == 'vus').mean(),
        'pct_conf':   100 * sub['is_conflicting_any'].mean(),
    }
    prof_rows.append(row)
    tag = '  ← HR' if c == HR_CLUSTER else ''
    print(f'  C{c}: n={len(sub):>6,}  AM={row["am_med"]:.3f}  REVEL={row["revel_med"]:.3f}'
          f'  %P={row["pct_path"]:.1f}  %VUS={row["pct_vus"]:.1f}'
          f'  %Conf={row["pct_conf"]:.1f}{tag}')

prof_df = pd.DataFrame(prof_rows)
prof_df.to_csv(f'{OUTDIR}/cluster_profiles.csv', index=False)

# Gradiente Spearman cluster ~ proporzione patogenici
grad_rows = []
for c in cluster_ids:
    sub  = df_gmm[df_gmm['cluster'] == c]
    cats = sub['clinvar_category'].value_counts().to_dict()
    n    = len(sub)
    grad_rows.append({
        'cluster':         c,
        'n':               n,
        'pathogenic_prop': cats.get('pathogenic', 0) / n,
        'path_group_prop': (cats.get('pathogenic', 0) + cats.get('likely_pathogenic', 0)) / n,
        'conflicting_prop': cats.get('conflicting', 0) / n,
    })
grad_df = pd.DataFrame(grad_rows)
rho_g,  p_g  = spearmanr(grad_df['cluster'], grad_df['pathogenic_prop'])
rho_gg, p_gg = spearmanr(grad_df['cluster'], grad_df['path_group_prop'])
print(f'\n  Gradiente Spearman  (P puro):   rho={rho_g:.4f}  p={p_g:.2e}')
print(f'  Gradiente Spearman  (P+LP):     rho={rho_gg:.4f}  p={p_gg:.2e}')
grad_df.to_csv(f'{OUTDIR}/clinvar_gradient.csv', index=False)

# Test monotonia: Spearman cluster ~ score (su varianti individuali)
rho_mono_comp, p_mono_comp = spearmanr(df_gmm['cluster'], df_gmm['score_composite'])
rho_mono_am,   p_mono_am   = spearmanr(df_gmm['cluster'], df_gmm[AM_COL])
rho_mono_rv,   p_mono_rv   = spearmanr(df_gmm['cluster'], df_gmm[REVEL_COL])
print(f'\n  Monotonia cluster ~ score (varianti individuali, n={len(df_gmm):,}):')
print(f'    score_composite  rho={rho_mono_comp:.4f}  p={p_mono_comp:.2e}')
print(f'    AlphaMissense    rho={rho_mono_am:.4f}    p={p_mono_am:.2e}')
print(f'    REVEL            rho={rho_mono_rv:.4f}    p={p_mono_rv:.2e}')

# Verifica monotonia stretta sulle mediane per cluster
med_by_cl = df_gmm.groupby('cluster')[['score_composite', AM_COL, REVEL_COL]].median()
is_strict_mono = all(med_by_cl['score_composite'].diff().dropna() > 0)
print(f'    Mediane score_composite crescenti: {med_by_cl["score_composite"].values.round(3).tolist()}')
print(f'    Monotonia stretta: {"SI" if is_strict_mono else "NO"}')

mono_df = pd.DataFrame({
    'score':       ['score_composite', AM_COL, REVEL_COL],
    'spearman_rho': [rho_mono_comp, rho_mono_am, rho_mono_rv],
    'p_value':      [p_mono_comp, p_mono_am, p_mono_rv],
    'strict_monotone_medians': [is_strict_mono,
        all(med_by_cl[AM_COL].diff().dropna() > 0),
        all(med_by_cl[REVEL_COL].diff().dropna() > 0)],
})
mono_df.to_csv(f'{OUTDIR}/cluster_monotonicity.csv', index=False)

# Enrichment OR Fisher HR vs benigno
hr  = df_gmm[df_gmm['cluster'] == HR_CLUSTER]
ben = df_gmm[df_gmm['cluster'] == 0]

def _fisher_2x2(col, g1, g2):
    t = np.array([[g1[col].sum(), len(g1) - g1[col].sum()],
                  [g2[col].sum(), len(g2) - g2[col].sum()]])
    return fisher_exact(t)

OR_p,   p_p   = _fisher_2x2('is_pathogenic_clean', hr, ben)
OR_pg,  p_pg  = _fisher_2x2('is_path_group',       hr, ben)
print(f'\n  OR Fisher P puro  (HR vs C0): {OR_p:.2f}  p={p_p:.2e}')
print(f'  OR Fisher P+LP    (HR vs C0): {OR_pg:.2f}  p={p_pg:.2e}')
pd.DataFrame({'test':       ['pathogenic_puro', 'path_group'],
              'odds_ratio': [OR_p,  OR_pg],
              'p_value':    [p_p,   p_pg]
              }).to_csv(f'{OUTDIR}/enrichment_HR.csv', index=False)

# =============================================================================
# SEZIONE 6 — ANALISI 3: KRUSKAL-WALLIS + POST-HOC BONFERRONI
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 6 — KRUSKAL-WALLIS + POST-HOC BONFERRONI')
print('='*70)

# Test di normalità Shapiro-Wilk (giustifica la scelta di Kruskal-Wallis)
from scipy.stats import shapiro

print('\n  Test normalità Shapiro-Wilk (campione ≤5000):')
for col, lbl in [(AM_COL, 'AlphaMissense'), (REVEL_COL, 'REVEL')]:
    samp = df_gmm[col].dropna()
    samp = samp.sample(min(5000, len(samp)), random_state=42)
    _, p_norm = shapiro(samp)
    verdict = 'NON normale' if p_norm < 0.05 else 'normale'
    print(f'    Shapiro {lbl}: p={p_norm:.2e} → {verdict}')
print('  → Kruskal-Wallis appropriato (distribuzioni non normali)\n')

groups_am = [df_gmm[df_gmm['cluster'] == c][AM_COL].dropna().values
             for c in cluster_ids]
groups_rv = [df_gmm[df_gmm['cluster'] == c][REVEL_COL].dropna().values
             for c in cluster_ids]

H_am, p_kw_am = kruskal(*groups_am)
H_rv, p_kw_rv = kruskal(*groups_rv)
print(f'\n  Kruskal-Wallis AM    H={H_am:.1f}  p={p_kw_am:.2e}')
print(f'  Kruskal-Wallis REVEL H={H_rv:.1f}  p={p_kw_rv:.2e}')

coppie  = list(combinations(cluster_ids, 2))
n_comp  = len(coppie)
ph_rows = []
print(f'\n  Post-hoc Mann-Whitney + Bonferroni (n_comp={n_comp}):')
print(f'  {"Coppia":<14} {"AM p_bonf":>12}    {"REVEL p_bonf":>14}')
for c1, c2 in coppie:
    _, p_am = mannwhitneyu(groups_am[c1], groups_am[c2], alternative='two-sided')
    _, p_rv = mannwhitneyu(groups_rv[c1], groups_rv[c2], alternative='two-sided')
    p_am_b  = min(p_am * n_comp, 1.0)
    p_rv_b  = min(p_rv * n_comp, 1.0)
    sig_am  = '***' if p_am_b < 1e-3 else ('*' if p_am_b < 0.05 else 'ns')
    sig_rv  = '***' if p_rv_b < 1e-3 else ('*' if p_rv_b < 0.05 else 'ns')
    print(f'  C{c1} vs C{c2}          {p_am_b:>12.2e} {sig_am:>4}  {p_rv_b:>14.2e} {sig_rv}')
    ph_rows.append({'c1': c1, 'c2': c2, 'p_am_bonf': p_am_b, 'p_revel_bonf': p_rv_b})

pd.DataFrame(ph_rows).to_csv(f'{OUTDIR}/posthoc_bonferroni.csv', index=False)

# =============================================================================
# SEZIONE 7 — ANALISI 4: VUS RISK FORMULA CORRETTA
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 7 — VUS RISK PER CLUSTER (no conflicting al denominatore)')
print('='*70)

vus_rows = []
for c in cluster_ids:
    sub   = df_gmm[df_gmm['cluster'] == c]
    cats  = sub['clinvar_category'].value_counts().to_dict()
    n_p   = cats.get('pathogenic', 0)
    n_lp  = cats.get('likely_pathogenic', 0)
    n_b   = cats.get('benign', 0)
    n_lb  = cats.get('likely_benign', 0)
    n_vus = cats.get('vus', 0)
    n_cf  = cats.get('conflicting', 0)
    n_cls = n_p + n_lp + n_b + n_lb   # denominatore corretto (no conflicting, no vus)
    p_est = n_p / n_cls if n_cls > 0 else np.nan
    vus_rows.append({'cluster': c, 'n_total': len(sub),
                     'n_pathogenic': n_p, 'n_likely_path': n_lp,
                     'n_benign': n_b, 'n_likely_benign': n_lb,
                     'n_vus': n_vus, 'n_conflicting': n_cf,
                     'n_classified': n_cls, 'p_pathogenic_estimate': p_est})

vus_df = pd.DataFrame(vus_rows)
print(f'\n  {"Cl":>3} {"n":>7} {"n_P":>5} {"n_Conf":>7} {"n_VUS":>7}'
      f'  {"n_class":>8}  {"p_path_est":>12}')
for _, r in vus_df.iterrows():
    print(f'  {int(r.cluster):>3} {int(r.n_total):>7,} {int(r.n_pathogenic):>5}'
          f' {int(r.n_conflicting):>7} {int(r.n_vus):>7}'
          f'  {int(r.n_classified):>8}  {r.p_pathogenic_estimate:>12.3f}')
vus_df.to_csv(f'{OUTDIR}/vus_risk_by_cluster.csv', index=False)

# --- VUS nello spazio score: evidenza della zona grigia ---
print('\n  Analisi VUS nello spazio score:')
vus_mask = df_gmm['clinvar_category'] == 'vus'
non_vus_mask = df_gmm['clinvar_category'].isin(['pathogenic', 'likely_pathogenic',
                                                  'benign', 'likely_benign'])
vus_sub  = df_gmm[vus_mask]
nvus_sub = df_gmm[non_vus_mask]

# Distribuzione VUS per cluster
vus_per_cl = vus_sub['cluster'].value_counts().sort_index()
vus_pct_cl = (vus_sub.groupby('cluster').size()
              / df_gmm.groupby('cluster').size() * 100).round(1)

print(f'\n  {"Cl":>3}  {"n_VUS":>7}  {"% VUS nel cl":>13}  {"AM med VUS":>11}  {"REVEL med VUS":>14}')
vus_score_rows = []
for c in cluster_ids:
    v = vus_sub[vus_sub['cluster'] == c]
    n_cl = len(df_gmm[df_gmm['cluster'] == c])
    pct  = 100 * len(v) / n_cl if n_cl > 0 else 0
    am_m = float(v[AM_COL].median()) if len(v) > 0 else np.nan
    rv_m = float(v[REVEL_COL].median()) if len(v) > 0 else np.nan
    sc_m = float(v['score_composite'].median()) if len(v) > 0 else np.nan
    print(f'  {c:>3}  {len(v):>7,}  {pct:>12.1f}%  {am_m:>11.3f}  {rv_m:>14.3f}')
    vus_score_rows.append({
        'cluster': c, 'n_vus': len(v), 'pct_vus_in_cluster': round(pct, 1),
        'am_median': am_m, 'revel_median': rv_m, 'composite_median': sc_m,
    })

# Test: i VUS si concentrano nei cluster intermedi?
if best_k >= 3:
    intermediate_cls = [c for c in cluster_ids if c != 0 and c != HR_CLUSTER]
    n_vus_inter = int(vus_sub[vus_sub['cluster'].isin(intermediate_cls)].shape[0])
    n_vus_extr  = len(vus_sub) - n_vus_inter
    pct_inter   = 100 * n_vus_inter / len(vus_sub) if len(vus_sub) > 0 else 0
    print(f'\n  VUS in cluster intermedi ({intermediate_cls}): {n_vus_inter:,}'
          f' / {len(vus_sub):,} ({pct_inter:.1f}%)')
    print(f'  VUS in cluster estremi (C0, C{HR_CLUSTER}):     {n_vus_extr:,}'
          f' ({100-pct_inter:.1f}%)')

    # Fisher: VUS sovra-rappresentati nei cluster intermedi?
    in_inter = df_gmm['cluster'].isin(intermediate_cls).astype(int)
    tab_vus  = pd.crosstab(vus_mask.astype(int), in_inter)
    if tab_vus.shape == (2, 2):
        OR_vus_inter, p_vus_inter = fisher_exact(tab_vus)
        print(f'  OR VUS in intermedi vs estremi: {OR_vus_inter:.2f}  p={p_vus_inter:.2e}')
    else:
        OR_vus_inter, p_vus_inter = np.nan, np.nan

# Kruskal-Wallis score_composite VUS vs classificati
if len(vus_sub) > 10 and len(nvus_sub) > 10:
    _, p_vus_vs = mannwhitneyu(
        vus_sub['score_composite'].dropna(),
        nvus_sub['score_composite'].dropna(), alternative='two-sided')
    print(f'  Mann-Whitney score_composite VUS vs classificati: p={p_vus_vs:.2e}')
    print(f'    VUS mediana={vus_sub["score_composite"].median():.3f}'
          f'  classificati mediana={nvus_sub["score_composite"].median():.3f}')

vus_score_df = pd.DataFrame(vus_score_rows)
vus_score_df.to_csv(f'{OUTDIR}/vus_score_space.csv', index=False)
print('  → vus_score_space.csv salvato')

del vus_sub, nvus_sub
gc.collect()

# =============================================================================
# SEZIONE 8 — ANALISI 5: FREQUENZA ALLELICA (AF_grpmax)
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 8 — FREQUENZA ALLELICA PER CLUSTER')
print('='*70)

df_af = df_gmm[df_gmm[AF_COL].notna()].copy()
print(f'  Varianti con AF reale: {len(df_af):,}')

if len(df_af) > 0:
    if LOG_AF not in df_af.columns:
        df_af[LOG_AF] = np.log10(df_af[AF_COL].clip(lower=1e-6))

    # Kruskal-Wallis
    af_groups = [df_af[df_af['cluster'] == c][AF_COL].values for c in cluster_ids]
    _, p_kw_af = kruskal(*[g for g in af_groups if len(g) > 1])
    rho_af, p_af = spearmanr(df_af['cluster'], df_af[LOG_AF])
    _, p_mw = mannwhitneyu(df_af[df_af['is_hr'] == 1][AF_COL],
                           df_af[df_af['is_hr'] == 0][AF_COL], alternative='less')
    print(f'  Kruskal-Wallis AF  p={p_kw_af:.2e}')
    print(f'  Spearman cluster~log10AF  rho={rho_af:.4f}  p={p_af:.2e}')
    print(f'  Mann-Whitney HR<nonHR     p={p_mw:.2e}')

    # Mediana AF + IEV per cluster
    df_af['is_ultrarare'] = (df_af[AF_COL] < 1e-4).astype('int8')
    df_af['is_iev']       = ((df_af[AF_COL] >= 1e-5) & (df_af[AF_COL] <= 1e-3)).astype('int8')

    print(f'\n  {"Cl":>3}  {"n AF":>7}  {"log10AF med":>12}  {"% no AF":>9}  {"n IEV":>7}  {"% IEV":>7}')
    af_rows = []
    for c in cluster_ids:
        sub_all = df_gmm[df_gmm['cluster'] == c]
        sub_af  = df_af[df_af['cluster'] == c]
        pct_miss = 100 * sub_all[AF_COL].isna().sum() / len(sub_all)
        med      = float(sub_af[LOG_AF].median()) if len(sub_af) > 0 else np.nan
        n_iev    = int(sub_af['is_iev'].sum())
        pct_iev  = 100 * n_iev / len(sub_all)
        tag      = ' ←HR' if c == HR_CLUSTER else ''
        print(f'  {c:>3}  {len(sub_af):>7,}  {med:>12.2f}  {pct_miss:>8.1f}%'
              f'  {n_iev:>7,}  {pct_iev:>6.1f}%{tag}')
        af_rows.append({'cluster': c, 'n_af': len(sub_af), 'log10af_med': med,
                        'pct_missing': pct_miss, 'n_iev': n_iev, 'pct_iev': pct_iev})

    pd.DataFrame(af_rows).to_csv(f'{OUTDIR}/af_by_cluster.csv', index=False)

    # OR ultrarare HR vs rest
    tab_hr = pd.crosstab(df_af['is_hr'], df_af['is_ultrarare'])
    OR_ur, p_ur = (np.nan, np.nan)
    if tab_hr.shape == (2, 2):
        OR_ur, p_ur = fisher_exact(tab_hr)
        print(f'\n  OR ultrarare HR vs nonHR: {OR_ur:.2f}  p={p_ur:.2e}')

    # AF per categoria ClinVar
    af_cat = df_af.groupby('clinvar_category')[AF_COL].agg(['median','mean','count'])
    af_cat.to_csv(f'{OUTDIR}/af_by_clinvar_category.csv')
    print('\n  AF per categoria ClinVar:')
    print(af_cat.round(6).to_string())

    # Cluster intermedi vs estremi (da gradiente)
    path_af   = df_af[df_af['clinvar_category'] == 'pathogenic'][LOG_AF]
    benign_af = df_af[df_af['clinvar_category'] == 'benign'][LOG_AF]
    print(f'\n  Riferimento Pathogenic: log10(AF) med={path_af.median():.2f} (n={len(path_af)})')
    print(f'  Riferimento Benign:     log10(AF) med={benign_af.median():.2f} (n={len(benign_af)})')
    for c in [1, 2]:
        cl_af = df_af[df_af['cluster'] == c][LOG_AF]
        if len(cl_af) > 0 and len(path_af) > 0 and len(benign_af) > 0:
            _, p_vs_p = mannwhitneyu(cl_af, path_af, alternative='two-sided')
            _, p_vs_b = mannwhitneyu(cl_af, benign_af, alternative='two-sided')
            print(f'  C{c} (med={cl_af.median():.2f})  vs Path p={p_vs_p:.2e} | vs Benign p={p_vs_b:.2e}')

    del path_af, benign_af

# =============================================================================
# SEZIONE 9 — ANALISI 6: CONFLICTING PER CLUSTER
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 9 — CONFLICTING PER CLUSTER')
print('='*70)

conf_by_cl = df_gmm.groupby('cluster')['is_conflicting_any'].agg(
    n_conflicting='sum', prop_conflicting='mean'
)
conf_by_cl.to_csv(f'{OUTDIR}/conflicting_by_cluster.csv')
print(conf_by_cl.round(3).to_string())

conf_groups = [df_gmm[df_gmm['cluster'] == c]['is_conflicting_any'].values
               for c in cluster_ids]
_, p_conf_kw = kruskal(*conf_groups)
print(f'\n  Kruskal-Wallis conflicting tra cluster: p={p_conf_kw:.2e}')

# =============================================================================
# SEZIONE 10 — ANALISI 7: DOMAIN ENRICHMENT (Fisher + FDR)
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 10 — DOMAIN ENRICHMENT')
print('='*70)

_DOMAIN_SEARCH = ['interpro', 'pfam', 'domain', 'struct',
                  'is_uniprot', 'is_repeat', 'is_coiled', 'is_disordered',
                  'is_motif', 'is_family', 'is_regulatory', 'is_ptm',
                  'is_binding', 'is_transmembrane', 'is_structural']
# Escludi colonne di testo (note, nomi, strutture) e conteggi non binari
_DOMAIN_EXCLUDE = {'is_structural_domain', 'is_domain_any', 'is_active_site'}
_TEXT_HINTS = ['_note', '_name', '_names', 'subunit', 'structure']
domain_cols_raw = [c for c in df_gmm.columns
                   if any(x in c.lower() for x in _DOMAIN_SEARCH)
                   and c not in _DOMAIN_EXCLUDE
                   and not any(h in c.lower() for h in _TEXT_HINTS)]
# Tieni solo colonne convertibili a numerico con valori in {0,1} o conteggi
domain_cols = []
for col in domain_cols_raw:
    vals = pd.to_numeric(df_gmm[col], errors='coerce').dropna().unique()
    if len(vals) <= 10 and set(vals).issubset({0, 1, 0.0, 1.0}):
        domain_cols.append(col)
    elif col == 'n_interpro_domains':
        # Binarizza: ha almeno 1 dominio InterPro
        df_gmm['has_interpro'] = (pd.to_numeric(df_gmm[col], errors='coerce').fillna(0) > 0).astype('int8')
        domain_cols.append('has_interpro')
print(f'  Domini rilevati (binari): {len(domain_cols)}  —  {domain_cols}')

# Converti a int
for col in domain_cols:
    df_gmm[col] = pd.to_numeric(df_gmm[col], errors='coerce').fillna(0).astype('int8')

enr_rows = []
for col in domain_cols:
    n_tot = int(df_gmm[col].sum())
    if n_tot < 5:
        continue
    tab = pd.crosstab(df_gmm['is_hr'], df_gmm[col])
    if tab.shape != (2, 2):
        continue
    OR, p = fisher_exact(tab)

    OR_pc, p_pc = np.nan, np.nan
    if 'is_pathogenic_clean' in df_gmm.columns:
        t2 = pd.crosstab(df_gmm['is_pathogenic_clean'], df_gmm[col])
        if t2.shape == (2, 2):
            OR_pc, p_pc = fisher_exact(t2)

    OR_cf, p_cf = np.nan, np.nan
    t3 = pd.crosstab(df_gmm['is_conflicting_any'], df_gmm[col])
    if t3.shape == (2, 2):
        OR_cf, p_cf = fisher_exact(t3)

    enr_rows.append({'domain': col, 'OR': OR, 'pvalue': p, 'N_total': n_tot,
                     'OR_path_clean': OR_pc, 'p_path_clean': p_pc,
                     'OR_conflicting': OR_cf, 'p_conflicting': p_cf})

enr_df = pd.DataFrame(enr_rows)
if len(enr_df) > 0:
    enr_df = _bh_safe(enr_df.copy(), 'pvalue').sort_values('OR', ascending=False)
    print(f'\n  {"Domain":<30} {"OR":>8} {"p":>12} {"padj":>10} {"OR_path":>9}')
    for _, r in enr_df.iterrows():
        sig = '*' if r.pvalue < 0.05 else ''
        print(f'  {r.domain:<30} {r.OR:>8.2f} {r.pvalue:>12.2e}'
              f' {r.padj:>10.2e} {r.OR_path_clean:>9.2f} {sig}')
    enr_df.to_csv(f'{OUTDIR}/domain_enrichment.csv', index=False)

    # Gradient Spearman cluster ~ domain
    grad_dom_rows = []
    for dom in domain_cols:
        rho_d, p_d = spearmanr(df_gmm['cluster'], df_gmm[dom])
        grad_dom_rows.append({'domain': dom, 'spearman_rho': rho_d, 'pvalue': p_d})
    grad_dom_df = pd.DataFrame(grad_dom_rows).sort_values('pvalue')
    grad_dom_df.to_csv(f'{OUTDIR}/domain_gradient.csv', index=False)

    # Conflicting per dominio
    conf_dom_rows = []
    for dom in domain_cols:
        sub_dom = df_gmm[df_gmm[dom] == 1]
        if len(sub_dom) < 5:
            continue
        sub_no  = df_gmm[df_gmm[dom] == 0]
        prop_d  = float(sub_dom['is_conflicting_any'].mean())
        prop_no = float(sub_no['is_conflicting_any'].mean())
        tab_c   = pd.crosstab(df_gmm['is_conflicting_any'], df_gmm[dom])
        OR_c2, p_c2 = (np.nan, np.nan)
        if tab_c.shape == (2, 2):
            OR_c2, p_c2 = fisher_exact(tab_c)
        conf_dom_rows.append({'domain': dom, 'n_in_domain': len(sub_dom),
                               'prop_conf_in': prop_d, 'prop_conf_out': prop_no,
                               'OR_conflicting': OR_c2, 'p_conflicting': p_c2})
    pd.DataFrame(conf_dom_rows).sort_values('prop_conf_in', ascending=False
                 ).to_csv(f'{OUTDIR}/conflicting_by_domain.csv', index=False)

gc.collect()

# =============================================================================
# SEZIONE 11 — ANALISI 8: ROBUSTEZZA (Bootstrap n=100, Seed n=10, CV n=20)
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 11 — ROBUSTEZZA GMM (Bootstrap/Seed/CV)')
print('='*70)

# Bootstrap n=100
print('  Bootstrap (n=100)...')
b_ari = []
for i in range(100):
    idx    = resample(range(len(X_scaled)), random_state=i)
    c_boot = GaussianMixture(n_components=best_k, random_state=i).fit_predict(X_scaled[idx])
    b_ari.append(adjusted_rand_score(clusters_original[idx], c_boot))
bm, bs = np.mean(b_ari), np.std(b_ari)
print(f'    ARI mean={bm:.4f}  sd={bs:.4f}')

# Seed stability n=10
print('  Seed stability (n=10)...')
s_ari = []
for s in range(10):
    c_s = GaussianMixture(n_components=best_k, random_state=s).fit_predict(X_scaled)
    s_ari.append(adjusted_rand_score(clusters_original, c_s))
sm, ss = np.mean(s_ari), np.std(s_ari)
print(f'    ARI mean={sm:.4f}  sd={ss:.4f}')

# Cross-validazione n=20
print('  Cross-validation (n=20, test_size=20%)...')
cv_ari = []
for i in range(20):
    tr, te = train_test_split(np.arange(len(X_scaled)), test_size=0.2, random_state=i)
    gmm_cv = GaussianMixture(n_components=best_k, random_state=i)
    gmm_cv.fit(X_scaled[tr])
    c_te = gmm_cv.predict(X_scaled[te])
    cv_ari.append(adjusted_rand_score(clusters_original[te], c_te))
    del gmm_cv
cm, cs = np.mean(cv_ari), np.std(cv_ari)
print(f'    ARI mean={cm:.4f}  sd={cs:.4f}')

# AM-only per confronto
X_am     = X_scaled[:, 0].reshape(-1, 1)
c_am     = GaussianMixture(n_components=best_k, random_state=42).fit_predict(X_am)
sil_am   = silhouette_score(X_am, c_am)
ari_vs   = adjusted_rand_score(clusters_original, c_am)
print(f'\n  AM-only silhouette={sil_am:.4f}  ARI vs 2D={ari_vs:.4f}')

rob_df = pd.DataFrame({'metric': [
    'baseline_silhouette', 'am_only_silhouette', 'am_vs_2d_ari',
    'bootstrap_ari_mean', 'bootstrap_ari_sd',
    'seed_ari_mean', 'seed_ari_sd',
    'cv_ari_mean', 'cv_ari_sd'],
    'value': [sil, sil_am, ari_vs, bm, bs, sm, ss, cm, cs]})
rob_df.to_csv(f'{OUTDIR}/robustness_metrics.csv', index=False)
pd.DataFrame({'iter': range(len(b_ari)), 'bootstrap_ari': b_ari}
             ).to_csv(f'{OUTDIR}/bootstrap_ari.csv', index=False)

del c_am, X_am
gc.collect()

# =============================================================================
# SEZIONE 12 — PREDITTIVITÀ PER GENE (AUC-ROC)
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 12 — PREDITTIVITÀ PER GENE')
print('='*70)

from sklearn.metrics import roc_auc_score

gene_pred_rows = []
for gene in sorted(df_gmm['Gene'].unique()):
    sub = df_gmm[df_gmm['Gene'] == gene]
    cl = sub[sub['clinvar_category'].isin(['pathogenic','likely_pathogenic','benign','likely_benign'])]
    y = cl['is_path_group'].astype(int).values
    if y.sum() < 3 or (1-y).sum() < 3:
        continue
    auc_am = roc_auc_score(y, cl[AM_COL].values)
    auc_rv = roc_auc_score(y, cl[REVEL_COL].values)
    auc_co = roc_auc_score(y, cl['score_composite'].values)

    p_sub = sub[sub['clinvar_category'] == 'pathogenic']
    am_med_p = float(p_sub[AM_COL].median()) if len(p_sub) >= 3 else np.nan
    rv_med_p = float(p_sub[REVEL_COL].median()) if len(p_sub) >= 3 else np.nan

    pct_hr = 100 * (sub['cluster'] == HR_CLUSTER).sum() / len(sub)
    pct_vus = 100 * (sub['clinvar_category'] == 'vus').sum() / len(sub)

    gene_pred_rows.append({
        'gene': gene, 'n_total': len(sub), 'n_classified': len(cl),
        'n_plp': int(y.sum()), 'n_blb': int((1-y).sum()),
        'auc_am': auc_am, 'auc_rv': auc_rv, 'auc_composite': auc_co,
        'am_med_path': am_med_p, 'rv_med_path': rv_med_p,
        'pct_hr': pct_hr, 'pct_vus': pct_vus,
    })
    print(f'  {gene:10} n={len(sub):>5}  cls={len(cl):>4}  P+LP={y.sum():>4}  '
          f'AUC_AM={auc_am:.3f}  AUC_RV={auc_rv:.3f}  AUC_co={auc_co:.3f}  %HR={pct_hr:.1f}')

gene_pred_df = pd.DataFrame(gene_pred_rows)
gene_pred_df.to_csv(f'{OUTDIR}/gene_predictability.csv', index=False)
print(f'  Geni analizzati: {len(gene_pred_df)}')

del gene_pred_rows
gc.collect()

# =============================================================================
# SEZIONE 13 — BANDA IEV (AF 1e-4 — 1e-2)
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 13 — BANDA DI PENETRANZA INTERMEDIA (AF 1e-4 — 1e-2)')
print('='*70)

AF_COL_L = AF_COL  # reminder: 'AF_grpmax'
iev_mask = (df_gmm[AF_COL_L] >= 1e-4) & (df_gmm[AF_COL_L] <= 1e-2)
iev = df_gmm[iev_mask].copy()
ultra = df_gmm[(df_gmm[AF_COL_L] < 1e-4) & df_gmm[AF_COL_L].notna()].copy()

print(f'  IEV (AF 1e-4 — 1e-2): {len(iev):,}')
print(f'  Ultrarare (AF < 1e-4): {len(ultra):,}')
print(f'  Senza AF:              {df_gmm[AF_COL_L].isna().sum():,}')

print('\n  Composizione ClinVar:')
print(f'  {"Categoria":22}  {"IEV":>6} {"% IEV":>7}   {"Ultra":>6} {"% Ultra":>7}')
for cat in ['pathogenic','likely_pathogenic','vus','conflicting','likely_benign','benign']:
    n_i = int((iev['clinvar_category'] == cat).sum())
    p_i = 100 * n_i / len(iev) if len(iev) > 0 else 0
    n_u = int((ultra['clinvar_category'] == cat).sum())
    p_u = 100 * n_u / len(ultra) if len(ultra) > 0 else 0
    print(f'  {cat:22}  {n_i:>6} {p_i:>6.1f}%   {n_u:>6} {p_u:>6.1f}%')

print(f'\n  Score mediani:')
print(f'  {"Fascia":22}  {"AM":>8} {"REVEL":>8} {"AF med":>12}')
for label, sub_f in [('IEV', iev), ('Ultrarare', ultra),
                      ('Comuni (>1e-2)', df_gmm[(df_gmm[AF_COL_L]>1e-2)&df_gmm[AF_COL_L].notna()])]:
    if len(sub_f) == 0:
        continue
    print(f'  {label:22}  {sub_f[AM_COL].median():>8.3f} {sub_f[REVEL_COL].median():>8.3f} '
          f'{sub_f[AF_COL_L].median():>12.2e}')

# IEV per cluster
print(f'\n  IEV per cluster:')
for c in range(best_k):
    n_c = int((iev['cluster'] == c).sum())
    pct_c = 100 * n_c / len(iev) if len(iev) > 0 else 0
    n_plp = int(iev[iev['cluster'] == c]['is_path_group'].sum())
    print(f'    C{c}: {n_c:>5} ({pct_c:.1f}%)  P+LP={n_plp}')

# IEV per gene
iev_gene = iev.groupby('Gene').agg(
    n_iev=('Gene', 'size'),
    am_med=(AM_COL, 'median'), rv_med=(REVEL_COL, 'median'),
    n_path=('is_path_group', 'sum'),
    n_vus=('clinvar_category', lambda x: (x == 'vus').sum()),
    pct_hr=('is_hr', 'mean'),
).sort_values('n_iev', ascending=False)
iev_gene['pct_hr'] = (iev_gene['pct_hr'] * 100).round(1)
iev_gene.to_csv(f'{OUTDIR}/iev_by_gene.csv')
print(f'\n  IEV top 10 geni:')
print(iev_gene.head(10).to_string())

del iev, ultra
gc.collect()

# =============================================================================
# SEZIONE 14 — IDENTIKIT GENE-SPECIFICO + VUS CANDIDATE
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 14 — IDENTIKIT PATOGENICO GENE-SPECIFICO')
print('='*70)

domain_bin_cols = ['is_binding_site', 'is_uniprot_domain', 'is_disordered_region',
                   'is_repeat_region', 'is_coiled_coil', 'is_transmembrane_gene',
                   'is_ptm_site', 'in_interpro']
domain_bin_cols = [c for c in domain_bin_cols if c in df_gmm.columns]

# Build identikits for genes with >= 5 P+LP
identikit = {}
for gene in sorted(df_gmm['Gene'].unique()):
    sub = df_gmm[df_gmm['Gene'] == gene]
    path = sub[sub['clinvar_category'].isin(['pathogenic', 'likely_pathogenic'])]
    if len(path) < 5:
        continue

    prof = {
        'gene': gene, 'n_plp': len(path),
        'am_lo': float(path[AM_COL].quantile(0.10)),
        'am_hi': float(path[AM_COL].quantile(0.90)),
        'am_med': float(path[AM_COL].median()),
        'rv_lo': float(path[REVEL_COL].quantile(0.10)),
        'rv_hi': float(path[REVEL_COL].quantile(0.90)),
        'rv_med': float(path[REVEL_COL].median()),
        'pct_hr_path': 100 * (path['cluster'] == HR_CLUSTER).sum() / len(path),
        'af_q90_path': float(path[AF_COL].quantile(0.90)) if path[AF_COL].notna().sum() >= 3 else np.nan,
    }
    # Domain enrichment in pathogenics vs gene total
    for d in domain_bin_cols:
        pct_p = float(path[d].mean())
        pct_a = float(sub[d].mean())
        prof[f'{d}_enriched'] = pct_p > pct_a + 0.05

    identikit[gene] = prof
    enr = [d for d in domain_bin_cols if prof.get(f'{d}_enriched', False)]
    print(f'  {gene:10} P+LP={len(path):>4}  AM=[{prof["am_lo"]:.3f}–{prof["am_hi"]:.3f}]  '
          f'RV=[{prof["rv_lo"]:.3f}–{prof["rv_hi"]:.3f}]  %HR={prof["pct_hr_path"]:.0f}%'
          + (f'  enr={enr}' if enr else ''))

print(f'\n  Identikit costruiti: {len(identikit)} geni')
pd.DataFrame(identikit.values()).to_csv(f'{OUTDIR}/gene_identikits.csv', index=False)

# Score VUS against gene-specific identikit
print('\n  Scoring VUS contro identikit...')
vus_all = df_gmm[df_gmm['clinvar_category'] == 'vus'].copy()
all_candidates = []

for gene, prof in identikit.items():
    sub_vus = vus_all[vus_all['Gene'] == gene].copy()
    if len(sub_vus) == 0:
        continue

    score = pd.Series(0.0, index=sub_vus.index)
    max_score = 0

    # 1. AM in pathogenic range
    max_score += 1
    score += ((sub_vus[AM_COL] >= prof['am_lo']) & (sub_vus[AM_COL] <= prof['am_hi'])).astype(float)

    # 2. REVEL in pathogenic range
    max_score += 1
    score += ((sub_vus[REVEL_COL] >= prof['rv_lo']) & (sub_vus[REVEL_COL] <= prof['rv_hi'])).astype(float)

    # 3. Cluster HR (if >50% of pathogenics are HR)
    if prof['pct_hr_path'] > 50:
        max_score += 1
        score += (sub_vus['cluster'] == HR_CLUSTER).astype(float)

    # 4. AF ≤ q90 of pathogenics (or no AF)
    if not np.isnan(prof['af_q90_path']):
        max_score += 1
        score += (sub_vus[AF_COL].isna() | (sub_vus[AF_COL] <= prof['af_q90_path'])).astype(float)

    # 5. Enriched domain criteria
    for d in domain_bin_cols:
        if prof.get(f'{d}_enriched', False):
            max_score += 1
            score += sub_vus[d].astype(float)

    sub_vus['identikit_score'] = score
    sub_vus['identikit_max'] = max_score
    sub_vus['identikit_pct'] = (score / max_score * 100).round(1)

    threshold = max(0.80 * max_score, max_score - 1)
    cands = sub_vus[sub_vus['identikit_score'] >= threshold].copy()
    cands['gene_profile'] = gene

    if len(cands) > 0:
        all_candidates.append(cands)

if all_candidates:
    cand_df = pd.concat(all_candidates, ignore_index=True)
    n_full = int((cand_df['identikit_score'] == cand_df['identikit_max']).sum())

    out_cols = ['Gene', 'Name', 'gene_profile', AM_COL, REVEL_COL, 'score_composite',
                AF_COL, 'cluster', 'is_hr', 'clinvar_category',
                'identikit_score', 'identikit_max', 'identikit_pct',
                'aa_pos', 'Protein_Length', 'aa_pos_norm', 'ClinicalSignificance']
    out_cols += [d for d in domain_bin_cols if d in cand_df.columns]
    out_cols = [c for c in out_cols if c in cand_df.columns]

    cand_out = cand_df[out_cols].sort_values(['Gene', 'identikit_pct'], ascending=[True, False])
    cand_out.to_csv(f'{OUTDIR}/vus_reclassification_candidates.csv', index=False)

    # Confronto con soglia globale
    global_thr = vus_all[(vus_all[AM_COL] >= 0.564) & (vus_all[REVEL_COL] >= 0.75) &
                          (vus_all[AF_COL].isna() | (vus_all[AF_COL] < 1e-4))]
    overlap = len(set(cand_df.index) & set(global_thr.index))
    only_idk = len(cand_df) - overlap
    only_glob = len(global_thr) - overlap

    print(f'\n  VUS candidate (identikit):  {len(cand_df):,}  (match perfetto: {n_full:,})')
    print(f'  VUS soglia globale:         {len(global_thr):,}')
    print(f'  Overlap:                    {overlap:,}')
    print(f'  Solo identikit (rescued):   {only_idk:,}')
    print(f'  Solo globale:               {only_glob:,}')

    n_idk_total = len(cand_df)
    n_idk_full  = n_full

    print(f'\n  Per gene:')
    for gene in sorted(cand_df['gene_profile'].unique()):
        sub_c = cand_df[cand_df['gene_profile'] == gene]
        n_f = int((sub_c['identikit_score'] == sub_c['identikit_max']).sum())
        n_vus_gene = int((vus_all['Gene'] == gene).sum())
        print(f'    {gene:10}: {len(sub_c):>5} / {n_vus_gene:>5} VUS  '
              f'(match perfetto: {n_f})')
else:
    n_idk_total = 0
    n_idk_full  = 0
    print('  Nessuna VUS candidata trovata')

del vus_all, all_candidates
gc.collect()

# =============================================================================
# SEZIONE 15 — SALVA METRICHE MASTER
# =============================================================================

master_metrics = pd.DataFrame({'metric': [
    'n_variants_total', 'n_gmm', 'n_propagated', 'best_k', 'HR_cluster',
    'silhouette', 'davies_bouldin',
    'spearman_AM_REVEL_rho', 'spearman_AM_REVEL_p',
    'gradient_spearman_rho', 'gradient_spearman_p',
    'OR_path_HR_vs_C0', 'p_path_HR_vs_C0',
    'OR_pathgroup_HR_vs_C0', 'p_pathgroup_HR_vs_C0',
    'kruskal_AM_H', 'kruskal_AM_p', 'kruskal_REVEL_H', 'kruskal_REVEL_p',
    'bootstrap_ari_mean', 'bootstrap_ari_sd',
    'seed_ari_mean', 'seed_ari_sd',
    'cv_ari_mean', 'cv_ari_sd',
    'spearman_rank_AM_rho', 'spearman_rank_AM_p',
    'spearman_rank_REVEL_rho', 'spearman_rank_REVEL_p',
    'monotonicity_composite_rho', 'monotonicity_composite_p',
    'monotonicity_AM_rho', 'monotonicity_REVEL_rho',
    'n_genes_with_auc', 'n_identikit_genes',
    'n_vus_identikit_candidates', 'n_vus_identikit_full_match',
],
'value': [
    len(df_full), len(df_gmm), len(df_full)-len(df_gmm), best_k, HR_CLUSTER,
    sil, db, rho_feat, p_feat,
    rho_g, p_g, OR_p, p_p, OR_pg, p_pg,
    H_am, p_kw_am, H_rv, p_kw_rv,
    bm, bs, sm, ss, cm, cs,
    rho_am, p_am, rho_rv, p_rv,
    rho_mono_comp, p_mono_comp, rho_mono_am, rho_mono_rv,
    len(gene_pred_df), len(identikit),
    n_idk_total, n_idk_full,
]})
master_metrics.to_csv(f'{OUTDIR}/master_metrics.csv', index=False)

# =============================================================================
# SEZIONE 16 — FIGURE 1–6
# =============================================================================

print('\n' + '='*70)
print('SEZIONE 16 — FIGURE')
print('='*70)

# Subset ridotto per scatter (risparmio RAM)
_samp = df_gmm.dropna(subset=[AM_COL, REVEL_COL]).sample(
    min(5000, len(df_gmm)), random_state=42)

# ---- FIGURA 1 — Boxplot gradiente AM + REVEL per classe ----
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, score, lbl in zip(axes, [REVEL_COL, AM_COL], ['REVEL', 'AlphaMissense']):
    plot_d = df_gmm[df_gmm['class_clean'].isin(ORDINE_CLASSI)].copy()
    plot_d['class_clean'] = pd.Categorical(
        plot_d['class_clean'], categories=ORDINE_CLASSI, ordered=True)
    sns.boxplot(data=plot_d.sort_values('class_clean'),
                x='class_clean', y=score, order=ORDINE_CLASSI,
                hue='class_clean', palette=PALETTE_CLASS,
                showfliers=False, width=0.6, linewidth=1.2, legend=False, ax=ax)
    ax.set_xlabel('')
    ax.set_ylabel(lbl, fontsize=12)
    ax.set_title(f'Gradiente — {lbl}', fontsize=13)
    ax.set_xticks(range(len(ORDINE_CLASSI)))
    ax.set_xticklabels([c.replace('/', '/\n') for c in ORDINE_CLASSI],
                        rotation=40, ha='right', fontsize=9)
    del plot_d
plt.suptitle('Distribuzione score per classe ClinVar (categoria esclusiva)', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/fig1_gradiente_scores.png', dpi=300, bbox_inches='tight')
plt.close(); gc.collect()
print('  fig1_gradiente_scores.png')

# ---- FIGURA 2 — Scatter AM × REVEL per cluster ----
cluster_labels_map = {
    0: 'C0 — Basso rischio',
    1: 'C1 — Intermedio REVEL',
    2: 'C2 — Intermedio AM',
    3: 'C3 — Alto rischio',
}
fig, ax = plt.subplots(figsize=(8, 7))
for c in cluster_ids:
    sub = _samp[_samp['cluster'] == c]
    ax.scatter(sub[REVEL_COL], sub[AM_COL],
               c=PALETTE_CLUSTER.get(c, '#999'), alpha=0.35, s=8,
               label=cluster_labels_map.get(c, f'Cluster {c}'))
ax.axhline(0.564, color='gray', linestyle='--', linewidth=0.8)
ax.axvline(0.75,  color='gray', linestyle='--', linewidth=0.8)
ax.text(0.76, 0.02, 'soglia REVEL', fontsize=8, color='gray')
ax.text(0.02, 0.575, 'soglia AM',   fontsize=8, color='gray')
ax.set_xlabel('REVEL', fontsize=12); ax.set_ylabel('AlphaMissense', fontsize=12)
ax.set_title('Spazio AM × REVEL — cluster GMM k=4', fontsize=13)
ax.legend(markerscale=3, fontsize=10)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/fig2_cluster_score_space.png', dpi=300, bbox_inches='tight')
plt.close(); gc.collect()
print('  fig2_cluster_score_space.png')

# ---- FIGURA 3 — Scatter AF × REVEL per classe ----
plot_classes = ['Pathogenic','Likely_pathogenic','VUS','Likely_benign','Benign']
plot_af_fig  = df_gmm[df_gmm['class_clean'].isin(plot_classes)].dropna(
    subset=[AF_COL, REVEL_COL]).copy()
if len(plot_af_fig) > 5000:
    plot_af_fig = plot_af_fig.sample(5000, random_state=42)
fig, ax = plt.subplots(figsize=(9, 6))
for cat in plot_classes:
    sub = plot_af_fig[plot_af_fig['class_clean'] == cat]
    if len(sub) == 0:
        continue
    ax.scatter(sub[AF_COL], sub[REVEL_COL],
               c=PALETTE_CLASS.get(cat, '#999'), alpha=0.35, s=8, label=cat)
ax.set_xscale('log')
ax.set_xlabel('AF_grpmax (log)', fontsize=12); ax.set_ylabel('REVEL', fontsize=12)
ax.set_title('Spazio AF × REVEL per classe ClinVar', fontsize=13)
ax.legend(markerscale=3, fontsize=10)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/fig3_af_revel.png', dpi=300, bbox_inches='tight')
plt.close()
del plot_af_fig; gc.collect()
print('  fig3_af_revel.png')

# ---- FIGURA 4 — Heatmap P(Pathogenic) nello spazio AF × REVEL ----
df_hm = df_gmm.dropna(subset=[AF_COL, REVEL_COL]).copy()
df_hm = df_hm[df_hm[AF_COL] > 0]
df_hm['is_path_hm'] = (df_hm['clinvar_category'] == 'pathogenic').astype(int)
df_hm['af_bin_hm']  = pd.cut(np.log10(df_hm[AF_COL].clip(lower=1e-7)), bins=25)
df_hm['rv_bin_hm']  = pd.cut(df_hm[REVEL_COL], bins=25)
grid = (df_hm.groupby(['af_bin_hm', 'rv_bin_hm'], observed=True)['is_path_hm']
         .agg(['mean', 'count']).reset_index())
grid = grid[grid['count'] >= 5]
grid_piv = grid.pivot(index='rv_bin_hm', columns='af_bin_hm', values='mean'
                      ).sort_index(ascending=False)
print(f'  Heatmap celle con dati (≥5 var): {grid_piv.notna().sum().sum()}')
print(f'  P(path) max nel grid: {grid_piv.max().max():.3f}')
fig, ax = plt.subplots(figsize=(12, 7))
sns.heatmap(grid_piv, cmap='viridis', vmin=0, vmax=grid_piv.max().max(),
            cbar_kws={'label': 'P(Pathogenic | ClinVar consenso)'},
            ax=ax, xticklabels=5, yticklabels=5)
ax.set_title('Probability landscape of pathogenicity (AF × REVEL, ≥5 obs/cella)',
             fontsize=12)
ax.set_xlabel('log₁₀(AF_grpmax)', fontsize=11)
ax.set_ylabel('REVEL', fontsize=11)
ax.tick_params(axis='x', rotation=45, labelsize=7)
ax.tick_params(axis='y', rotation=0,  labelsize=7)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/fig4_prob_pathogenicity_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
del df_hm, grid, grid_piv; gc.collect()
print('  fig4_prob_pathogenicity_heatmap.png')

# ---- FIGURA 5 — Domain enrichment (volcano + heatmap cluster×domain) ----
if len(enr_df) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Volcano
    ax = axes[0]
    colors_v = ['crimson' if p < 0.05 else 'grey' for p in enr_df['pvalue']]
    ax.scatter(enr_df['OR'], -np.log10(enr_df['pvalue'] + 1e-10),
               c=colors_v, alpha=0.7, s=40)
    ax.axvline(1, linestyle='--', color='black', linewidth=0.8)
    ax.axhline(-np.log10(0.05), linestyle='--', color='red', linewidth=0.8)
    for _, row in enr_df[enr_df['pvalue'] < 0.05].iterrows():
        ax.annotate(row['domain'][:20],
                    (row['OR'], -np.log10(row['pvalue'] + 1e-10)), fontsize=6, alpha=0.8,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.3))
    ax.set_xlabel('Odds Ratio'); ax.set_ylabel('-log10(p)')
    ax.set_title('Domain Enrichment in HR Cluster'); ax.grid(alpha=0.3)

    # Cluster-domain heatmap
    ax = axes[1]
    cd_rows = []
    for c in cluster_ids:
        sub = df_gmm[df_gmm['cluster'] == c]
        for dom in domain_cols:
            cd_rows.append({'cluster': c, 'domain': dom,
                             'proportion': float(sub[dom].mean())})
    cd_df  = pd.DataFrame(cd_rows)
    pivot  = cd_df.pivot(index='cluster', columns='domain', values='proportion')
    sns.heatmap(pivot, cmap='viridis', ax=ax, annot=True, fmt='.2f',
                linewidths=0.3, cbar_kws={'label': 'Proportion'})
    ax.set_title('Domain Burden per Cluster')
    ax.tick_params(axis='x', rotation=45)
    del cd_rows, cd_df

    plt.suptitle('Domain Enrichment Analysis', fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig5_domain_enrichment.png', dpi=300, bbox_inches='tight')
    plt.close()
    print('  fig5_domain_enrichment.png')

# ---- FIGURA 6 — Bootstrap distribution ----
fig, ax = plt.subplots(figsize=(8, 5))
sns.histplot(b_ari, bins=15, kde=True, alpha=0.7, label=f'Bootstrap n=100  mean={bm:.3f}', ax=ax)
ax.axvline(bm, color='steelblue', linestyle='--')
sns.histplot(cv_ari, bins=15, kde=True, alpha=0.5, color='orange',
             label=f'CV n=20  mean={cm:.3f}', ax=ax)
ax.axvline(cm, color='orange', linestyle='--')
ax.set_xlabel('Adjusted Rand Index'); ax.set_title('Robustezza GMM k=4')
ax.legend()
plt.tight_layout()
plt.savefig(f'{OUTDIR}/fig6_robustness.png', dpi=300, bbox_inches='tight')
plt.close()
print('  fig6_robustness.png')

del _samp; gc.collect()

# ---- FIGURA 7 — Probability Surface AM × REVEL ----
print('\n  Generating fig7: Probability Surface AM × REVEL...')
_classified = df_gmm[df_gmm['clinvar_category'].isin(
    ['pathogenic','likely_pathogenic','benign','likely_benign'])].copy()
_classified['_y'] = _classified['clinvar_category'].isin(
    ['pathogenic','likely_pathogenic']).astype(int)

_N_BINS = 20
_am_bins = np.linspace(0, 1, _N_BINS+1)
_rv_bins = np.linspace(0, 1, _N_BINS+1)
_classified['_am_bin'] = pd.cut(_classified[AM_COL], bins=_am_bins, include_lowest=True)
_classified['_rv_bin'] = pd.cut(_classified[REVEL_COL], bins=_rv_bins, include_lowest=True)
_grid = _classified.groupby(['_am_bin','_rv_bin'], observed=True).agg(
    n_path=('_y','sum'), n_total=('_y','size')
).reset_index()
_grid['p_path'] = _grid['n_path'] / _grid['n_total']
_grid = _grid[_grid['n_total'] >= 3]
_pivot7 = _grid.pivot(index='_am_bin', columns='_rv_bin', values='p_path').sort_index(ascending=False)

fig, ax = plt.subplots(figsize=(10, 9))
sns.heatmap(_pivot7, cmap='RdYlBu_r', vmin=0, vmax=1,
            cbar_kws={'label': 'P(Pathogenic | ClinVar classificato)', 'shrink': 0.8},
            ax=ax, xticklabels=4, yticklabels=4, linewidths=0.1, linecolor='white')
ax.set_xlabel('REVEL', fontsize=12); ax.set_ylabel('AlphaMissense', fontsize=12)
ax.set_title('Probability landscape: P(Pathogenic) nello spazio AM × REVEL\n'
             '(solo varianti classificate P/LP/B/LB, ≥3 obs/cella)', fontsize=12)
ax.tick_params(axis='x', rotation=45, labelsize=7)
ax.tick_params(axis='y', rotation=0, labelsize=7)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/fig7_probability_surface_AM_REVEL.png', dpi=300, bbox_inches='tight')
plt.close(); gc.collect()
print('  fig7_probability_surface_AM_REVEL.png')

# ---- FIGURA 8 — Probability Surface stratificata per AF ----
print('  Generating fig8: Probability Surface × AF...')
_class_af = _classified[_classified[AF_COL].notna()].copy()
_af_bands = [
    ('Ultrarare (AF < 10⁻⁵)', _class_af[_class_af[AF_COL] < 1e-5]),
    ('Rare (10⁻⁵ ≤ AF < 10⁻³)', _class_af[(_class_af[AF_COL] >= 1e-5) & (_class_af[AF_COL] < 1e-3)]),
    ('Comuni (AF ≥ 10⁻³)', _class_af[_class_af[AF_COL] >= 1e-3]),
]
fig, axes = plt.subplots(1, 3, figsize=(22, 7))
for idx, (label, sub_af) in enumerate(_af_bands):
    ax = axes[idx]
    if len(sub_af) < 10:
        ax.set_title(f'{label}\n(n={len(sub_af)}, insufficiente)')
        continue
    sub_af = sub_af.copy()
    sub_af['_am_bin'] = pd.cut(sub_af[AM_COL], bins=_am_bins, include_lowest=True)
    sub_af['_rv_bin'] = pd.cut(sub_af[REVEL_COL], bins=_rv_bins, include_lowest=True)
    _gr = sub_af.groupby(['_am_bin','_rv_bin'], observed=True).agg(
        n_path=('_y','sum'), n_total=('_y','size')
    ).reset_index()
    _gr['p_path'] = _gr['n_path'] / _gr['n_total']
    _gr = _gr[_gr['n_total'] >= 3]
    _piv = _gr.pivot(index='_am_bin', columns='_rv_bin', values='p_path').sort_index(ascending=False)
    n_p = int(sub_af['_y'].sum()); n_b = len(sub_af) - n_p
    sns.heatmap(_piv, cmap='RdYlBu_r', vmin=0, vmax=1,
                cbar_kws={'label': 'P(Path)', 'shrink': 0.7},
                ax=ax, xticklabels=4, yticklabels=4, linewidths=0.1, linecolor='white')
    ax.set_title(f'{label}\nn={len(sub_af):,} (P+LP={n_p}, B+LB={n_b})', fontsize=10)
    ax.set_xlabel('REVEL', fontsize=10)
    ax.set_ylabel('AlphaMissense' if idx == 0 else '', fontsize=10)
    ax.tick_params(axis='x', rotation=45, labelsize=6)
    ax.tick_params(axis='y', rotation=0, labelsize=6)
plt.suptitle('Probability landscape stratificato per frequenza allelica', fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/fig8_probability_surface_by_AF.png', dpi=300, bbox_inches='tight')
plt.close()
del _classified, _class_af, _grid, _pivot7; gc.collect()
print('  fig8_probability_surface_by_AF.png')

# ---- ANALISI: Enrichment classe funzionale × cluster ----
print('\n  Enrichment classe funzionale per cluster...')
_SARC = {'ACTC1','MYBPC3','MYH7','MYL2','MYL3','TNNC1','TNNI3','TNNT2','TPM1'}
_DESM = {'DSG2','DSP','PKP2'}
_CYTO = {'DES','FLNC','LMNA','VCL','NEXN'}
_ION  = {'SCN5A'}
_RNA  = {'RBM20'}

def _fc(gene):
    if gene in _SARC: return 'Sarcomerica'
    elif gene in _DESM: return 'Desmosomiale'
    elif gene in _CYTO: return 'Citoscheletrica'
    elif gene in _ION: return 'Canale ionico'
    elif gene in _RNA: return 'RNA-binding'
    else: return 'Altra'

df_gmm['func_class'] = df_gmm['Gene'].apply(_fc)
_enr_rows = []
for cls in ['Sarcomerica','Desmosomiale','Citoscheletrica','Canale ionico','RNA-binding','Altra']:
    sub_fc = df_gmm[df_gmm['func_class'] == cls]
    if len(sub_fc) < 20: continue
    pct_hr = 100 * (sub_fc['cluster'] == HR_CLUSTER).sum() / len(sub_fc)
    _is_cls = (df_gmm['func_class'] == cls).astype(int)
    _is_hr  = (df_gmm['cluster'] == HR_CLUSTER).astype(int)
    _tab = pd.crosstab(_is_cls, _is_hr)
    _or, _p = fisher_exact(_tab) if _tab.shape == (2,2) else (np.nan, np.nan)
    _enr_rows.append({'class': cls, 'n': len(sub_fc), 'pct_hr': pct_hr, 'OR': _or, 'p': _p})
    print(f'    {cls:20} n={len(sub_fc):>6}  %HR={pct_hr:.1f}%  OR={_or:.2f}  p={_p:.2e}')
pd.DataFrame(_enr_rows).to_csv(f'{OUTDIR}/functional_class_enrichment.csv', index=False)

# C1 vs C2 domini
_c1c2_rows = []
_c1 = df_gmm[df_gmm['cluster'] == 1]
_c2 = df_gmm[df_gmm['cluster'] == 2]
for d in domain_cols:
    _tab = pd.crosstab(df_gmm[df_gmm['cluster'].isin([1,2])]['cluster'],
                        df_gmm[df_gmm['cluster'].isin([1,2])][d])
    _or_d, _p_d = fisher_exact(_tab) if _tab.shape == (2,2) else (np.nan, np.nan)
    _c1c2_rows.append({'domain': d, 'pct_c1': 100*_c1[d].mean(),
                        'pct_c2': 100*_c2[d].mean(), 'OR_c1_vs_c2': _or_d, 'p': _p_d})
pd.DataFrame(_c1c2_rows).to_csv(f'{OUTDIR}/c1_vs_c2_domains.csv', index=False)
print('  c1_vs_c2_domains.csv salvato')

gc.collect()

# =============================================================================
# RIEPILOGO FINALE
# =============================================================================

print('\n' + '='*70)
print('PIPELINE COMPLETATA')
print('='*70)
print(f'  Output directory  : {OUTDIR}/')
print(f'  Varianti totali   : {len(df_full):,}')
print(f'  Subset GMM        : {len(df_gmm):,}')
print(f'  Cluster k         : {best_k}  (HR = Cluster {HR_CLUSTER})')
print(f'  Silhouette        : {sil:.4f}')
print(f'  Davies-Bouldin    : {db:.4f}')
print(f'  Bootstrap ARI     : {bm:.4f} ± {bs:.4f}')
print(f'  Seed ARI          : {sm:.4f} ± {ss:.4f}')
print(f'  CV ARI            : {cm:.4f} ± {cs:.4f}')
print(f'  OR Fisher HR vs C0: {OR_p:.1f} (P puro) / {OR_pg:.1f} (P+LP)')
print(f'  ClinVar fix       : is_conflicting priorità massima (n_both={n_both:,})')
print()
print('  File salvati:')
for f in sorted(os.listdir(OUTDIR)):
    size = os.path.getsize(f'{OUTDIR}/{f}')
    print(f'    {f:<45} {size/1024:>6.0f} KB')
