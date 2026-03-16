#!/usr/bin/env python3
# =============================================================================
# VALIDAZIONE IDENTIKIT GENE-SPECIFICO — Multi-anno v1.0
#
# Per ogni release ClinVar storica:
#   1. Trova le VUS di quell'anno
#   2. Scorale contro l'identikit del proprio gene
#   3. Verifica se le riclassificate P/LP matchano di più
#   4. Confronta identikit vs soglia globale vs cluster HR
#   5. Breakdown per gene
#
# INPUT:
#   - PIPELINE_RESULTS/clustered_dataset.csv.gz
#   - PIPELINE_RESULTS/gene_identikits.csv
#   - variant_summary_YYYY*.txt.gz (tutte le release disponibili)
#
# OUTPUT (in IDENTIKIT_VALIDATION/):
#   - identikit_temporal_summary.csv       (metriche per anno)
#   - identikit_per_gene_YYYY.csv          (per gene, per anno)
#   - identikit_method_comparison.csv      (identikit vs globale vs HR)
#   - identikit_pool_analysis.csv          (pool tutte le riclassificate)
#   - fig_identikit_validation.png
# =============================================================================

import os, sys, re, glob, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
import warnings
warnings.filterwarnings('ignore')

try:
    sys.stdout.reconfigure(encoding='utf-8')
except:
    pass

sns.set_theme(style='whitegrid', font_scale=1.0)

OUTDIR = 'IDENTIKIT_VALIDATION'
os.makedirs(OUTDIR, exist_ok=True)

print('=' * 70)
print('VALIDAZIONE IDENTIKIT GENE-SPECIFICO — Multi-anno v1.0')
print('=' * 70)

# =============================================================================
# 1. CARICAMENTO PIPELINE + IDENTIKIT
# =============================================================================

print('\n[1/6] Caricamento dataset e identikit...')

# Pipeline dataset
PIPELINE_FILE = None
for c in ['PIPELINE_RESULTS/clustered_dataset.csv.gz',
          'PIPELINE_RESULTS/clustered_dataset.xlsx',
          'FINAL_HYBRID_PIPELINE_RESULTS/clustered_dataset.xlsx',
          'clustered_dataset.csv.gz']:
    if os.path.exists(c):
        PIPELINE_FILE = c; break

if PIPELINE_FILE is None:
    print('  ERRORE: clustered_dataset non trovato'); sys.exit(1)

if PIPELINE_FILE.endswith('.gz'):
    df = pd.read_csv(PIPELINE_FILE, low_memory=False)
elif PIPELINE_FILE.endswith('.xlsx'):
    df = pd.read_excel(PIPELINE_FILE)
else:
    df = pd.read_csv(PIPELINE_FILE, low_memory=False)

HR = int(df['cluster'].max())
df['is_hr'] = (df['cluster'] == HR).astype(int)
df['VariationID'] = pd.to_numeric(df.get('VariationID'), errors='coerce')
df['AlleleID'] = pd.to_numeric(df.get('AlleleID'), errors='coerce')
gene_col = 'Gene' if 'Gene' in df.columns else 'GeneSymbol'
AM = 'am_pathogenicity'; RV = 'REVEL_max'; AF = 'AF_grpmax'

def current_class(row):
    cat = row.get('clinvar_category', '')
    if pd.notna(cat) and str(cat).strip(): return str(cat)
    if row.get('is_pathogenic', 0) == 1: return 'pathogenic'
    if row.get('is_likely_pathogenic', 0) == 1: return 'likely_pathogenic'
    if row.get('is_vus', 0) == 1: return 'vus'
    if row.get('is_likely_benign', 0) == 1: return 'likely_benign'
    if row.get('is_benign', 0) == 1: return 'benign'
    return 'other'
df['class_current'] = df.apply(current_class, axis=1)

print(f'  Dataset: {len(df):,} varianti | HR={HR}')

# Identikit
IDK_FILE = None
for c in ['PIPELINE_RESULTS/gene_identikits.csv', 'gene_identikits.csv']:
    if os.path.exists(c):
        IDK_FILE = c; break

if IDK_FILE is None:
    print('  ERRORE: gene_identikits.csv non trovato'); sys.exit(1)

idf = pd.read_csv(IDK_FILE)
domain_cols = [c for c in idf.columns if c.endswith('_enriched')]

# Build identikit lookup
idk = {}
for _, row in idf.iterrows():
    g = row['gene']
    enr = [c.replace('_enriched', '') for c in domain_cols
           if str(row.get(c)) == 'True' or row.get(c) == True]

    n_crit = 2  # AM + REVEL always
    af_max = pd.to_numeric(row.get('af_q90_path'), errors='coerce')
    if pd.notna(af_max): n_crit += 1
    hr_req = row.get('pct_hr_path', 0) >= 50
    if hr_req: n_crit += 1
    if enr: n_crit += 1

    idk[g] = {
        'am_lo': row['am_lo'], 'am_hi': row['am_hi'],
        'rv_lo': row['rv_lo'], 'rv_hi': row['rv_hi'],
        'af_max': af_max, 'hr_constrained': hr_req,
        'enriched_domains': enr, 'n_criteria': n_crit,
    }

print(f'  Identikit: {len(idk)} geni')
for g, k in idk.items():
    print(f'    {g:10} AM=[{k["am_lo"]:.2f}-{k["am_hi"]:.2f}]  '
          f'RV=[{k["rv_lo"]:.2f}-{k["rv_hi"]:.2f}]  '
          f'HR={"si" if k["hr_constrained"] else "no"}  '
          f'domini={k["enriched_domains"]}  criteri={k["n_criteria"]}')

# =============================================================================
# 2. SCORING FUNCTION
# =============================================================================

def score_vus(row, identikit_dict, gene_column, hr_cluster):
    g = row.get(gene_column, '')
    if g not in identikit_dict:
        return np.nan, np.nan, np.nan
    k = identikit_dict[g]
    hits = 0
    total = k['n_criteria']
    details = []

    # AM range
    am = row.get(AM, np.nan)
    if pd.notna(am) and k['am_lo'] <= am <= k['am_hi']:
        hits += 1; details.append('AM')

    # REVEL range
    rv = row.get(RV, np.nan)
    if pd.notna(rv) and k['rv_lo'] <= rv <= k['rv_hi']:
        hits += 1; details.append('RV')

    # AF
    if pd.notna(k['af_max']):
        af = row.get(AF, np.nan)
        if pd.notna(af) and af <= k['af_max']:
            hits += 1; details.append('AF')
        elif pd.isna(af):
            hits += 1; details.append('AF_na')

    # HR cluster
    if k['hr_constrained']:
        if row.get('cluster') == hr_cluster:
            hits += 1; details.append('HR')

    # Enriched domains
    if k['enriched_domains']:
        for dom in k['enriched_domains']:
            if dom in row.index and row.get(dom, 0) == 1:
                hits += 1; details.append(dom); break

    pct = hits / total if total > 0 else 0
    return hits, pct, ';'.join(details)

# =============================================================================
# 3. DETECTION FILE CLINVAR STORICI
# =============================================================================

print('\n[2/6] Ricerca file ClinVar storici...')

def classify_clinvar(sig):
    if pd.isna(sig): return 'other'
    s = str(sig).lower()
    if 'conflict' in s: return 'conflicting'
    if 'pathogenic' in s and 'likely' not in s and 'benign' not in s: return 'pathogenic'
    if 'likely pathogenic' in s: return 'likely_pathogenic'
    if 'benign' in s and 'likely' not in s and 'pathogenic' not in s: return 'benign'
    if 'likely benign' in s: return 'likely_benign'
    if 'uncertain' in s: return 'vus'
    return 'other'

def find_clinvar_files():
    patterns = ['variant_summary_*.txt.gz', 'variant_summary_*.txt',
                'variant_summary*.txt.gz', 'variant_summary*.txt']
    found = {}
    for pattern in patterns:
        for fp in glob.glob(pattern):
            match = re.search(r'(20\d{2})', fp)
            if match:
                year = int(match.group(1))
                if year not in found: found[year] = fp
    return dict(sorted(found.items()))

clinvar_files = find_clinvar_files()
if not clinvar_files:
    print('  ERRORE: nessun variant_summary trovato'); sys.exit(1)

for year, path in clinvar_files.items():
    print(f'  {year}: {path} ({os.path.getsize(path)/1e6:.0f} MB)')

# =============================================================================
# 4. LOOP ANNUALE
# =============================================================================

print('\n[3/6] Analisi per anno...')

GENES = set(df[gene_col].dropna().unique())
all_yearly = []
all_yearly_genes = []
all_method_comp = []
all_reclass_pool = []

for year, filepath in clinvar_files.items():
    print(f'\n  === {year} ===')
    t0 = time.time()

    try:
        cv = pd.read_csv(filepath, sep='\t',
                          compression='gzip' if filepath.endswith('.gz') else None,
                          low_memory=False, on_bad_lines='skip')
    except TypeError:
        cv = pd.read_csv(filepath, sep='\t', low_memory=False, error_bad_lines=False)

    if cv.columns[0].startswith('#'):
        cv = cv.rename(columns={cv.columns[0]: cv.columns[0].lstrip('#')})

    # Find ID column (VariantID or AlleleID)
    vid_col = None; merge_key = 'VariationID'
    for c in cv.columns:
        if 'variantid' in c.lower().replace(' ', '').replace('#', ''):
            vid_col = c; merge_key = 'VariationID'; break
    if vid_col is None:
        for c in cv.columns:
            if 'alleleid' in c.lower().replace(' ', '').replace('#', ''):
                vid_col = c; merge_key = 'AlleleID'; break
    if vid_col is None:
        print(f'    SKIP: no ID column'); continue

    # Gene filter
    gene_cv = None
    for c in cv.columns:
        if c.lower() in ('genesymbol', 'gene symbol', 'gene'):
            gene_cv = c; break
    if gene_cv: cv = cv[cv[gene_cv].isin(GENES)].copy()

    # Sig column
    sig_col = None
    for c in cv.columns:
        if 'clinicalsignificance' in c.lower().replace(' ', ''):
            sig_col = c; break
    if sig_col is None:
        print(f'    SKIP: no ClinicalSignificance'); continue

    # Dedup + classify
    cv[vid_col] = pd.to_numeric(cv[vid_col], errors='coerce')
    cv = cv.dropna(subset=[vid_col]).drop_duplicates(vid_col, keep='first')
    cv['_class'] = cv[sig_col].apply(classify_clinvar)

    # Merge
    merged = df.merge(
        cv[[vid_col, '_class']].rename(columns={vid_col: merge_key, '_class': f'class_{year}'}),
        on=merge_key, how='inner'
    )
    print(f'    Match: {len(merged):,} ({merge_key})  [{time.time()-t0:.0f}s]')

    # VUS baseline
    col_y = f'class_{year}'
    vus_y = merged[merged[col_y] == 'vus'].copy()
    n_vus = len(vus_y)
    if n_vus < 20:
        print(f'    Troppo poche VUS ({n_vus}), skip'); continue

    # Reclassificazioni
    vus_y['reclass_path'] = vus_y['class_current'].isin(['pathogenic', 'likely_pathogenic']).astype(int)
    vus_y['reclass_benign'] = vus_y['class_current'].isin(['benign', 'likely_benign']).astype(int)

    n_path = int(vus_y['reclass_path'].sum())
    n_benign = int(vus_y['reclass_benign'].sum())
    print(f'    VUS: {n_vus:,}  ->P/LP: {n_path}  ->B/LB: {n_benign}')

    if n_path < 3:
        print(f'    Troppo poche riclassif ({n_path}), skip dettagli')
        all_yearly.append({'year': year, 'n_vus': n_vus, 'n_path': n_path,
                            'OR_idk80': np.nan, 'p_idk80': np.nan,
                            'OR_hr': np.nan, 'OR_glob': np.nan})
        continue

    # Score identikit
    vus_y['has_idk'] = vus_y[gene_col].isin(idk.keys()).astype(int)
    vi = vus_y[vus_y['has_idk'] == 1].copy()

    scores = vi.apply(lambda r: score_vus(r, idk, gene_col, HR), axis=1)
    vi['idk_hits'] = [s[0] for s in scores]
    vi['idk_pct'] = [s[1] for s in scores]
    vi['idk_details'] = [s[2] for s in scores]
    vi['idk_80'] = (vi['idk_pct'] >= 0.8).astype(int)
    vi['idk_100'] = (vi['idk_pct'] >= 1.0).astype(int)

    # Soglia globale
    vi['glob_thr'] = ((vi[AM] >= 0.564) & (vi[RV] >= 0.75) &
                       (vi[AF].fillna(0) < 1e-4)).astype(int)

    n_vi = len(vi)
    n_rec_vi = int(vi['reclass_path'].sum())
    rec = vi[vi['reclass_path'] == 1]
    stab = vi[vi['reclass_path'] == 0]

    print(f'    Con identikit: {n_vi:,}  riclass: {n_rec_vi}')

    # --- TEST IDENTIKIT ---
    OR_idk80, p_idk80 = np.nan, np.nan
    OR_idk100, p_idk100 = np.nan, np.nan
    if n_rec_vi >= 2 and len(rec) > 0 and len(stab) > 0:
        for col_name, or_name, p_name in [('idk_80', 'OR_idk80', 'p_idk80'),
                                             ('idk_100', 'OR_idk100', 'p_idk100')]:
            tab = np.array([
                [int(rec[col_name].sum()), len(rec) - int(rec[col_name].sum())],
                [int(stab[col_name].sum()), len(stab) - int(stab[col_name].sum())]
            ])
            try:
                _or, _p = fisher_exact(tab)
            except:
                _or, _p = np.nan, np.nan
            if or_name == 'OR_idk80': OR_idk80, p_idk80 = _or, _p
            else: OR_idk100, p_idk100 = _or, _p

    # TEST HR
    OR_hr, p_hr = np.nan, np.nan
    if n_rec_vi >= 2:
        tab_hr = np.array([
            [int(rec['is_hr'].sum()), len(rec) - int(rec['is_hr'].sum())],
            [int(stab['is_hr'].sum()), len(stab) - int(stab['is_hr'].sum())]
        ])
        try: OR_hr, p_hr = fisher_exact(tab_hr)
        except: pass

    # TEST GLOBALE
    OR_glob, p_glob = np.nan, np.nan
    if n_rec_vi >= 2:
        tab_gl = np.array([
            [int(rec['glob_thr'].sum()), len(rec) - int(rec['glob_thr'].sum())],
            [int(stab['glob_thr'].sum()), len(stab) - int(stab['glob_thr'].sum())]
        ])
        try: OR_glob, p_glob = fisher_exact(tab_gl)
        except: pass

    # Sensitivity comparison
    sens_idk = rec['idk_80'].mean() if len(rec) > 0 else 0
    sens_hr = rec['is_hr'].mean() if len(rec) > 0 else 0
    sens_glob = rec['glob_thr'].mean() if len(rec) > 0 else 0
    sens_idk100 = rec['idk_100'].mean() if len(rec) > 0 else 0

    ppv_idk = int(rec['idk_80'].sum()) / int(vi['idk_80'].sum()) if int(vi['idk_80'].sum()) > 0 else 0
    ppv_hr = int(rec['is_hr'].sum()) / int(vi['is_hr'].sum()) if int(vi['is_hr'].sum()) > 0 else 0
    ppv_glob = int(rec['glob_thr'].sum()) / int(vi['glob_thr'].sum()) if int(vi['glob_thr'].sum()) > 0 else 0

    sig80 = '***' if p_idk80 < 0.001 else '**' if p_idk80 < 0.01 else '*' if p_idk80 < 0.05 else 'ns'
    print(f'    Identikit >=80%: OR={OR_idk80:.2f} p={p_idk80:.2e} {sig80}  '
          f'sens={100*sens_idk:.0f}%  ppv={100*ppv_idk:.1f}%')
    print(f'    Cluster HR:      OR={OR_hr:.2f} p={p_hr:.2e}      '
          f'sens={100*sens_hr:.0f}%  ppv={100*ppv_hr:.1f}%')
    print(f'    Soglia globale:  OR={OR_glob:.2f} p={p_glob:.2e}      '
          f'sens={100*sens_glob:.0f}%  ppv={100*ppv_glob:.1f}%')

    # Tasso per cluster
    print(f'    Tasso riclass per cluster:')
    for c in sorted(vi['cluster'].unique()):
        sub = vi[vi['cluster'] == c]
        r = sub['reclass_path'].mean()
        m = sub['idk_80'].mean()
        nr = int(sub['reclass_path'].sum())
        tag = ' <- HR' if c == HR else ''
        print(f'      C{c}: {nr:>3}/{len(sub):>5} = {100*r:.1f}%  '
              f'idk_match={100*m:.0f}%{tag}')

    # Store
    all_yearly.append({
        'year': year, 'n_vus': n_vus, 'n_vus_with_idk': n_vi,
        'n_path': n_path, 'n_path_idk': n_rec_vi,
        'OR_idk80': OR_idk80, 'p_idk80': p_idk80,
        'OR_idk100': OR_idk100, 'p_idk100': p_idk100,
        'OR_hr': OR_hr, 'p_hr': p_hr,
        'OR_glob': OR_glob, 'p_glob': p_glob,
        'sens_idk80': sens_idk, 'sens_hr': sens_hr, 'sens_glob': sens_glob,
        'ppv_idk80': ppv_idk, 'ppv_hr': ppv_hr, 'ppv_glob': ppv_glob,
    })

    all_method_comp.append({
        'year': year, 'n_rec': n_rec_vi,
        'sens_idk80': sens_idk, 'sens_idk100': sens_idk100,
        'sens_hr': sens_hr, 'sens_glob': sens_glob,
        'ppv_idk80': ppv_idk, 'ppv_hr': ppv_hr, 'ppv_glob': ppv_glob,
    })

    # Gene-level
    gres = []
    for g, sub in vi.groupby(gene_col):
        if len(sub) < 3: continue
        rg = sub[sub['reclass_path'] == 1]
        sg = sub[sub['reclass_path'] == 0]
        nr = int(sub['reclass_path'].sum())
        pm_r = rg['idk_80'].mean() if len(rg) > 0 else 0
        pm_s = sg['idk_80'].mean()
        OR_g, p_g = np.nan, np.nan
        if len(rg) >= 1 and len(sg) >= 1:
            t = np.array([[int(rg['idk_80'].sum()), len(rg) - int(rg['idk_80'].sum())],
                           [int(sg['idk_80'].sum()), len(sg) - int(sg['idk_80'].sum())]])
            try: OR_g, p_g = fisher_exact(t)
            except: pass
        gres.append({'year': year, 'gene': g, 'n_vus': len(sub),
                      'n_rec': nr, 'match_rec_80': 100 * pm_r,
                      'match_stab_80': 100 * pm_s, 'OR': OR_g, 'p': p_g})
    gdf = pd.DataFrame(gres).sort_values('n_rec', ascending=False)
    gdf.to_csv(f'{OUTDIR}/identikit_per_gene_{year}.csv', index=False)
    all_yearly_genes.append(gdf)

    # Pool riclassificate
    rec_export = rec.copy()
    rec_export['source_year'] = year
    all_reclass_pool.append(rec_export)

# =============================================================================
# 5. POOL ANALYSIS
# =============================================================================

print(f'\n{"="*70}')
print('[4/6] POOL ANALYSIS (tutte le riclassificate)')
print(f'{"="*70}')

if len(all_reclass_pool) > 0:
    pool = pd.concat(all_reclass_pool, ignore_index=True)
    pool = pool.drop_duplicates(subset='VariationID', keep='first')  # dedup cross-anno
    print(f'\n  Riclassificate uniche: {len(pool)}')

    # Score distribution
    print(f'\n  Score identikit nelle riclassificate:')
    for pct_thr in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        n = (pool['idk_pct'] >= pct_thr).sum()
        print(f'    >= {pct_thr:.0%}: {n:>4} ({100*n/len(pool):.1f}%)')

    # Top genes
    print(f'\n  Top geni riclassificati con identikit:')
    gene_pool = pool.groupby(gene_col).agg(
        n=('reclass_path', 'size'),
        idk_80_pct=('idk_80', 'mean'),
        idk_med=('idk_pct', 'median'),
    ).sort_values('n', ascending=False)
    for g, r in gene_pool.head(10).iterrows():
        print(f'    {g:10} n={int(r["n"]):>3}  '
              f'match_80={100*r["idk_80_pct"]:.0f}%  '
              f'score_med={r["idk_med"]:.2f}')

    pool.to_csv(f'{OUTDIR}/identikit_pool_analysis.csv', index=False)

# =============================================================================
# 6. RIEPILOGO + FIGURE
# =============================================================================

print(f'\n{"="*70}')
print('[5/6] RIEPILOGO')
print(f'{"="*70}')

summary = pd.DataFrame(all_yearly)
summary.to_csv(f'{OUTDIR}/identikit_temporal_summary.csv', index=False)

method_df = pd.DataFrame(all_method_comp)
method_df.to_csv(f'{OUTDIR}/identikit_method_comparison.csv', index=False)

print(f'\n  {"Anno":>4} {"VUS":>6} {"rec":>4} '
      f'{"OR_idk80":>9} {"p":>10} '
      f'{"OR_HR":>7} {"OR_glob":>8} '
      f'{"S_idk":>6} {"S_HR":>5} {"S_glob":>6}')
print(f'  {"-"*75}')

for _, r in summary.iterrows():
    if pd.isna(r.get('OR_idk80')): continue
    sig = '***' if r['p_idk80'] < 0.001 else '**' if r['p_idk80'] < 0.01 else '*' if r['p_idk80'] < 0.05 else ''
    or_hr = f"{r['OR_hr']:.1f}" if pd.notna(r['OR_hr']) else '-'
    or_gl = f"{r['OR_glob']:.1f}" if pd.notna(r['OR_glob']) else '-'
    print(f'  {int(r["year"]):>4} {int(r["n_vus"]):>6} {int(r["n_path"]):>4} '
          f'{r["OR_idk80"]:>8.2f}{sig:1} {r["p_idk80"]:>10.2e} '
          f'{or_hr:>7} {or_gl:>8} '
          f'{100*r.get("sens_idk80",0):>5.0f}% {100*r.get("sens_hr",0):>4.0f}% '
          f'{100*r.get("sens_glob",0):>5.0f}%')

# Conteggio anni significativi
n_sig = sum(1 for _, r in summary.iterrows()
            if pd.notna(r.get('p_idk80')) and r['p_idk80'] < 0.05)
n_tested = sum(1 for _, r in summary.iterrows() if pd.notna(r.get('OR_idk80')))
print(f'\n  Identikit significativo in {n_sig}/{n_tested} anni')

# Media sensibilita
if len(method_df) > 0:
    print(f'\n  Sensibilita media (tutti gli anni):')
    print(f'    Identikit >=80%: {100*method_df["sens_idk80"].mean():.1f}%')
    print(f'    Cluster HR:      {100*method_df["sens_hr"].mean():.1f}%')
    print(f'    Soglia globale:  {100*method_df["sens_glob"].mean():.1f}%')
    print(f'  PPV media:')
    print(f'    Identikit >=80%: {100*method_df["ppv_idk80"].mean():.1f}%')
    print(f'    Cluster HR:      {100*method_df["ppv_hr"].mean():.1f}%')
    print(f'    Soglia globale:  {100*method_df["ppv_glob"].mean():.1f}%')

# FIGURE
print(f'\n[6/6] Figure...')

if len(summary[summary['OR_idk80'].notna()]) >= 2:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    s = summary[summary['OR_idk80'].notna()].copy()

    # 6a. OR per anno
    ax = axes[0, 0]
    ax.plot(s['year'], s['OR_idk80'], 'o-', color='#d62728', label='Identikit >=80%', linewidth=2)
    if 'OR_hr' in s.columns:
        ax.plot(s['year'], s['OR_hr'], 's--', color='#1f77b4', label='Cluster HR')
    if 'OR_glob' in s.columns:
        ax.plot(s['year'], s['OR_glob'], '^--', color='#2ca02c', label='Soglia globale')
    ax.axhline(1, color='grey', linestyle=':', alpha=0.5)
    ax.set_xlabel('Anno ClinVar'); ax.set_ylabel('Odds Ratio')
    ax.set_title('OR per riclassificazione VUS->P/LP')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # 6b. Sensitivity per anno
    ax = axes[0, 1]
    if len(method_df) > 0:
        ax.plot(method_df['year'], 100 * method_df['sens_idk80'], 'o-',
                color='#d62728', label='Identikit >=80%', linewidth=2)
        ax.plot(method_df['year'], 100 * method_df['sens_hr'], 's--',
                color='#1f77b4', label='Cluster HR')
        ax.plot(method_df['year'], 100 * method_df['sens_glob'], '^--',
                color='#2ca02c', label='Soglia globale')
    ax.set_xlabel('Anno ClinVar'); ax.set_ylabel('Sensitivity (%)')
    ax.set_title('Sensitivity (% riclassificate catturate)')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # 6c. PPV per anno
    ax = axes[1, 0]
    if len(method_df) > 0:
        ax.plot(method_df['year'], 100 * method_df['ppv_idk80'], 'o-',
                color='#d62728', label='Identikit >=80%', linewidth=2)
        ax.plot(method_df['year'], 100 * method_df['ppv_hr'], 's--',
                color='#1f77b4', label='Cluster HR')
        ax.plot(method_df['year'], 100 * method_df['ppv_glob'], '^--',
                color='#2ca02c', label='Soglia globale')
    ax.set_xlabel('Anno ClinVar'); ax.set_ylabel('PPV (%)')
    ax.set_title('Positive Predictive Value')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # 6d. Gene heatmap (pool)
    ax = axes[1, 1]
    if len(all_yearly_genes) > 0:
        all_genes_df = pd.concat(all_yearly_genes)
        top_genes = all_genes_df.groupby('gene')['n_rec'].sum().nlargest(10).index
        pivot = all_genes_df[all_genes_df['gene'].isin(top_genes)].pivot_table(
            index='gene', columns='year', values='match_rec_80', aggfunc='first'
        )
        pivot = pivot.reindex(top_genes)
        sns.heatmap(pivot, cmap='RdYlGn', vmin=0, vmax=100, annot=True, fmt='.0f',
                    ax=ax, cbar_kws={'label': '% match identikit'}, linewidths=0.3)
        ax.set_title('Match identikit (%) nelle riclassificate')
    else:
        ax.text(0.5, 0.5, 'Dati insufficienti', ha='center', va='center')

    plt.suptitle('Validazione identikit gene-specifico: serie temporale ClinVar',
                 fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_identikit_validation.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  fig_identikit_validation.png')

# Lista file
print(f'\n  File in {OUTDIR}/:')
for f in sorted(os.listdir(OUTDIR)):
    print(f'    {f}')

print(f'\n{"="*70}')
print('VALIDAZIONE IDENTIKIT COMPLETATA')
print(f'{"="*70}')
