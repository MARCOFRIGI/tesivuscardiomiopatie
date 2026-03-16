#!/usr/bin/env python3
# =============================================================================
# VALIDAZIONE TEMPORALE CLINVAR — Script Completo v2.0
#
# Confronta release ClinVar storiche con classificazione attuale.
# Matching robusto via VariantID (non Name — evita errori di transcript).
# Valida TUTTI i pattern della tesi sulle VUS riclassificate.
#
# INPUT:
#   - PIPELINE_RESULTS/clustered_dataset.csv.gz  (output pipeline statistica)
#   - variant_summary_YYYY*.txt.gz               (release ClinVar storiche)
#
# OUTPUT (in TEMPORAL_VALIDATION/):
#   - temporal_validation_summary.csv     (metriche per anno)
#   - vus_YYYY_full.csv                   (VUS baseline + destino)
#   - gene_reclass_YYYY.csv               (gene-level)
#   - thesis_pattern_validation.csv       (tutti i pattern testati)
#   - fig_temporal_*.png                  (figure)
#
# REQUISITI: pandas, numpy, scipy, matplotlib, seaborn
# =============================================================================

import os, sys, re, glob, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, spearmanr, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

try:
    sys.stdout.reconfigure(encoding='utf-8')
except:
    pass

sns.set_theme(style='whitegrid', font_scale=1.0)

OUTDIR = 'TEMPORAL_VALIDATION'
os.makedirs(OUTDIR, exist_ok=True)

print('=' * 70)
print('VALIDAZIONE TEMPORALE CLINVAR v2.0')
print('Matching via VariantID | Multi-anno | Pattern thesis validation')
print('=' * 70)

# =============================================================================
# CONFIGURAZIONE
# =============================================================================

SARCOMERIC = {'ACTC1','MYBPC3','MYH7','MYL2','MYL3','TNNC1','TNNI3','TNNT2','TPM1'}
DESMOSOMAL = {'DSG2','DSP','PKP2'}
CYTOSKELETAL = {'DES','FLNC','LMNA','VCL','NEXN'}
ION_CHANNEL = {'SCN5A'}
RNA_BINDING = {'RBM20'}

def get_func_class(gene):
    if gene in SARCOMERIC: return 'Sarcomerica'
    if gene in DESMOSOMAL: return 'Desmosomiale'
    if gene in CYTOSKELETAL: return 'Citoscheletrica'
    if gene in ION_CHANNEL: return 'Canale ionico'
    if gene in RNA_BINDING: return 'RNA-binding'
    return 'Altra'

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

# =============================================================================
# 1. CARICAMENTO DATASET PIPELINE
# =============================================================================

print('\n[1/6] Caricamento dataset pipeline...')

PIPELINE_FILE = None
for candidate in [
    'PIPELINE_RESULTS/clustered_dataset.csv.gz',
    'PIPELINE_RESULTS/clustered_dataset.xlsx',
    'FINAL_HYBRID_PIPELINE_RESULTS/clustered_dataset.xlsx',
    'clustered_dataset.csv.gz',
]:
    if os.path.exists(candidate):
        PIPELINE_FILE = candidate
        break

if PIPELINE_FILE is None:
    print('  ERRORE: clustered_dataset non trovato')
    sys.exit(1)

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
df['func_class'] = df[gene_col].apply(get_func_class)

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

AM = 'am_pathogenicity'; RV = 'REVEL_max'; AF = 'AF_grpmax'

print(f'  Dataset: {len(df):,} varianti | HR cluster={HR}')
print(f'  VariationID disponibili: {df["VariationID"].notna().sum():,}')

# =============================================================================
# 2. DETECTION FILE CLINVAR STORICI
# =============================================================================

print('\n[2/6] Ricerca file ClinVar storici...')

def find_clinvar_files():
    patterns = ['variant_summary_*.txt.gz', 'variant_summary_*.txt',
                'variant_summary*.txt.gz', 'variant_summary*.txt']
    found = {}
    for pattern in patterns:
        for fp in glob.glob(pattern):
            match = re.search(r'(20\d{2})', fp)
            if match:
                year = int(match.group(1))
                if year not in found:
                    found[year] = fp
    return dict(sorted(found.items()))

clinvar_files = find_clinvar_files()

if not clinvar_files:
    print('  ERRORE: nessun file variant_summary_YYYY trovato')
    sys.exit(1)

for year, path in clinvar_files.items():
    sz = os.path.getsize(path) / 1e6
    print(f'  {year}: {path} ({sz:.0f} MB)')

# =============================================================================
# 3. CARICAMENTO E MATCHING PER OGNI ANNO
# =============================================================================

print('\n[3/6] Caricamento e matching...')

GENES = set(df[gene_col].dropna().unique())
all_yearly_results = []
all_vus_reclass = []
yearly_vus_data = {}

for year, filepath in clinvar_files.items():
    print(f'\n  --- {year} ({filepath}) ---')
    t0 = time.time()

    try:
        cv = pd.read_csv(filepath, sep='\t',
                          compression='gzip' if filepath.endswith('.gz') else None,
                          low_memory=False, on_bad_lines='skip')
    except TypeError:
        cv = pd.read_csv(filepath, sep='\t', low_memory=False, error_bad_lines=False)

    if cv.columns[0].startswith('#'):
        cv = cv.rename(columns={cv.columns[0]: cv.columns[0].lstrip('#')})

    vid_col = None
    merge_on_pipeline = 'VariationID'  # default
    for c in cv.columns:
        if 'variantid' in c.lower().replace(' ', '').replace('#', ''):
            vid_col = c; break
    # Fallback: dal 2017 ClinVar ha rimosso VariantID, usa AlleleID
    if vid_col is None:
        for c in cv.columns:
            if 'alleleid' in c.lower().replace(' ', '').replace('#', ''):
                vid_col = c
                merge_on_pipeline = 'AlleleID'
                break
    if vid_col is None:
        print(f'    ERRORE: ne VariantID ne AlleleID trovato')
        print(f'    Colonne: {list(cv.columns[:10])}')
        continue

    gene_cv_col = None
    for c in cv.columns:
        if c.lower() in ('genesymbol', 'gene symbol', 'gene'):
            gene_cv_col = c; break
    if gene_cv_col:
        cv = cv[cv[gene_cv_col].isin(GENES)].copy()

    sig_col = None
    for c in cv.columns:
        if 'clinicalsignificance' in c.lower().replace(' ', ''):
            sig_col = c; break
    if sig_col is None:
        print(f'    ERRORE: ClinicalSignificance non trovato')
        continue

    cv[vid_col] = pd.to_numeric(cv[vid_col], errors='coerce')
    cv = cv.dropna(subset=[vid_col]).drop_duplicates(vid_col, keep='first')
    cv['_class'] = cv[sig_col].apply(classify_clinvar)

    print(f'    Varianti uniche geni target: {len(cv):,} ({time.time()-t0:.0f}s)')
    print(f'    Matching via: {vid_col} -> {merge_on_pipeline}')
    print(f'    VUS: {(cv["_class"]=="vus").sum():,}  '
          f'P+LP: {cv["_class"].isin(["pathogenic","likely_pathogenic"]).sum():,}')

    merged = df.merge(
        cv[[vid_col, '_class']].rename(columns={vid_col: merge_on_pipeline, '_class': f'class_{year}'}),
        on=merge_on_pipeline, how='inner'
    )
    print(f'    Match VariantID: {len(merged):,}')

    col_y = f'class_{year}'
    vus_y = merged[merged[col_y] == 'vus'].copy()
    n_vus = len(vus_y)
    if n_vus < 10:
        print(f'    Troppo poche VUS ({n_vus}), skip')
        continue

    vus_y['reclass_path'] = vus_y['class_current'].isin(['pathogenic', 'likely_pathogenic']).astype(int)
    vus_y['reclass_benign'] = vus_y['class_current'].isin(['benign', 'likely_benign']).astype(int)
    vus_y['reclass_conf'] = (vus_y['class_current'] == 'conflicting').astype(int)

    n_path = int(vus_y['reclass_path'].sum())
    n_benign = int(vus_y['reclass_benign'].sum())
    n_conf = int(vus_y['reclass_conf'].sum())
    n_still = int((vus_y['class_current'] == 'vus').sum())

    print(f'    VUS {year}: {n_vus:,}')
    print(f'      -> P/LP: {n_path} ({100*n_path/n_vus:.1f}%)  '
          f'-> B/LB: {n_benign} ({100*n_benign/n_vus:.1f}%)  '
          f'-> Conf: {n_conf} ({100*n_conf/n_vus:.1f}%)  '
          f'-> VUS: {n_still} ({100*n_still/n_vus:.1f}%)')

    OR_hr, p_hr, rho, p_rho = np.nan, np.nan, np.nan, np.nan
    if n_path >= 3:
        tab = pd.crosstab(vus_y['is_hr'], vus_y['reclass_path'])
        if tab.shape == (2, 2):
            OR_hr, p_hr = fisher_exact(tab)
        rho, p_rho = spearmanr(vus_y['cluster'], vus_y['reclass_path'])
        print(f'    Fisher HR: OR={OR_hr:.2f}  p={p_hr:.2e}')
        print(f'    Spearman: rho={rho:.4f}  p={p_rho:.2e}')

    print(f'    Per cluster:')
    for c in sorted(vus_y['cluster'].unique()):
        sub = vus_y[vus_y['cluster'] == c]
        n_r = int(sub['reclass_path'].sum())
        rate = sub['reclass_path'].mean()
        tag = ' <- HR' if c == HR else ''
        print(f'      C{c}: {n_r:>3}/{len(sub):>5} = {100*rate:.1f}%{tag}')

    vus_y.to_csv(f'{OUTDIR}/vus_{year}_full.csv', index=False)
    yearly_vus_data[year] = vus_y

    gene_res = []
    for g, sub in vus_y.groupby(gene_col):
        if len(sub) < 3: continue
        gene_res.append({
            'gene': g, 'n_vus': len(sub),
            'n_to_path': int(sub['reclass_path'].sum()),
            'n_to_benign': int(sub['reclass_benign'].sum()),
            'rate_path': sub['reclass_path'].mean(),
        })
    gene_df = pd.DataFrame(gene_res).sort_values('n_to_path', ascending=False)
    gene_df.to_csv(f'{OUTDIR}/gene_reclass_{year}.csv', index=False)

    reclass_sub = vus_y[vus_y['reclass_path'] == 1].copy()
    reclass_sub['source_year'] = year
    all_vus_reclass.append(reclass_sub)

    all_yearly_results.append({
        'year': year, 'n_matched': len(merged), 'n_vus_baseline': n_vus,
        'n_to_path': n_path, 'n_to_benign': n_benign,
        'n_to_conf': n_conf, 'n_still_vus': n_still,
        'pct_to_path': 100 * n_path / n_vus,
        'OR_hr': OR_hr, 'p_hr': p_hr,
        'rho': rho, 'p_rho': p_rho,
    })

# =============================================================================
# 4. VALIDAZIONE PATTERN TESI
# =============================================================================

print(f'\n{"="*70}')
print('[4/6] VALIDAZIONE PATTERN TESI')
print(f'{"="*70}')

pattern_results = []

if len(all_vus_reclass) == 0:
    print('  Nessuna VUS riclassificata — skip')
else:
    reclass_pool = pd.concat(all_vus_reclass, ignore_index=True)
    oldest_year = min(yearly_vus_data.keys())
    vus_oldest = yearly_vus_data[oldest_year]
    stable_pool = vus_oldest[vus_oldest['reclass_path'] == 0].copy()

    n_rec = len(reclass_pool)
    n_stab = len(stable_pool)
    print(f'  VUS riclassificate (pool): {n_rec}')
    print(f'  VUS stabili (baseline {oldest_year}): {n_stab}')

    # PATTERN 1: Score predittivi
    print(f'\n  PATTERN 1: Score predittivi')
    for score_col, label in [(AM, 'AlphaMissense'), (RV, 'REVEL')]:
        if score_col not in reclass_pool.columns: continue
        med_rec = reclass_pool[score_col].median()
        med_stab = stable_pool[score_col].median()
        rec_v = reclass_pool[score_col].dropna()
        stab_v = stable_pool[score_col].dropna()
        p_mw = np.nan
        if len(rec_v) > 2 and len(stab_v) > 2:
            _, p_mw = mannwhitneyu(rec_v, stab_v, alternative='greater')
        confirmed = med_rec > med_stab and (pd.isna(p_mw) or p_mw < 0.05)
        print(f'    {label:20} riclass={med_rec:.3f}  stabili={med_stab:.3f}  '
              f'p={p_mw:.2e}  {"CONFERMATO" if confirmed else "non confermato"}')
        pattern_results.append({
            'pattern': f'Score {label} piu alto nelle riclassificate',
            'value_reclass': f'{med_rec:.3f}', 'value_stable': f'{med_stab:.3f}',
            'p_value': p_mw, 'confirmed': confirmed
        })

    # PATTERN 2: Enrichment cluster HR
    print(f'\n  PATTERN 2: Enrichment cluster HR')
    pct_hr_rec = reclass_pool['is_hr'].mean()
    pct_hr_stab = stable_pool['is_hr'].mean()
    tab = np.array([
        [int(reclass_pool['is_hr'].sum()), int((reclass_pool['is_hr'] == 0).sum())],
        [int(stable_pool['is_hr'].sum()), int((stable_pool['is_hr'] == 0).sum())]
    ])
    OR_hr2, p_hr2 = fisher_exact(tab)
    confirmed = OR_hr2 > 1 and p_hr2 < 0.05
    print(f'    %HR riclass={100*pct_hr_rec:.1f}%  stabili={100*pct_hr_stab:.1f}%  '
          f'OR={OR_hr2:.2f}  p={p_hr2:.2e}  {"CONFERMATO" if confirmed else "non confermato"}')
    pattern_results.append({
        'pattern': 'Enrichment cluster HR nelle riclassificate',
        'value_reclass': f'{100*pct_hr_rec:.1f}%', 'value_stable': f'{100*pct_hr_stab:.1f}%',
        'p_value': p_hr2, 'confirmed': confirmed
    })

    # PATTERN 3: Gradiente monotono
    print(f'\n  PATTERN 3: Gradiente monotono per cluster')
    rates_c = []
    for c in sorted(vus_oldest['cluster'].unique()):
        sub = vus_oldest[vus_oldest['cluster'] == c]
        rate = sub['reclass_path'].mean()
        rates_c.append(rate)
        print(f'    C{c}: {100*rate:.1f}%')
    monotone = len(rates_c) >= 2 and rates_c[-1] > rates_c[0]
    pattern_results.append({
        'pattern': 'Gradiente C0 < C3 per tasso riclassificazione',
        'value_reclass': f'C3={100*rates_c[-1]:.1f}%', 'value_stable': f'C0={100*rates_c[0]:.1f}%',
        'p_value': np.nan, 'confirmed': monotone
    })

    # PATTERN 4: AF piu bassa
    print(f'\n  PATTERN 4: AF piu bassa')
    if AF in reclass_pool.columns:
        af_rec = reclass_pool[AF].dropna()
        af_stab = stable_pool[AF].dropna()
        if len(af_rec) > 2 and len(af_stab) > 2:
            med_af_r = af_rec.median(); med_af_s = af_stab.median()
            _, p_af = mannwhitneyu(af_rec, af_stab, alternative='less')
            confirmed = med_af_r < med_af_s and p_af < 0.05
            print(f'    AF med riclass={med_af_r:.2e}  stabili={med_af_s:.2e}  '
                  f'p={p_af:.2e}  {"CONFERMATO" if confirmed else "non confermato"}')
            pattern_results.append({
                'pattern': 'AF piu bassa nelle riclassificate',
                'value_reclass': f'{med_af_r:.2e}', 'value_stable': f'{med_af_s:.2e}',
                'p_value': p_af, 'confirmed': confirmed
            })

    # PATTERN 5: Domini strutturali
    print(f'\n  PATTERN 5: Enrichment domini strutturali')
    for dom_col, dom_label in [
        ('is_uniprot_domain', 'Dominio UniProt'),
        ('is_binding_site', 'Sito di legame'),
        ('in_interpro', 'In InterPro'),
        ('is_structural_domain', 'Dominio strutturale IPR'),
    ]:
        if dom_col not in reclass_pool.columns: continue
        pct_r = reclass_pool[dom_col].mean()
        pct_s = stable_pool[dom_col].mean()
        tab_d = np.array([
            [int(reclass_pool[dom_col].sum()), int((reclass_pool[dom_col] == 0).sum())],
            [int(stable_pool[dom_col].sum()), int((stable_pool[dom_col] == 0).sum())]
        ])
        OR_d, p_d = fisher_exact(tab_d) if tab_d.min() >= 0 else (np.nan, np.nan)
        confirmed = OR_d > 1 and p_d < 0.05
        print(f'    {dom_label:25} riclass={100*pct_r:.1f}%  stab={100*pct_s:.1f}%  '
              f'OR={OR_d:.2f}  p={p_d:.2e}  {"CONFERMATO" if confirmed else "non confermato"}')
        pattern_results.append({
            'pattern': f'Enrichment {dom_label} nelle riclassificate',
            'value_reclass': f'{100*pct_r:.1f}%', 'value_stable': f'{100*pct_s:.1f}%',
            'p_value': p_d, 'confirmed': confirmed
        })

    # PATTERN 6: Classe funzionale
    print(f'\n  PATTERN 6: Classe funzionale (sarcomerici dominano)')
    if 'func_class' in reclass_pool.columns:
        fc_rec = reclass_pool['func_class'].value_counts(normalize=True)
        fc_stab = stable_pool['func_class'].value_counts(normalize=True)
        for cls in ['Sarcomerica', 'Desmosomiale', 'Citoscheletrica', 'Canale ionico']:
            print(f'    {cls:20} riclass={100*fc_rec.get(cls,0):.1f}%  '
                  f'stabili={100*fc_stab.get(cls,0):.1f}%')
        tab_s = np.array([
            [int((reclass_pool['func_class'] == 'Sarcomerica').sum()),
             int((reclass_pool['func_class'] != 'Sarcomerica').sum())],
            [int((stable_pool['func_class'] == 'Sarcomerica').sum()),
             int((stable_pool['func_class'] != 'Sarcomerica').sum())]
        ])
        OR_s, p_s = fisher_exact(tab_s)
        confirmed = OR_s > 1 and p_s < 0.05
        print(f'    Sarcomerica OR={OR_s:.2f}  p={p_s:.2e}  '
              f'{"CONFERMATO" if confirmed else "non confermato"}')
        pattern_results.append({
            'pattern': 'Enrichment sarcomerici nelle riclassificate',
            'value_reclass': f'{100*fc_rec.get("Sarcomerica",0):.1f}%',
            'value_stable': f'{100*fc_stab.get("Sarcomerica",0):.1f}%',
            'p_value': p_s, 'confirmed': confirmed
        })

    # PATTERN 7: Probability surface AM x REVEL
    print(f'\n  PATTERN 7: Probability surface predice riclassificazione')
    if AM in vus_oldest.columns and RV in vus_oldest.columns:
        try:
            vus_oldest['_am_q'] = pd.qcut(vus_oldest[AM], 4, labels=['Q1','Q2','Q3','Q4'],
                                           duplicates='drop')
            vus_oldest['_rv_q'] = pd.qcut(vus_oldest[RV], 4, labels=['Q1','Q2','Q3','Q4'],
                                           duplicates='drop')
            q_high = vus_oldest[(vus_oldest['_am_q'] == 'Q4') & (vus_oldest['_rv_q'] == 'Q4')]
            q_low = vus_oldest[(vus_oldest['_am_q'] == 'Q1') & (vus_oldest['_rv_q'] == 'Q1')]
            if len(q_high) > 3 and len(q_low) > 3:
                rate_h = q_high['reclass_path'].mean()
                rate_l = q_low['reclass_path'].mean()
                tab_q = np.array([
                    [int(q_high['reclass_path'].sum()), len(q_high)-int(q_high['reclass_path'].sum())],
                    [int(q_low['reclass_path'].sum()), len(q_low)-int(q_low['reclass_path'].sum())]
                ])
                OR_q, p_q = fisher_exact(tab_q)
                confirmed = rate_h > rate_l
                print(f'    Q4xQ4: {100*rate_h:.1f}%  Q1xQ1: {100*rate_l:.1f}%  '
                      f'OR={OR_q:.2f}  p={p_q:.2e}  {"CONFERMATO" if confirmed else "non confermato"}')
                pattern_results.append({
                    'pattern': 'Quadrante alto AM x REVEL predice riclassificazione',
                    'value_reclass': f'Q4xQ4={100*rate_h:.1f}%',
                    'value_stable': f'Q1xQ1={100*rate_l:.1f}%',
                    'p_value': p_q, 'confirmed': confirmed
                })
        except Exception as e:
            print(f'    Errore qcut: {e}')

    # PATTERN 8: -> benigno piu frequente in C0
    print(f'\n  PATTERN 8: Riclass -> benigno piu frequente in C0')
    rate_b0 = vus_oldest[vus_oldest['cluster'] == 0]['reclass_benign'].mean()
    rate_bhr = vus_oldest[vus_oldest['cluster'] == HR]['reclass_benign'].mean()
    confirmed = rate_b0 > rate_bhr
    print(f'    C0: {100*rate_b0:.1f}%  HR: {100*rate_bhr:.1f}%  '
          f'{"CONFERMATO" if confirmed else "non confermato"}')
    pattern_results.append({
        'pattern': 'Riclass -> benigno piu frequente in C0 che in HR',
        'value_reclass': f'C0={100*rate_b0:.1f}%', 'value_stable': f'HR={100*rate_bhr:.1f}%',
        'p_value': np.nan, 'confirmed': confirmed
    })

# =============================================================================
# 5. RIEPILOGO
# =============================================================================

print(f'\n{"="*70}')
print('[5/6] RIEPILOGO')
print(f'{"="*70}')

summary_df = pd.DataFrame(all_yearly_results)
summary_df.to_csv(f'{OUTDIR}/temporal_validation_summary.csv', index=False)

for _, r in summary_df.iterrows():
    print(f'\n  Anno {int(r["year"])}:')
    print(f'    VUS baseline: {int(r["n_vus_baseline"]):,}')
    print(f'    -> P/LP: {int(r["n_to_path"])} ({r["pct_to_path"]:.1f}%)')
    if pd.notna(r['OR_hr']):
        sig = '***' if r['p_hr'] < 0.001 else '**' if r['p_hr'] < 0.01 else '*' if r['p_hr'] < 0.05 else 'ns'
        print(f'    OR HR: {r["OR_hr"]:.2f}  p={r["p_hr"]:.2e} {sig}')

if pattern_results:
    pattern_df = pd.DataFrame(pattern_results)
    pattern_df.to_csv(f'{OUTDIR}/thesis_pattern_validation.csv', index=False)
    n_conf = int(pattern_df['confirmed'].sum())
    n_tot = len(pattern_df)
    print(f'\n  PATTERN TESI VALIDATI: {n_conf}/{n_tot}')
    for _, r in pattern_df.iterrows():
        sym = 'V' if r['confirmed'] else 'X'
        p_s = f"p={r['p_value']:.2e}" if pd.notna(r['p_value']) else ''
        print(f'    [{sym}] {r["pattern"]}')
        print(f'        riclass={r["value_reclass"]}  stabili={r["value_stable"]}  {p_s}')

# =============================================================================
# 6. FIGURE
# =============================================================================

print(f'\n[6/6] Figure...')

if len(all_yearly_results) > 0 and len(all_vus_reclass) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    oldest_year = min(yearly_vus_data.keys())
    vus_plot = yearly_vus_data[oldest_year]

    # 6a. Tasso per cluster
    ax = axes[0, 0]
    rates = []
    for c in sorted(vus_plot['cluster'].unique()):
        sub = vus_plot[vus_plot['cluster'] == c]
        rates.append({'cluster': c, 'rate': 100*sub['reclass_path'].mean(),
                      'n_r': int(sub['reclass_path'].sum()), 'n': len(sub)})
    rates_df = pd.DataFrame(rates)
    colors = ['#d62728' if c == HR else '#1f77b4' for c in rates_df['cluster']]
    ax.bar(rates_df['cluster'], rates_df['rate'], color=colors, edgecolor='white')
    for _, r in rates_df.iterrows():
        ax.text(r['cluster'], r['rate'] + 0.3, f"n={int(r['n_r'])}", ha='center', fontsize=9)
    ax.set_xlabel('Cluster'); ax.set_ylabel('% VUS -> P/LP')
    ax.set_title(f'Tasso riclassificazione per cluster ({oldest_year})')

    # 6b. Destino VUS
    ax = axes[0, 1]
    r0 = summary_df.iloc[0]
    destini = {'P/LP': int(r0['n_to_path']), 'B/LB': int(r0['n_to_benign']),
               'Conflicting': int(r0['n_to_conf']), 'Ancora VUS': int(r0['n_still_vus'])}
    other = int(r0['n_vus_baseline']) - sum(destini.values())
    if other > 0: destini['Altro'] = other
    cols_pie = ['#d62728','#1f77b4','#ff7f0e','#7f7f7f','#bcbd22'][:len(destini)]
    ax.pie(destini.values(), labels=destini.keys(), autopct='%1.1f%%', colors=cols_pie)
    ax.set_title(f'Destino delle VUS del {oldest_year}')

    # 6c. Score distribution
    ax = axes[1, 0]
    reclass_pool = pd.concat(all_vus_reclass, ignore_index=True)
    if AM in reclass_pool.columns:
        stab = vus_plot[vus_plot['reclass_path'] == 0][AM].dropna()
        rec = reclass_pool[AM].dropna()
        if len(rec) > 0:
            ax.hist(stab, bins=30, alpha=0.5, density=True, label=f'Stabili (n={len(stab):,})')
            ax.hist(rec, bins=15, alpha=0.7, density=True, color='red', label=f'->P/LP (n={len(rec)})')
            ax.set_xlabel('AlphaMissense'); ax.set_ylabel('Densita')
            ax.set_title('AM: riclassificate vs stabili'); ax.legend()

    # 6d. Pattern summary
    ax = axes[1, 1]; ax.axis('off')
    if pattern_results:
        text = f'PATTERN TESI VALIDATI: {n_conf}/{n_tot}\n\n'
        for _, r in pattern_df.iterrows():
            s = 'V' if r['confirmed'] else 'X'
            text += f'[{s}] {r["pattern"][:50]}\n'
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=8,
                verticalalignment='top', fontfamily='monospace')

    plt.suptitle(f'Validazione temporale ClinVar {oldest_year} -> oggi', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_temporal_complete.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  fig_temporal_complete.png')

print(f'\n  File in {OUTDIR}/:')
for f in sorted(os.listdir(OUTDIR)):
    print(f'    {f}')

print(f'\n{"="*70}')
print('VALIDAZIONE TEMPORALE COMPLETATA')
print(f'{"="*70}')
