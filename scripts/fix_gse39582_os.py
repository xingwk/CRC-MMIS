#!/usr/bin/env python3
"""Fix GSE39582 OS endpoint only — minimal script"""
import pandas as pd, numpy as np
from lifelines import CoxPHFitter
import sys, os

BASE = r'C:\Users\xingw\WorkBuddy\2026-06-04-22-17-17'

# Load v1 results (already correct for GSE17536/14333/38832/33113)
v1 = pd.read_csv(os.path.join(BASE, 'results/revisions_r2/geo_validation_fixed.csv'))

# Load GSE39582 expression — use Python's faster approach
print('Loading expression...')
expr = pd.read_table(os.path.join(BASE, 'data/geo/GSE39582_series_matrix.txt.gz'), 
                      comment='!', index_col=0)
surv = pd.read_csv(os.path.join(BASE, 'data/geo/GSE39582_survival.csv'))
probes = pd.read_csv(os.path.join(BASE, 'data/geo/GPL570_probe_to_symbol.csv'))

# Map probes to genes using merge (faster and handles missing probes)
pm = probes.set_index('probe_id')['gene_symbol']
expr_genes = expr.copy()
expr_genes['gene'] = expr_genes.index.map(pm)
expr_genes = expr_genes.dropna(subset=['gene'])
expr_genes = expr_genes[expr_genes['gene'] != '']
expr_genes = expr_genes.groupby('gene').mean()
print(f'Gene matrix: {expr_genes.shape}')

# Load coefficients
for run_label, coef_file in [('Run1_seed42', 'results/v31_run01_seed42/selected_genes.csv'),
                               ('Run2_seed123', 'results/v31_run02_seed123/selected_genes.csv')]:
    coefs = pd.read_csv(os.path.join(BASE, coef_file), index_col=0)['coefficient']
    common = [g for g in coefs.index if g in expr_genes.index]
    print(f'\n{run_label}: {len(common)}/{len(coefs)} genes matched')
    
    # Risk score
    risk = np.zeros(expr_genes.shape[1])
    for g in common:
        risk += coefs[g] * expr_genes.loc[g].values
    risk = (risk - risk.mean()) / risk.std()
    
    # Match samples and use OS_time + OS_event
    surv_ids = surv['sample'].tolist()
    expr_ids = expr_genes.columns.tolist()
    common_ids = [s for s in surv_ids if s in expr_ids]
    surv_idx = [surv_ids.index(s) for s in common_ids]
    expr_idx = [expr_ids.index(s) for s in common_ids]
    
    t = surv.iloc[surv_idx]['os_time'].values.astype(float)  # OS time
    e = surv.iloc[surv_idx]['os_event'].values.astype(float)  # OS event
    r = risk[expr_idx]
    
    valid = ~np.isnan(t) & (t > 0) & ~np.isnan(e) & ~np.isnan(r)
    df = pd.DataFrame({'time': t[valid], 'event': e[valid], 'risk': r[valid]})
    
    cp = CoxPHFitter(penalizer=0.1)
    cp.fit(df, duration_col='time', event_col='event')
    hr = np.exp(cp.summary.loc['risk', 'coef'])
    p = cp.summary.loc['risk', 'p']
    ci = np.exp(cp.summary.loc['risk', ['coef lower 95%', 'coef upper 95%']].values)
    
    print(f'  N={valid.sum()}, events={int(e[valid].sum())}')
    print(f'  HR={hr:.3f} [{ci[0]:.3f}-{ci[1]:.3f}], p={p:.4f}')
    
    # Update v1
    mask = (v1['cohort'] == 'GSE39582') & (v1['run'] == run_label)
    v1.loc[mask, 'HR'] = hr
    v1.loc[mask, 'p'] = p
    v1.loc[mask, 'CI_low'] = ci[0]
    v1.loc[mask, 'CI_high'] = ci[1]
    v1.loc[mask, 'events'] = int(e[valid].sum())

# Save corrected
out = os.path.join(BASE, 'results/revisions_r2/geo_validation_corrected.csv')
v1.to_csv(out, index=False)
print(f'\n=== CORRECTED GSE39582 RESULTS ===')
print(v1[['cohort','run','n','HR','p','CI_low','CI_high','events']].to_string(index=False))
print(f'\nSaved: {out}')

# Quick meta
from lifelines.statistics import multivariate_logrank_test
# Simple RE meta (DerSimonian-Laird)
log_hr = np.log(v1['HR'].values)
se_hr = (np.log(v1['CI_high'].values) - np.log(v1['CI_low'].values)) / 3.92

# Weighted summary
w = 1.0 / se_hr**2
pooled_log_hr = np.sum(w * log_hr) / np.sum(w)
re_se = np.sqrt(1.0 / np.sum(w))
# I-squared
q = np.sum(w * (log_hr - pooled_log_hr)**2)
df = len(log_hr) - 1
c = sum(w) - sum(w**2)/sum(w)
i2 = max(0, (q - df) / q * 100) if q > df else 0

print(f'\n=== META-ANALYSIS (all valid cohorts) ===')
print(f'RE pooled HR: {np.exp(pooled_log_hr):.3f} [{np.exp(pooled_log_hr-1.96*re_se):.3f}-{np.exp(pooled_log_hr+1.96*re_se):.3f}]')
print(f'I2 = {i2:.1f}%')
print(f'Q = {q:.1f}, df = {df}')
