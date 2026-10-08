"""Run the complete, offline GSE116773 count-level reanalysis.

Usage: python -m scripts.analyze --out results
The study's 24 files are technical observations, collapsed to six biological units.
"""
from __future__ import annotations

import argparse
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import warnings

# Bound native parallelism before importing numerical libraries.
for variable in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(variable, '1')

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from joblib import parallel_config
from scipy.spatial.distance import pdist, squareform
from scipy.stats import false_discovery_control
from sklearn.decomposition import PCA

from scripts.data import ROOT, MARKERS, SAMPLES, audit_legacy, load_data, sha256, verify_preservation
from scripts.figures import make_figures


def write_table(frame: pd.DataFrame, path: Path, index: bool = True):
    frame.to_csv(path, index=index, float_format='%.12g', lineterminator='\n', na_rep='NA')


def classify(result: pd.DataFrame) -> pd.Series:
    status = pd.Series('FDR >= 0.05', index=result.index)
    status.loc[result.padj.isna()] = 'Not assigned an FDR'
    status.loc[result.padj.lt(.05) & result.log2FoldChange.gt(0)] = 'Higher in CKO'
    status.loc[result.padj.lt(.05) & result.log2FoldChange.lt(0)] = 'Lower in CKO'
    return status


def run(out: Path, threads: int = 2):
    if threads < 1:
        raise ValueError('--threads must be positive')
    out.mkdir(parents=True, exist_ok=True)
    tables = out / 'tables'; tables.mkdir(exist_ok=True)
    figures = out / 'figures'; figures.mkdir(exist_ok=True)
    # A stale success report must never survive a failed rerun.
    summary_path = out / 'summary.json'
    summary_path.unlink(missing_ok=True)
    (out / '.running').write_text('Analysis is incomplete until summary.json is written.\n')
    preserved = verify_preservation()
    meta, tech, htseq, bio, bio_meta = load_data()
    write_table(bio, tables / 'biological_counts.csv')
    write_table(bio_meta, tables / 'biological_samples.csv')
    write_table(htseq.T, tables / 'htseq_summary_counts.csv')
    library = meta.set_index('title')[['geo_accession', 'run_accession', 'biological_sample',
                                      'condition', 'library', 'lane']].copy()
    library['assigned_gene_counts'] = tech.sum()
    library['htseq_total'] = tech.sum() + htseq.sum()
    library['assigned_fraction'] = library.assigned_gene_counts / library.htseq_total
    write_table(library, tables / 'technical_library_qc.csv')
    tech_log = np.log2(tech.div(tech.sum(), axis=1) * 1e6 + 1)
    correlation = tech_log.corr(method='pearson')
    write_table(correlation, tables / 'technical_correlations.csv')

    # Condition-independent expression filter; three is the smaller biological group size.
    keep = bio.ge(10).sum(axis=1).ge(3)
    write_table(pd.DataFrame({'samples_with_count_ge_10': bio.ge(10).sum(axis=1),
                             'retained': keep}), tables / 'gene_filter.csv')
    counts = bio.loc[keep].T
    bio_meta['condition'] = pd.Categorical(bio_meta.condition, categories=['WildType', 'Notch2CKO'])
    inference = DefaultInference(n_cpus=threads)
    print(f'Fitting {counts.shape[1]:,} genes in {counts.shape[0]} biological samples.', flush=True)
    with warnings.catch_warnings(record=True) as caught, parallel_config(backend='threading'):
        warnings.simplefilter('always')
        dds = DeseqDataSet(counts=counts, metadata=bio_meta, design='~condition',
                           refit_cooks=True, min_replicates=7, inference=inference, quiet=False)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=['condition','Notch2CKO','WildType'],
                           alpha=.05, cooks_filter=True, independent_filter=True,
                           inference=inference, quiet=True)
        stats.summary()
        result = stats.results_df.copy()
        # MLE effect estimates and Wald SE are retained together: no posterior/Wald mixing.
        result['ci95_low'] = result.log2FoldChange - 1.959963984540054 * result.lfcSE
        result['ci95_high'] = result.log2FoldChange + 1.959963984540054 * result.lfcSE
        result['status'] = classify(result)
        result['fdr_and_abs_lfc_ge_1'] = result.padj.lt(.05) & result.log2FoldChange.abs().ge(1)
        # Sensitivity to independent filtering, using exactly the finite tested P-values.
        finite = result.pvalue.notna()
        result['padj_without_independent_filter'] = np.nan
        result.loc[finite, 'padj_without_independent_filter'] = false_discovery_control(
            result.loc[finite, 'pvalue'].to_numpy(), method='bh')
        normalized = pd.DataFrame(dds.layers['normed_counts'].T, index=counts.columns, columns=SAMPLES)
        model_diagnostics = dds.var.copy()
        size_factors = dds.obs['size_factors'].copy()
        dds.vst(use_design=False)
        vst = pd.DataFrame(dds.layers['vst_counts'].T, index=counts.columns, columns=SAMPLES)
    messages = sorted({f'{w.category.__name__}: {w.message}' for w in caught})
    for message in messages:
        print('Recorded warning:', message, flush=True)
    if not np.isfinite(vst.to_numpy()).all() or not np.isfinite(normalized.to_numpy()).all():
        raise ValueError('Nonfinite transformed values; cannot produce valid exploratory figures.')
    if not np.isfinite(result[['baseMean','log2FoldChange','lfcSE']].to_numpy()).all():
        raise ValueError('Nonfinite model estimates; inspect model diagnostics.')
    result.index.name = 'gene'
    ordered = result.sort_values(['padj','pvalue'], kind='stable', na_position='last')
    write_table(ordered, tables / 'differential_expression.csv')
    write_table(ordered[ordered.padj.lt(.05)], tables / 'significant_genes.csv')
    marker_results = result.reindex(MARKERS).copy()
    marker_results['retained_for_test'] = marker_results.index.isin(result.index)
    write_table(marker_results, tables / 'candidate_genes.csv')
    write_table(normalized, tables / 'normalized_counts.csv')
    write_table(vst, tables / 'vst_counts.csv')
    write_table(model_diagnostics, tables / 'model_diagnostics.csv')
    sample_qc = bio_meta.copy()
    sample_qc['total_assigned_counts'] = bio.sum()
    sample_qc['size_factor'] = size_factors
    sample_qc['detected_genes'] = bio.gt(0).sum()
    write_table(sample_qc, tables / 'biological_sample_qc.csv')

    # PCA genes are selected by variance without consulting genotype or DE P-values.
    pca_genes = vst.var(axis=1).sort_values(ascending=False, kind='stable').head(500).index
    pca = PCA(n_components=2, svd_solver='full')
    coords = pd.DataFrame(pca.fit_transform(vst.loc[pca_genes].T), index=SAMPLES, columns=['PC1','PC2'])
    write_table(coords, tables / 'pca_coordinates.csv')
    write_table(pd.DataFrame({'gene':pca_genes}), tables / 'pca_genes.csv', index=False)
    distances = pd.DataFrame(squareform(pdist(vst.T)), index=SAMPLES, columns=SAMPLES)
    write_table(distances, tables / 'sample_distances.csv')
    legacy = audit_legacy()
    (out / 'legacy_table_audit.json').write_text(json.dumps(legacy, indent=2)+'\n')
    convergence = {c: int((~model_diagnostics[c].astype(bool)).sum())
                   for c in ['_genewise_converged','_MAP_converged','_LFC_converged']
                   if c in model_diagnostics}
    summary = {
        'analysis': 'GSE116773 public HTSeq counts; new count-level reanalysis',
        'reference': 'GRCm38', 'design': '~ condition', 'contrast': 'Notch2CKO / WildType',
        'effect_estimator': 'Unshrunk maximum-likelihood log2 fold change',
        'source_genes': len(bio), 'technical_files': len(meta), 'biological_samples': 6,
        'biological_samples_per_condition': {'WildType':3,'Notch2CKO':3},
        'filter': 'At least 10 raw counts in at least 3 biological samples',
        'retained_genes': int(keep.sum()), 'removed_low_count_genes': int((~keep).sum()),
        'finite_pvalues': int(result.pvalue.notna().sum()), 'finite_padj': int(result.padj.notna().sum()),
        'missing_pvalues': int(result.pvalue.isna().sum()),
        'missing_padj': int(result.padj.isna().sum()),
        'fdr_below_0_05': int(result.padj.lt(.05).sum()),
        'higher_in_cko': int(result.status.eq('Higher in CKO').sum()),
        'lower_in_cko': int(result.status.eq('Lower in CKO').sum()),
        'fdr_and_abs_lfc_ge_1': int(result.fdr_and_abs_lfc_ge_1.sum()),
        'fdr_without_independent_filter': int(result.padj_without_independent_filter.lt(.05).sum()),
        'pca_variance_fraction': pca.explained_variance_ratio_.tolist(),
        'nonconverged_flags': convergence, 'warnings': messages,
        'original_files_preserved': preserved,
        'source_archive_sha256': sha256(ROOT / 'data/source/GSE116773_RAW.tar'),
        'legacy_audit': legacy,
    }
    make_figures(figures, tables, meta, tech, library, correlation, normalized, vst, result,
                 model_diagnostics, coords, pca.explained_variance_ratio_, distances)
    environment = {'python':platform.python_version(), 'system':platform.system(),
                   'packages':{p:version(p) for p in ['numpy','pandas','scipy','matplotlib','pydeseq2',
                                 'anndata','scikit-learn','formulaic','formulaic-contrasts','statsmodels']}}
    (out / 'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    (out / 'warnings.txt').write_text('\n'.join(messages) + ('\n' if messages else 'No Python warnings captured.\n'))
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    (out / '.running').unlink()
    print(json.dumps(summary,indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'results')
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    run(args.out.resolve(), args.threads)


if __name__ == '__main__':
    main()
