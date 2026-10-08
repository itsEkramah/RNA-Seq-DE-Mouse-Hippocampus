import json
from pathlib import Path

import pandas as pd
import fitz

from scripts.data import ROOT, MARKERS, SAMPLES


def test_full_run_outputs():
    path=ROOT/'results'
    assert (path/'summary.json').exists(), 'Run python -m scripts.analyze first.'
    summary=json.loads((path/'summary.json').read_text())
    assert not (path/'.running').exists()
    assert summary['biological_samples']==6
    assert summary['technical_files']==24
    assert summary['retained_genes']==len(pd.read_csv(path/'tables/differential_expression.csv'))
    result=pd.read_csv(path/'tables/differential_expression.csv')
    assert result.gene.is_unique
    assert result.padj.lt(.05).sum()==summary['fdr_below_0_05']
    assert result.pvalue.notna().sum()==summary['finite_pvalues']
    diagnostics=pd.read_csv(path/'tables/model_diagnostics.csv').set_index('gene')
    nonconverged=diagnostics.index[~diagnostics['_MAP_converged']]
    assert not set(nonconverged).intersection(result.loc[result.padj.lt(.05),'gene'])
    assert set(MARKERS)==set(pd.read_csv(path/'tables/candidate_genes.csv').gene)
    assert list(pd.read_csv(path/'tables/biological_counts.csv').columns[1:])==SAMPLES
    assert len(list((path/'figures').glob('*.png')))==7
    assert len(list((path/'figures').glob('*.pdf')))==7
    for pdf in (path/'figures').glob('*.pdf'):
        with fitz.open(pdf) as document:
            assert document.page_count == 1, pdf.name
            assert document[0].rect.width > 400 and document[0].rect.height > 300
