from pathlib import Path
import json
import tarfile

import pandas as pd
import pytest

from scripts.data import ROOT, audit_legacy, load_data, validate_counts, validate_metadata, verify_preservation


def test_source_integrity_and_experimental_units():
    meta, technical, htseq, biological, bio_meta = load_data()
    assert verify_preservation() == 31
    assert technical.shape[1] == 24
    assert list(biological.columns) == ['WT1','WT2','WT3','KO1','KO2','KO3']
    assert bio_meta.condition.value_counts().to_dict() == {'WildType':3,'Notch2CKO':3}
    assert technical.sum(axis=1).equals(biological.sum(axis=1))
    assert len(htseq) == 5
    assert set(htseq.index) == {'__no_feature','__ambiguous','__too_low_aQual',
                                '__not_aligned','__alignment_not_unique'}
    assert biological.ge(10).sum(axis=1).ge(3).sum() > 10000


def test_bad_counts_cannot_be_silently_coerced():
    base = pd.DataFrame({'WT1':[3,4], 'KO1':[0,5]},index=['Id4','Notch2'])
    for altered in [base.assign(WT1=[3,None]), base.assign(WT1=[3,-1]),
                    base.assign(WT1=[3,4.5]), base.set_axis(['Id4','Id4'])]:
        with pytest.raises(ValueError):
            validate_counts(altered)


def test_replicate_assignment_cannot_be_reused_or_swapped():
    meta=pd.read_csv(ROOT/'data/sample_sheet.csv')
    wrong=meta.copy(); wrong.loc[wrong.biological_sample.eq('KO1'),'condition']='WildType'
    with pytest.raises(ValueError, match='genotype'):
        validate_metadata(wrong)
    wrong=meta.copy(); wrong.loc[0,'lane']=2
    with pytest.raises(ValueError, match='libraries'):
        validate_metadata(wrong)


def test_original_table_audit():
    audit=audit_legacy()
    assert audit['rows']==29970
    assert audit['fdr_and_abs_lfc_at_least_1']==29
    assert audit['missing_padj'] > 0
