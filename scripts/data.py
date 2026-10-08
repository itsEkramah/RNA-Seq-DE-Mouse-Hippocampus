"""Validate published counts, resolve experimental units, and preserve provenance."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ['WT1', 'WT2', 'WT3', 'KO1', 'KO2', 'KO3']
MARKERS = ['Notch2', 'Id4', 'Hes5', 'Hopx', 'Ascl1', 'Egfr', 'Eomes', 'Mki67']


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_counts(counts: pd.DataFrame) -> None:
    if counts.empty or counts.index.has_duplicates or counts.columns.has_duplicates:
        raise ValueError('Counts must be nonempty with unique genes and sample names.')
    if counts.index.isna().any() or counts.columns.isna().any():
        raise ValueError('Missing gene or sample identifier.')
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Counts must be finite and nonnegative; missing values are not zero.')
    if (values != np.floor(values)).any():
        raise ValueError('Raw integer counts are required; do not round normalized values.')
    if (values > np.iinfo(np.int64).max / 24).any():
        raise ValueError('Counts exceed the safe aggregation range.')
    if (counts.sum(axis=0) == 0).any():
        raise ValueError('A sample contains no assigned gene counts.')


def validate_metadata(meta: pd.DataFrame) -> None:
    required = {'geo_accession', 'run_accession', 'title', 'biological_sample',
                'condition', 'library', 'lane', 'filename', 'sha256'}
    if not required.issubset(meta.columns) or meta[list(required)].isna().any().any():
        raise ValueError('Sample metadata are missing required fields.')
    if len(meta) != 24 or any(meta[c].duplicated().any() for c in
                             ['geo_accession', 'run_accession', 'title', 'filename']):
        raise ValueError('Expected 24 distinct technical files and accessions.')
    if set(meta.biological_sample) != set(SAMPLES):
        raise ValueError('Expected WT1-WT3 and KO1-KO3 biological units.')
    for sample, block in meta.groupby('biological_sample'):
        expected = 'WildType' if sample.startswith('WT') else 'Notch2CKO'
        if set(block.condition) != {expected}:
            raise ValueError('Inconsistent genotype within a biological sample.')
        if set(zip(block.library, block.lane)) != {(1,1), (1,2), (2,1), (2,2)} or len(block) != 4:
            raise ValueError('Each biological sample must have two libraries and two lanes.')
        expected_titles = {f'{sample}_rep{lib}_{lane}' for lib in (1,2) for lane in (1,2)}
        if set(block.title) != expected_titles:
            raise ValueError('GEO titles do not agree with library and biological sample mapping.')


def validate_against_geo(meta: pd.DataFrame, excerpt: Path, ena_file: Path) -> None:
    """Cross-check the curated sheet against independent GEO and ENA fields."""
    text = excerpt.read_text(encoding='utf8')
    blocks = re.split(r'(?=^\^SAMPLE = )', text, flags=re.M)[1:]
    if len(blocks) != len(meta):
        raise ValueError('GEO excerpt has an unexpected number of sample records.')
    source = {}
    for block in blocks:
        accession = re.search(r'^\^SAMPLE = (\S+)',block,re.M).group(1)
        source[accession] = {
            'title': re.search(r'^!Sample_title = (\S+)',block,re.M).group(1),
            'genotype': re.search(r'^!Sample_characteristics_ch1 = genotype: (.+)$',block,re.M).group(1).strip(),
            'experiment': re.search(r'SRA: .*term=(SRX\d+)',block).group(1),
            'filename': re.search(r'^!Sample_supplementary_file_1 = .*/([^/\r\n]+)$',block,re.M).group(1),
        }
    ena = pd.read_csv(ena_file,sep='\t').set_index('experiment_accession')
    if set(meta.experiment_accession) != set(ena.index) or len(ena) != 24:
        raise ValueError('ENA and curated experiment accessions differ.')
    for row in meta.itertuples():
        record=source[row.geo_accession]
        genotype='Wild type' if row.condition == 'WildType' else 'Notch2 CKO'
        if (record['title'] != row.title or record['genotype'] != genotype or
            record['experiment'] != row.experiment_accession or record['filename'] != row.filename or
            ena.loc[row.experiment_accession,'run_accession'] != row.run_accession or
            ena.loc[row.experiment_accession,'library_layout'] != 'SINGLE'):
            raise ValueError(f'GEO/ENA mismatch for {row.geo_accession}.')


def load_data(root: Path = ROOT):
    source = root / 'data/source'
    provenance = json.loads((source / 'provenance.json').read_text())
    for name, digest in provenance['files'].items():
        if sha256(source / name) != digest:
            raise ValueError(f'Source checksum mismatch: {name}')
    if sha256(root / 'data/sample_sheet.csv') != provenance['sample_sheet_sha256']:
        raise ValueError('Sample sheet differs from the reviewed source mapping.')
    meta = pd.read_csv(root / 'data/sample_sheet.csv')
    validate_metadata(meta)
    validate_against_geo(meta, source/'geo_metadata_excerpt.soft', source/'ena_runs.tsv')
    series = []
    with tarfile.open(source / 'GSE116773_RAW.tar') as archive:
        members = archive.getmembers()
        if len(members) != len(meta) or {m.name for m in members} != set(meta.filename):
            raise ValueError('Archive members differ from the sample sheet.')
        for row in meta.itertuples():
            raw = archive.extractfile(row.filename).read()
            if hashlib.sha256(raw).hexdigest() != row.sha256:
                raise ValueError(f'Count-file checksum mismatch: {row.filename}')
            frame = pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep='\t', header=None,
                                names=['gene', row.title], index_col=0)
            validate_counts(frame)
            series.append(frame.iloc[:, 0])
    technical = pd.concat(series, axis=1)
    validate_counts(technical)  # rejects mismatched gene sets rather than filling with zero
    technical = technical.astype('int64')
    special = technical.loc[technical.index.str.startswith('__')]
    genes = technical.loc[~technical.index.str.startswith('__')]
    biological = pd.DataFrame({s: genes[meta.loc[meta.biological_sample.eq(s), 'title']].sum(axis=1)
                               for s in SAMPLES})
    validate_counts(biological)
    if not np.array_equal(biological.sum(axis=1), genes.sum(axis=1)):
        raise ValueError('Technical aggregation failed to conserve counts.')
    bio_meta = pd.DataFrame({'condition': ['WildType']*3 + ['Notch2CKO']*3}, index=SAMPLES)
    bio_meta.index.name = 'sample'
    return meta, genes, special, biological, bio_meta


def audit_legacy(root: Path = ROOT) -> dict:
    table = pd.read_csv(root / '05_DE_analysis/deseq2_results.csv')
    if table.gene.duplicated().any() or not table.gene.equals(table.iloc[:, 0]):
        raise ValueError('Original table has inconsistent or duplicated gene identifiers.')
    for col in ['pvalue', 'padj']:
        if not table[col].dropna().between(0, 1).all():
            raise ValueError(f'Original {col} contains invalid probabilities.')
    hit = table.padj.lt(.05) & table.log2FoldChange.abs().ge(1)
    return {'rows': len(table), 'duplicate_gene_ids': int(table.gene.duplicated().sum()),
            'missing_pvalue': int(table.pvalue.isna().sum()),
            'missing_padj': int(table.padj.isna().sum()),
            'fdr_below_0_05': int(table.padj.lt(.05).sum()),
            'fdr_and_abs_lfc_at_least_1': int(hit.sum()),
            'positive_lfc_hits': int((hit & table.log2FoldChange.gt(0)).sum()),
            'negative_lfc_hits': int((hit & table.log2FoldChange.lt(0)).sum()),
            'missing_legacy_flag': int(table.is_significant.isna().sum()),
            'flag_yes_disagreement': int((table.is_significant.eq('Yes') != hit).sum()),
            'status': 'Table audit only; original model cannot be reproduced without original counts.'}


def verify_preservation(root: Path = ROOT) -> int:
    manifest = json.loads((root / 'legacy/original_manifest.json').read_text())
    for entry in manifest['files']:
        # Git normalizes text line endings on checkout; compare canonical original bytes.
        raw = (root / entry['preserved_path']).read_bytes()
        if entry['preserved_path'].endswith(('.sh', '.r', '.R', '.md', '.csv')) or entry['original_path'] in {'.gitattributes', 'LICENSE'}:
            raw = raw.replace(b'\r\n', b'\n')
        if hashlib.sha256(raw).hexdigest() != entry['sha256']:
            raise ValueError(f'Original file changed: {entry["preserved_path"]}')
    return len(manifest['files'])
