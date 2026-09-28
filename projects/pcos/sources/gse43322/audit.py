"""Audit historically pooled PCOS SuperSeries, preserving the old 24-vs-7 effect.

Exploratory baseline-only sensitivity; no independent validation or novel marker.
"""
import csv
import gzip
import hashlib
import io
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
from ubiomark.geo import parse_series_matrix,to_gene_level
from ubiomark.stats import hedges_g

SERIES='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43322'
MATRIX_HASH='b388e83c3715c6d4df4bf1e299569f0a9e84d9b877ed574804a0250f225cd160'


def run():
    path=ROOT/'data/geo/gse43322/GSE43322_series_matrix.txt.gz'
    platform=ROOT/'data/geo/gse43322/GPL15362_table.txt'
    if hashlib.sha256(path.read_bytes()).hexdigest()!=MATRIX_HASH:
        raise ValueError('series matrix changed')
    matrix,ann,_=parse_series_matrix(str(path))
    if matrix.shape!=(17126,31) or set(ann.Sample_platform_id)!={'GPL15362'}:
        raise ValueError('unexpected matrix/platform')
    table=pd.read_csv(io.StringIO(''.join(l for l in open(platform) if not l.startswith(('^','!','#')))),sep='\t',dtype=str)
    mapping=table.set_index('ID').ORF.dropna()
    gene=to_gene_level(matrix,mapping)
    labels=list(csv.DictReader(open(ROOT/'results/series/pcos__GSE43322.labels.csv')))
    old={r['gsm']:r['label'] for r in labels}
    if set(old)!=set(ann.index) or len(old)!=31:
        raise ValueError('old labels and source GSMs differ')
    rows=[]; baseline_case=[]; placebo_case=[]; treatment_case=[]; control=[]
    for gsm, a in ann.iterrows():
        title=a.Sample_title
        if re.fullmatch(r'Adipose tissue; PCOS_\d+',title) and old[gsm]=='case':
            group='baseline_PCOS';baseline_case.append(gsm)
        elif re.fullmatch(r'Adipose tissue; PCOS; placebo_\d+',title) and old[gsm]=='case':
            group='placebo_PCOS';placebo_case.append(gsm)
        elif re.fullmatch(r'Adipose tissue; PCOS; LC n-3 PUFA_\d+',title) and old[gsm]=='case':
            group='PUFA_PCOS';treatment_case.append(gsm)
        elif re.fullmatch(r'Adipose tissue; control_\d+',title,re.I) and old[gsm]=='control':
            group='healthy_control';control.append(gsm)
        else:
            raise ValueError(f'unexpected source label {gsm}: {title}/{old[gsm]}')
        url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
        file=HERE/f'{gsm}.soft.txt.gz'
        if not file.exists():
            with urllib.request.urlopen(url,timeout=35) as r:raw=r.read()
            if len(raw)<25000:raise ValueError(f'{gsm} short source record')
            with gzip.open(file,'wb',compresslevel=9) as f:f.write(raw)
        else:raw=gzip.open(file,'rb').read()
        text=raw.decode().replace('\r\n','\n')
        for line in [f'^SAMPLE = {gsm}',f'!Sample_geo_accession = {gsm}',
                     f'!Sample_title = {title}','!Sample_organism_ch1 = Homo sapiens',
                     '!Sample_platform_id = GPL15362','!Sample_series_id = GSE43322',
                     '!Sample_characteristics_ch1 = tissue: subcutaneous adipose tissue']:
            if line not in text.splitlines():raise ValueError(f'{gsm} source mismatch: {line}')
        if group in ('placebo_PCOS','PUFA_PCOS'):
            agent='placebo' if group=='placebo_PCOS' else 'LC n-3 PUFA'
            if '!Sample_series_id = GSE43266' not in text.splitlines() or f'!Sample_characteristics_ch1 = agent: {agent}' not in text.splitlines():
                raise ValueError(f'{gsm} subseries/intervention mismatch')
        else:
            if '!Sample_series_id = GSE43264' not in text.splitlines() or '!Sample_characteristics_ch1 = agent: placebo' in text.splitlines():
                raise ValueError(f'{gsm} baseline subseries mismatch')
        start=text.index('!sample_table_begin\n')+len('!sample_table_begin\n')
        end=text.index('!sample_table_end',start)
        t=pd.read_csv(io.StringIO(text[start:end]),sep='\t',index_col=0)
        t.index=t.index.astype(str)
        if len(t)!=17126 or set(t.index)!=set(matrix.index.astype(str)):
            raise ValueError(f'{gsm} probe mismatch')
        if not np.allclose(t.loc[matrix.index.astype(str),'VALUE'].to_numpy(float),
                           matrix[gsm].to_numpy(float),atol=1e-5,rtol=1e-6,equal_nan=True):
            raise ValueError(f'{gsm} GEO full record and matrix differ')
        rows.append({'gsm':gsm,'source_url':url,'source_sha256':hashlib.sha256(raw).hexdigest(),
                     'title':title,'old_label':old[gsm],'source_group':group,
                     'subseries':'GSE43266' if group in ('placebo_PCOS','PUFA_PCOS') else 'GSE43264',
                     'platform':'GPL15362','probes_verified':len(t),'matrix_sha256':MATRIX_HASH})
    if tuple(map(len,(baseline_case,placebo_case,treatment_case,control)))!=(8,8,8,7):
        raise ValueError('unexpected source strata')
    pd.DataFrame(rows).to_csv(HERE/'GSE43322_used_sample_crosswalk.csv',index=False)
    old_g,old_v=hedges_g(gene[baseline_case+placebo_case+treatment_case].to_numpy(float),gene[control].to_numpy(float))
    retained=pd.read_csv(ROOT/'results/series/pcos__GSE43322.csv.gz').set_index('gene').loc[gene.index]
    if np.nanmax(np.abs(old_g-retained.g.to_numpy(float)))>1e-4:
        raise ValueError('old effects do not replay from source')
    strict_g,strict_v=hedges_g(gene[baseline_case].to_numpy(float),gene[control].to_numpy(float))
    common=np.isfinite(old_g)&np.isfinite(strict_g)
    flips=common&(np.sign(old_g)!=np.sign(strict_g))
    outcomes=pd.DataFrame({'gene':gene.index,'old_pooled_g':old_g,'old_pooled_v':old_v,
                           'baseline_only_g':strict_g,'baseline_only_v':strict_v,'sign_flip':flips})
    outcomes.to_csv(HERE/'GSE43322_baseline_scope_gene_effects.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    result={'source':SERIES,'historical_pooled_case_arrays':24,'baseline_PCOS_arrays':8,
            'placebo_PCOS_arrays':8,'PUFA_PCOS_arrays':8,'control_arrays':7,'genes_common':int(common.sum()),
            'effect_sign_flips':int(flips.sum()),
            'old_effect_max_abs_replay_delta':float(np.nanmax(np.abs(old_g-retained.g.to_numpy(float)))),
            'caveat':'Post-exposure baseline-only case restriction in same SuperSeries. Eight placebo and eight LC n-3 PUFA biopsies are from a PCOS crossover intervention subseries; unknown participant overlap with eight baseline cases. No baseline-vs-placebo person matching, untouched validation, clinical marker, or published comparator win.'}
    (HERE/'GSE43322_audit_result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':print(json.dumps(run(),indent=2))
