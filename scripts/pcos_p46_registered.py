"""Locked GSE271363 cumulus-cell sign test; protocol at projects/judge_rounds/2026-09-27_pcos_new_cohort_registration.md."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / 'data/geo/p46/GSE271363_series_matrix.txt.gz'
DATA = ROOT / 'data/geo/p46/GSE271363_RPKM_PCOS_CONTROL.txt.gz'
SOURCE_SHA = 'c260e9e1a4465cf6cae6cacd1e31223acd2089f9227042c0dcaf202927da94f6'
META_SHA = 'fc3d45b681877e46b5b9f8d63248a5914f0eec8ca490ffc0768ddf18befe95c0'

def field(lines,key):
    entries=[row for row in csv.reader(lines,delimiter='\t') if row and row[0]==key]
    assert len(entries)==1,key
    return entries[0][1:]

def run():
    assert hashlib.sha256(META.read_bytes()).hexdigest()==META_SHA
    assert hashlib.sha256(DATA.read_bytes()).hexdigest()==SOURCE_SHA
    with gzip.open(META,'rt') as f:lines=f.readlines()
    gsm=field(lines,'!Sample_geo_accession');desc=field(lines,'!Sample_description');title=field(lines,'!Sample_title')
    assert len(gsm)==len(desc)==len(title)==16 and len(set(gsm))==16 and len(set(desc))==16
    assert all((d.startswith('PCOS-'))==(', PCOS' in t) and (d.startswith('CONTROL-'))==(', CONTROL' in t) for d,t in zip(desc,title))
    manifest={r['accession']:r for r in csv.DictReader((ROOT/'results/disease_tagged_accession_manifest.csv').open())}
    present=set(gsm)&set(manifest)
    assert not present or (present==set(gsm) and manifest.get('GSE271363') and all('P46' in manifest[g]['role'] and 'GSE271363' in manifest[g]['role'] for g in gsm)), 'unexpected source overlap'
    x=pd.read_csv(DATA,sep='\t',dtype=str)
    assert {'ENS','GENE_NAME'} <= set(x.columns)
    assert set(x.columns)-{'ENS','GENE_NAME'}==set(desc)
    assert len(x)>10000 and x.ENS.is_unique
    assert x.ENS.str.fullmatch(r'ENSG\d+(\.\d+)?').all()
    # P46 source is explicitly RPKM; decimal commas are in the raw text.
    numeric=x[desc].replace(',','.',regex=True).apply(pd.to_numeric,errors='raise')
    assert np.isfinite(numeric.to_numpy()).all() and (numeric.to_numpy()>=0).all()
    from ubiomark import geo
    h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna()
    ambiguous=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id'])
    h=h.loc[~h.ensembl_gene_id.isin(ambiguous)]
    sym=dict(zip(h.ensembl_gene_id,h.symbol))
    x['gene']=x.ENS.str.split('.').str[0].map(sym)
    valid=x.gene.notna();expr=numeric.loc[valid].copy();expr.index=x.loc[valid,'gene']
    expr=expr.groupby(level=0,sort=True).mean()
    log=np.log2(expr+1)
    case=[d for d in desc if d.startswith('PCOS-')];control=[d for d in desc if d.startswith('CONTROL-')]
    assert len(case)==5 and len(control)==11
    from ubiomark.stats import hedges_g
    g,v=hedges_g(log[case].to_numpy(float),log[control].to_numpy(float))
    effect=pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna()
    panel=pd.read_csv(ROOT/'results/pcos_exploratory_candidates.csv').head(20).set_index('gene')
    assert len(panel)==20 and panel.index.is_unique
    observed=panel.index.intersection(effect.index,sort=False)
    if len(observed)<15:raise RuntimeError(f'coverage below registered minimum: {len(observed)}/20')
    agreed=(np.sign(panel.loc[observed,'mu'])==np.sign(effect.loc[observed,'g']))
    D=pd.read_csv(ROOT/'results/meta_discovery/pcos.csv.gz',index_col=0)
    common=D.index.intersection(effect.index)
    pool=common[(D.loc[common,'k']>=4)&~common.isin(panel.index)]
    up=np.array(pool[D.loc[pool,'mu']>0]);down=np.array(pool[D.loc[pool,'mu']<0])
    cp=int((panel.loc[observed,'mu']>0).sum());cn=int((panel.loc[observed,'mu']<0).sum())
    assert cp+cn==len(observed) and len(up)>=cp and len(down)>=cn
    rng=np.random.default_rng(20260927)
    null=np.empty(10000,dtype=int)
    for i in range(len(null)):
        picks=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
        null[i]=int((np.sign(D.loc[picks,'mu'])==np.sign(effect.loc[picks,'g'])).sum())
    count=int(agreed.sum());p=(1+int((null>=count).sum()))/(len(null)+1)
    result={'gse':'GSE271363','source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271363','data_sha256':SOURCE_SHA,'meta_sha256':META_SHA,
        'n_case':5,'n_control':11,'n_panel_measured':len(observed),'n_positive':cp,'n_negative':cn,'n_agree':count,
        'null_mean':float(null.mean()),'empirical_p':p,'registered_support':bool(len(observed)>=15 and p<.025),
        'participant_identity':'No independently verified donor crosswalk; 16 distinct GSM/BioSample/title records only',
        'interpretation':'Exploratory frozen panel transport; not clinical prediction, new biomarker or published benchmark beat.'}
    out=ROOT/'results';(out/'pcos_p46_result.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame({'gene':observed,'discovery_mu':panel.loc[observed,'mu'].to_numpy(), 'test_g':effect.loc[observed,'g'].to_numpy(),'agrees':agreed.to_numpy()}).to_csv(out/'pcos_p46_genes.csv',index=False)
    pd.DataFrame({'gsm':gsm,'description':desc,'title':title,'label':['case' if d.startswith('PCOS-') else 'control' for d in desc]}).to_csv(out/'pcos_p46_samples.csv',index=False)
    pd.DataFrame({'matched_random_agreements':null}).to_csv(out/'pcos_p46_null.csv.gz',index=False)
    return result

if __name__=='__main__':
    import sys
    sys.path.insert(0,str(ROOT/'src'))
    print(json.dumps(run(),indent=2))
