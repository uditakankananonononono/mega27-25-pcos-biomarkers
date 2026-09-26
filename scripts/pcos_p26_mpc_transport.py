"""Registered P26 untreated iPSC-MPC PCOS fixed-panel transport."""
import csv,hashlib,json,sys,urllib.request
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p25';folder.mkdir(parents=True,exist_ok=True)
sources=[(folder/'GSE267287_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267287/matrix/GSE267287_series_matrix.txt.gz','93335008afe200eeb1f3f7edf7d77a996bc4ef69179763e9c441fe0fb6b5108f'),(folder/'GSE267287_PCOS_MPC_readcount_table.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267287/suppl/GSE267287_PCOS_MPC_readcount_table.txt.gz','38caa96acb45c5494e6c68af5d0667b750adcfcdc438e10257be726f9cc84ba9')]
for p,url,sha in sources:
 if not p.exists():urllib.request.urlretrieve(url,p)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
_,ann,_=geo.parse_series_matrix(str(sources[0][0]));assert len(ann)==14 and ann.index.is_unique
prior={r['accession'] for r in csv.DictReader((root/'results/dataset_manifest.csv').open())};assert set(ann.index).isdisjoint(prior)
assert ann.Sample_description.is_unique and ann.Sample_characteristics_ch1_2.isin(['disease: Control','disease: PCOS']).all() and ann.Sample_characteristics_ch1_3.isin(['treatment: No','treatment: Yes']).all()
for gsm,row in ann.iterrows():assert row.Sample_title.endswith(' w/ testosterone')==(row.Sample_characteristics_ch1_3=='treatment: Yes')
x=pd.read_csv(sources[1][0],sep='\t',low_memory=False)
assert x.columns[0]=='0' and len(x.columns)==15 and set(x.columns[1:])==set(ann.Sample_description)
assert x['0'].str.fullmatch(r'ENSG\d+\.\d+(?:_PAR_Y)?').all() and x['0'].is_unique
ids=list(ann.index[ann.Sample_characteristics_ch1_3=='treatment: No'])
case=list(ann.index[(ann.Sample_characteristics_ch1_3=='treatment: No')&(ann.Sample_characteristics_ch1_2=='disease: PCOS')]);ctrl=list(ann.index[(ann.Sample_characteristics_ch1_3=='treatment: No')&(ann.Sample_characteristics_ch1_2=='disease: Control')]);assert len(case)==len(ctrl)==4 and len(ids)==8
assert {ann.loc[gsm,'Sample_title'] for gsm in case}=={'P2','P5','P11','P12'}
assert {ann.loc[gsm,'Sample_title'] for gsm in ctrl}=={'C25','C24','C31','C32'}
values=x.iloc[:,1:].to_numpy(dtype=float)
assert np.isfinite(values).all() and (values>=0).all() and np.equal(values,np.floor(values)).all()
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna();amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)]
assert h.ensembl_gene_id.is_unique
mapping=dict(zip(h.ensembl_gene_id,h.symbol));symbols=x['0'].str.split('.').str[0].map(mapping) # _PAR_Y rows share a stable ENSG and sum at gene level
expr=pd.DataFrame(values,columns=[ann.index[ann.Sample_description==col][0] for col in x.columns[1:]])
expr['gene']=symbols.values;expr=expr.dropna(subset=['gene']).groupby('gene').sum()
assert len(expr)>10000
lib=expr.sum(axis=0);assert (lib>0).all()
cpm=expr.div(lib,axis=1)*1e6
keep=(cpm[ids]>1).mean(axis=1)>=.2
log=np.log2(cpm.loc[keep]+1)
g,v=stats.hedges_g(log[case].to_numpy(float),log[ctrl].to_numpy(float))
effect=pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna()
D=pd.read_csv(root/'results/meta_discovery/pcos.csv.gz',index_col=0)
C=pd.read_csv(root/'results/pcos_exploratory_candidates.csv').set_index('gene');assert len(C)==20
obs=C.index.intersection(effect.index)
cp=int((C.loc[obs,'mu']>0).sum());cn=len(obs)-cp
agreement=np.sign(C.loc[obs,'mu'])==np.sign(effect.loc[obs,'g'])
common=D.index.intersection(effect.index);pool=common[(D.loc[common,'k']>=4)&~common.isin(C.index)];up=np.array(pool[D.loc[pool,'mu']>0]);down=np.array(pool[D.loc[pool,'mu']<0]);assert len(up)>=cp and len(down)>=cn
rng=np.random.default_rng(20260925);null=np.zeros(10000,dtype=int)
for i in range(len(null)):
 picks=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
 null[i]=int((np.sign(D.loc[picks,'mu'])==np.sign(effect.loc[picks,'g'])).sum())
count=int(agreement.sum());p=float((1+(null>=count).sum())/(len(null)+1))
with (root/'results/pcos_p26_sample_map.csv').open('w',newline='') as f:
 fields=['gsm','title','sample_description','disease','treatment','primary'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows([dict(gsm=id,title=ann.loc[id,'Sample_title'],sample_description=ann.loc[id,'Sample_description'],disease=ann.loc[id,'Sample_characteristics_ch1_2'],treatment=ann.loc[id,'Sample_characteristics_ch1_3'],primary=int(id in ids)) for id in ann.index])
with (root/'results/pcos_p26_genes.csv').open('w',newline='') as f:
 fields=['gene','discovery_mu','cohort_g','cohort_v','agrees'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows([dict(gene=gene,discovery_mu=float(C.loc[gene,'mu']),cohort_g=float(effect.loc[gene,'g']),cohort_v=float(effect.loc[gene,'v']),agrees=int(agreement.loc[gene])) for gene in obs])
pd.DataFrame({'matched_random_agreements':null}).to_csv(root/'results/pcos_p26_null.csv.gz',index=False)
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267287',matrix_sha256=sources[0][2],count_sha256=sources[1][2],n_untreated_pcos=4,n_untreated_control=4,n_treated_excluded=6,n_unique_mapped_genes=len(expr),n_effect_genes=len(effect),n_fixed_measured=len(obs),n_positive=cp,n_negative=cn,n_sign_match=count,null_mean=float(null.mean()),empirical_p=p,registered_descriptive_support=bool(len(obs)>=15 and count>=15 and p<.025),ambiguous_ensembl_ids=len(amb))
(root/'results/pcos_p26_result.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
