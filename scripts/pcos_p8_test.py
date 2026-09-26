"""Audited GSE293353 follicular-granulosa PCOS sign-panel test registered in P8."""
import sys,gzip,io,urllib.request,hashlib,json,os,glob
import numpy as np,pandas as pd
sys.path.insert(0,'src')
from ubiomark import geo,stats
G='GSE293353';base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{G[:-3]}nnn/{G}/suppl/'
matfile='GSE293353_gene_count_matrix.txt.gz';clinfile='GSE293353_Sample_Information.csv.gz'
_,ann,_=geo.parse_series_matrix(geo.download_matrices(G)[0]); assert len(ann)==18
sample_info=gzip.decompress(urllib.request.urlopen(base+clinfile,timeout=20).read());clinic=pd.read_csv(io.BytesIO(sample_info))
assert len(clinic)==18 and clinic.group.is_unique
names=ann.Sample_title.str.extract(r'Follicular granulosa cells, (PCOS\d+|Control\d+)$',expand=False)
assert names.notna().all() and names.is_unique and set(names)==set(clinic.group)
path=geo._fetch(base+matfile,f'data/raw/rnaseq/{G}/{matfile}')
x=pd.read_csv(path,sep='\t',low_memory=False)
assert len(x.columns)==19 and set(x.columns)==set(clinic.group)|{'gene_id'}
assert x.gene_id.is_unique and not x.gene_id.isna().any()
assert x.gene_id.astype(str).str.match(r'^ENSG\d+(\.\d+)?$').all()
old=set()
for f in glob.glob('results/series/pcos__*.labels.csv'):
 old|=set(pd.read_csv(f).gsm.astype(str))
assert not(set(ann.index)&old)
# Map frozen Ensembl stable IDs using the HGNC complete set used by the original pipeline.
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna()
ambiguous=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id'])
# Do not assign one of two HGNC symbols to an ambiguous Ensembl ID.
h=h[~h.ensembl_gene_id.isin(ambiguous)]
assert h.ensembl_gene_id.is_unique
mapping=dict(zip(h.ensembl_gene_id,h.symbol))
x['gene']=x.gene_id.str.split('.').str[0].map(mapping)
n_unmapped=int(x.gene.isna().sum());x=x.dropna(subset=['gene']).drop(columns='gene_id').set_index('gene')
x=x.apply(pd.to_numeric,errors='raise')
assert (x.to_numpy()>=0).all()
x=x.groupby(level=0).sum()
cols=list(x.columns);cpm=x.div(x.sum(axis=0),axis=1)*1e6
x=np.log2(cpm.loc[(cpm>1).mean(axis=1)>=.2]+1)
case=[c for c in cols if c.startswith('PCOS')];ctrl=[c for c in cols if c.startswith('Control')]
assert len(case)==len(ctrl)==9
g,v=stats.hedges_g(x[case].to_numpy(float),x[ctrl].to_numpy(float))
E=pd.DataFrame({'g':g,'v':v},index=x.index).dropna()
D=pd.read_csv('results/meta_discovery/pcos.csv.gz',index_col=0)
C=pd.read_csv('results/pcos_exploratory_candidates.csv').set_index('gene')
obs=C.index.intersection(E.index)
cp=int(sum(C.loc[obs,'mu']>0));cn=len(obs)-cp
agreements=(np.sign(C.loc[obs,'mu'])==np.sign(E.loc[obs,'g']))
common=D.index.intersection(E.index)
pool=common[(D.loc[common,'k']>=4)&~common.isin(C.index)]
up=np.array(pool[D.loc[pool,'mu']>0]);down=np.array(pool[D.loc[pool,'mu']<0])
assert len(up)>=cp and len(down)>=cn
rng=np.random.default_rng(20260925);null=np.zeros(10000,int)
for i in range(len(null)):
 picks=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
 null[i]=int(sum(np.sign(D.loc[picks,'mu'])==np.sign(E.loc[picks,'g'])))
pd.DataFrame({'matched_random_agreements':null}).to_csv('results/pcos_p8_null_draws.csv.gz',index=False)
count=int(agreements.sum());p=(1+int(sum(null>=count)))/(len(null)+1)
row=dict(gse=G,n_case=9,n_control=9,n_genes_measured=len(obs),n_positive=cp,n_negative=cn,n_agree=count,fraction_agree=float(count/len(obs)),null_mean=float(null.mean()),null_sd=float(null.std(ddof=1)),empirical_p=p,passes_per_cohort=(len(obs)>=15 and p<.025),mapped_genes=len(E),unmapped_ensembl_rows=n_unmapped,ambiguous_hgnc_ensembl_ids=sorted(ambiguous),source=base+matfile,matrix_sha256=hashlib.sha256(open(path,'rb').read()).hexdigest(),sample_info_sha256=hashlib.sha256(sample_info).hexdigest())
with open('results/pcos_p8.json','w') as f:json.dump(row,f,indent=2)
pd.DataFrame({'gsm':ann.index,'group':names.to_numpy(),'label':np.where(names.str.startswith('PCOS'),'case','control')}).to_csv('results/pcos_p8_labels.csv',index=False)
pd.DataFrame({'gene':obs,'discovery_mu':C.loc[obs,'mu'].to_numpy(),'cohort_g':E.loc[obs,'g'].to_numpy(),'cohort_v':E.loc[obs,'v'].to_numpy(),'agrees':agreements.astype(int).to_numpy()}).to_csv('results/pcos_p8_genes.csv',index=False)
E.reset_index(names='gene').to_csv(f'results/series/pcos__{G}.csv.gz',index=False,float_format='%.6g')
print(json.dumps(row,indent=2))
