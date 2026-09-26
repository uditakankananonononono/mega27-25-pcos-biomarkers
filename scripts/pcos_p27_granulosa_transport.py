"""P27 pre-registered PCOS granulosa fixed-sign transport; TPM as deposited."""
import csv,gzip,hashlib,io,json,re,sys,tarfile,urllib.request
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,pandas as pd
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p27';folder.mkdir(parents=True,exist_ok=True)
sources=[('GSE294074_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294074/matrix/GSE294074_series_matrix.txt.gz','575fef37c3efe5404230e4dd041e51d07180c27c265e504fa366b6a6d233e548'),('GSE294074_RAW.tar','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294074/suppl/GSE294074_RAW.tar','aced9e17fed2f49b1cccf48d2463d46a4796322744920957f1499ee143496156'),('GSE294074_annotation.xls.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294074/suppl/GSE294074_annotation.xls.gz','f2ab5c472dfbe617d4d487d949f351f64c48d45c900635bfed1a442e33c568d7')]
for name,url,sha in sources:
 p=folder/name
 if not p.exists():urllib.request.urlretrieve(url,p)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
_,ann,_=geo.parse_series_matrix(str(folder/sources[0][0]));assert len(ann)==6 and ann.index.is_unique
case=['GSM8898322','GSM8898323','GSM8898324'];ctrl=['GSM8898319','GSM8898320','GSM8898321'];ids=ctrl+case
assert set(ann.index)==set(ids) and list(ann.loc[ctrl,'Sample_title'])==['Ctrl-1','Ctrl-2','Ctrl-3'] and list(ann.loc[case,'Sample_title'])==['PCOS-1','PCOS-2','PCOS-3']
prior={r['accession'] for r in csv.DictReader((root/'results/dataset_manifest.csv').open())};assert set(ids).isdisjoint(prior)
# Only unique human SwissProt GN assignments; do not guess from functional descriptions.
protein={};ambiguous=set();annotation_ids=set()
with gzip.open(folder/sources[2][0],'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  key=row['Unigene'].split('(')[0];annotation_ids.add(key);tokens=row['Swissprot'].split(' &gt;')
  genes=set()
  for t in tokens:
   if 'OS=Homo sapiens' in t:
    genes.update(re.findall(r'\bGN=([A-Za-z0-9_-]+)',t))
  if len(genes)==1:
   if key in protein and protein[key]!=next(iter(genes)):ambiguous.add(key)
   protein[key]=next(iter(genes))
  elif len(genes)>1:ambiguous.add(key)
for key in ambiguous:protein.pop(key,None)
assert len(protein)>50000
sample={}
with tarfile.open(folder/sources[1][0]) as tar:
 names=tar.getnames();assert len(names)==6
 for gsm in ids:
  matched=[n for n in names if n.startswith(gsm+'_')];assert len(matched)==1
  with gzip.open(tar.extractfile(matched[0]),'rt') as f:
   x=pd.read_csv(f,sep='\t',usecols=['gene','tpm']);assert x.gene.is_unique and x.tpm.notna().all() and (x.tpm>=0).all()
   sample[gsm]=x.set_index('gene').tpm
expr=pd.DataFrame(sample);assert set(expr.columns)==set(ids) and set(expr.index).issubset(annotation_ids)
# Per-sample PB transcript lists differ greatly. Zero-fill missing table entries solely for deposited-TPM aggregation; this does not establish biological absence or comparable transcript identities.
missing_cells=int(expr.isna().sum().sum());expr=expr.fillna(0)
expr['symbol']=expr.index.map(protein);expr=expr.dropna(subset=['symbol']).groupby('symbol').sum()
expr=expr.loc[(expr[ids]>1).sum(axis=1)>=2];log=np.log2(expr[ids]+1)
g,v=stats.hedges_g(log[case].to_numpy(float),log[ctrl].to_numpy(float));effect=pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna()
D=pd.read_csv(root/'results/meta_discovery/pcos.csv.gz',index_col=0);C=pd.read_csv(root/'results/pcos_exploratory_candidates.csv').set_index('gene');assert len(C)==20
obs=C.index.intersection(effect.index);cp=int((C.loc[obs,'mu']>0).sum());cn=len(obs)-cp
agreement=np.sign(C.loc[obs,'mu'])==np.sign(effect.loc[obs,'g'])
common=D.index.intersection(effect.index);pool=common[(D.loc[common,'k']>=4)&~common.isin(C.index)];up=np.array(pool[D.loc[pool,'mu']>0]);down=np.array(pool[D.loc[pool,'mu']<0]);assert len(up)>=cp and len(down)>=cn
rng=np.random.default_rng(20260925);null=np.zeros(10000,dtype=int)
for i in range(len(null)):
 picks=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
 null[i]=int((np.sign(D.loc[picks,'mu'])==np.sign(effect.loc[picks,'g'])).sum())
count=int(agreement.sum());p=float((1+(null>=count).sum())/(len(null)+1))
with (root/'results/pcos_p27_genes.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['gene','discovery_mu','cohort_g','cohort_v','agrees','status']);w.writeheader()
 for gene in C.index:w.writerow(dict(gene=gene,discovery_mu=float(C.loc[gene,'mu']),cohort_g=float(effect.loc[gene,'g']) if gene in obs else '',cohort_v=float(effect.loc[gene,'v']) if gene in obs else '',agrees=int(agreement.loc[gene]) if gene in obs else '',status='measured' if gene in obs else 'not uniquely mapped or filtered'))
pd.DataFrame({'matched_random_agreements':null}).to_csv(root/'results/pcos_p27_null.csv.gz',index=False)
rec=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294074',sha256={n:s for n,_,s in sources},n_pcos=3,n_control=3,n_mapped_transcripts=len(protein),n_omitted_transcript_sample_cells=missing_cells,n_filtered_genes=len(expr),n_effect_genes=len(effect),n_fixed_measured=len(obs),missing_fixed=list(C.index.difference(obs)),n_positive=cp,n_negative=cn,n_sign_match=count,null_mean=float(null.mean()),empirical_p=p,registered_descriptive_support=bool(len(obs)>=15 and p<.025),ambiguous_annotation_ids=len(ambiguous))
(root/'results/pcos_p27_result.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
