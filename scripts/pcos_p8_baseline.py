"""Exploratory matched-task benchmark of a frozen 19-gene signed score on P8 cases."""
import json,sys,os
import numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score
from sklearn.utils import resample
sys.path.insert(0,'src')
from ubiomark import geo
G='GSE293353';f='data/raw/rnaseq/GSE293353/GSE293353_gene_count_matrix.txt.gz'
x=pd.read_csv(f,sep='\t',low_memory=False)
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna()
amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)]
x['gene']=x.gene_id.str.split('.').str[0].map(dict(zip(h.ensembl_gene_id,h.symbol)))
x=x.dropna(subset=['gene']).drop(columns='gene_id').set_index('gene').apply(pd.to_numeric).groupby(level=0).sum()
cpm=x.div(x.sum(0),axis=1)*1e6
z=np.log2(cpm.loc[(cpm>1).mean(axis=1)>=.2]+1)
panel=pd.read_csv('results/pcos_exploratory_candidates.csv').set_index('gene'); genes=panel.index.intersection(z.index)
assert len(genes)==19
# Within-cohort z transform is intentionally descriptive, not a deployable predictor.
X=z.loc[genes];means=X.mean(axis=1);sd=X.std(axis=1,ddof=1).replace(0,np.nan)
score=((X.sub(means,axis=0).div(sd,axis=0)).mul(np.sign(panel.loc[genes,'mu']),axis=0)).mean(axis=0)
y=np.array([int(c.startswith('PCOS')) for c in score.index]);assert sum(y)==9
auroc=roc_auc_score(y,score);auprc=average_precision_score(y,score)
# Fixed random-sign baseline on the same selected genes, exploratory in this already seen cohort.
rng=np.random.default_rng(20260925)
null=np.array([roc_auc_score(y,(X.sub(means,axis=0).div(sd,axis=0)).mul(rng.choice([-1,1],len(genes)),axis=0).mean(axis=0)) for i in range(1000)])
out={'gse':G,'n_case':9,'n_control':9,'n_genes':len(genes),'signed_panel_within_cohort_auroc':float(auroc),'signed_panel_within_cohort_auprc':float(auprc),'random_sign_baseline_auroc_mean':float(null.mean()),'random_sign_baseline_auroc_p_ge_observed':float((1+sum(null>=auroc))/(len(null)+1)),'status':'EXPLORATORY after viewing P8 gene signs; within-cohort normalization uses held-out cohort itself. Not registered and not a valid clinical classifier or published benchmark.'}
with open('results/pcos_p8_score_exploratory.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
