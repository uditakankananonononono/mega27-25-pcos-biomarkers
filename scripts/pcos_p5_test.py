"""Registered P5 matched-sign test in fresh GSE193123 pooled-library cohort."""
import numpy as np
import pandas as pd
from scipy.stats import binomtest
D=pd.read_csv('results/meta_discovery/pcos.csv.gz',index_col=0)
C=pd.read_csv('results/pcos_exploratory_candidates.csv').set_index('gene')
v=pd.read_csv('results/series/pcos__GSE193123.csv.gz').set_index('gene')
test=C.index.intersection(v.index)
cp=int((C.loc[test,'mu']>0).sum());cn=len(test)-cp
obs=int((np.sign(C.loc[test,'mu'])==np.sign(v.loc[test,'g'])).sum())
candidates=D.index.intersection(v.index)
candidates=candidates[(D.loc[candidates,'k']>=4)&~candidates.isin(C.index)]
up=np.array(candidates[D.loc[candidates,'mu']>0]);down=np.array(candidates[D.loc[candidates,'mu']<0])
assert len(up)>=cp and len(down)>=cn
rng=np.random.default_rng(20260925); B=10000
null=np.zeros(B,int)
for i in range(B):
    genes=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
    null[i]=int((np.sign(D.loc[genes,'mu'])==np.sign(v.loc[genes,'g'])).sum())
row=dict(gse='GSE193123',n_genes=len(test),n_positive=cp,n_negative=cn,n_agree=obs,
         p_binomial_half=binomtest(obs,len(test),.5,alternative='greater').pvalue,
         null_mean=float(null.mean()),null_sd=float(null.std()),
         p_direction_matched=(1+int((null>=obs).sum()))/(B+1),
         n_libraries_per_group=3,pooled_individuals_per_library=2)
pd.DataFrame([row]).to_csv('results/pcos_p5_matched_null.csv',index=False)
print(row)
