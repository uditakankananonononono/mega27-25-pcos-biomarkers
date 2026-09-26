"""Pre-registered PCOS primary and secondary fresh cohorts against a matched null.

The cohort choices and threshold were registered before processing. Targeted genes
are mostly discovery-positive; 0.5 binomial null can be anti-conservative.
Each random set samples without replacement from genes shared with the cohort,
with k>=4 discovery cohorts, matching the 20-gene target's positive/negative
sign counts separately. Reports both cohorts and a shared-gene pooled test.
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import binomtest
sys.path.insert(0,'src')
from ubiomark import stats
rng=np.random.default_rng(20260925)
B=10000; D=pd.read_csv('results/meta_discovery/pcos.csv.gz',index_col=0)
C=pd.read_csv('results/pcos_exploratory_candidates.csv').set_index('gene')
V={g:pd.read_csv(f'results/series/pcos__{g}.csv.gz').set_index('gene') for g in ['GSE155489','GSE262735']}
rows=[]
for g,v in V.items():
    test=C.index.intersection(v.index)
    cp=int(sum(C.loc[test,'mu']>0));cn=len(test)-cp
    obs=int(sum(np.sign(C.loc[test,'mu'])==np.sign(v.loc[test,'g'])))
    candidates=D.index.intersection(v.index)
    candidates=candidates[(D.loc[candidates,'k']>=4) & ~candidates.isin(C.index)]
    up=np.array(candidates[D.loc[candidates,'mu']>0]); down=np.array(candidates[D.loc[candidates,'mu']<0])
    null=np.zeros(B,int)
    for i in range(B):
        genes=list(rng.choice(up,cp,replace=False))+list(rng.choice(down,cn,replace=False))
        null[i]=int(sum(np.sign(D.loc[genes,'mu'])==np.sign(v.loc[genes,'g'])))
    rows.append(dict(gse=g,n_genes=len(test),n_positive=cp,n_negative=cn,n_agree=obs,
                     p_binomial_half=binomtest(obs,len(test),.5,alternative='greater').pvalue,
                     null_mean=float(null.mean()),null_sd=float(null.std()),
                     p_direction_matched=(1+int(sum(null>=obs)))/(B+1)))
pd.DataFrame(rows).to_csv('results/pcos_prospective_matched_null.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
