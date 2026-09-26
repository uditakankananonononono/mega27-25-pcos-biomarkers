"""Post-audit PCOS GSE5090 effects on 8 cases/7 controls, averaging within EP tokens.
Not a new cohort or preregistered validation; historical sample effects are retained.
"""
import re,sys,json,collections
import pandas as pd,numpy as np
sys.path.insert(0,'src')
from ubiomark import geo,stats
G='GSE5090';_,ann,_=geo.parse_series_matrix(geo.download_matrices(G)[0]);X=pd.read_pickle(f'data/processed/pcos__{G}.pkl.gz');lab=pd.read_csv(f'results/series/pcos__{G}.labels.csv').set_index('gsm').label
assert set(X.columns)==set(ann.index)==set(lab.index) and len(X.columns)==17
person={g:re.match(r'(EP\d+)_',ann.loc[g,'Sample_title']).group(1) for g in X.columns}
groups={p:lab.loc[[g for g in X.columns if person[g]==p]].unique().tolist() for p in set(person.values())};assert all(len(v)==1 for v in groups.values())
assert collections.Counter(v[0] for v in groups.values())=={'case':8,'control':7}
assert {p:sum(person[g]==p for g in X.columns) for p in ['EP1','EP31']}=={'EP1':2,'EP31':2}
pooled=pd.DataFrame({p:X[[g for g in X.columns if person[g]==p]].mean(axis=1) for p in sorted(groups)})
cases=[p for p,v in groups.items() if v==['case']];ctrl=[p for p,v in groups.items() if v==['control']]
g,v=stats.hedges_g(pooled[cases].to_numpy(float),pooled[ctrl].to_numpy(float));new=pd.DataFrame({'gene':pooled.index,'g':g,'v':v}).replace([np.inf,-np.inf],np.nan).dropna();new.to_csv('results/pcos_GSE5090_donor_effects.csv.gz',index=False)
old=pd.read_csv('results/series/pcos__GSE5090.csv.gz').set_index('gene');q=new.set_index('gene');ix=old.index.intersection(q.index)
res={'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5090','type':'post-audit donor-level effect-only sensitivity','old_case_arrays':9,'old_control_arrays':8,'n_case_people':8,'n_control_people':7,'replicate_tokens':['EP1','EP31'],'genes_shared':len(ix),'n_sign_changed':int((np.sign(q.loc[ix,'g'])!=np.sign(old.loc[ix,'g'])).sum()),'median_absolute_hedges_g_shift':float((q.loc[ix,'g']-old.loc[ix,'g']).abs().median()),'downstream':'not rerun; old meta, panel and classifier cannot be called donor-corrected'}
with open('results/pcos_GSE5090_donor_sensitivity.json','w') as f:json.dump(res,f,indent=2)
print(json.dumps(res,indent=2))
