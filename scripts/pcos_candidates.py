"""Discovery-only PCOS exploratory shortlist; no validation data is read.

Rule locked before any *new* cohort validation: q<0.2, k>=4, not in the
Open Targets PCOS disease export, rank by |z| and retain first 20. These are
hypothesis-generating because q<0.05 yields no discovery hits.
"""
import pandas as pd
m = pd.read_csv('results/meta_discovery/pcos.csv.gz', index_col=0)
ot = pd.read_csv('data/meta/opentargets_pcos.csv')
c = m[(m.q < 0.2) & (m.k >= 4) & ~m.index.isin(ot.symbol)].copy()
c.index.name = 'gene'
c = c.assign(abs_z=c.z.abs()).sort_values(['abs_z', 'gene'], ascending=[False, True]).head(20)
c.to_csv('results/pcos_exploratory_candidates.csv', float_format='%.6g')
print(c[['mu', 'z', 'q', 'k']].to_string())
