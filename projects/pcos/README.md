# PCOS - proposed separate project

Status: project scaffold and gate audit only, not a finished paper or a positive benchmark. The parent repository is shared core; this directory must hold a disease-specific protocol, accession and external-service evidence ledger, reproducible results and a 50-page substantive paper before its gates can be claimed.

Current disease-tagged manifest records: 54 = 17 GSE studies + 36 nested GSM samples + 1 other. These are record units, not independent datasets or patients. The shared 40-service and 49-page PDF do not transfer as automatic per-project passes. Benchmark/discovery endpoint is open.

Disease-specific documents and source logs are not yet split from shared core; use the shared source code and result filenames by disease as leads, then verify original record attribution before copying.

## GSE5090 duplicate donor arrays

The old pooled discovery calls [GSE5090](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5090) nine PCOS versus eight control adipose **arrays**, but the series description and 17 individual sample titles bind two PCOS arrays to EP1 and two control arrays to EP31: eight PCOS versus seven control **people**. We fetched every GSM and checked title, GPL96 platform, source phenotype and expression column; `sources/GSE5090_used_sample_crosswalk.csv` records source URLs and hashes. A post-audit donor-mean effect sensitivity (`scripts/pcos_gse5090_donor_sensitivity.py`) averages each duplicated EP pair before the 8-versus-7 contrast. Across 13,047 shared genes, 1,346 effect signs change and median absolute Hedges-effect shift is 0.114. This is a correction to the independent statistical unit, not a new cohort, a new PCOS marker, or a preregistered validation. The earlier pooled meta and downstream CNN/replication were not rerun and cannot be cited as donor-corrected. The 17 historically analyzed GSM records remain outside the used-accession tally pending builder branch integration to avoid simultaneous manifest edits; they are nested records, not 17 people or studies.
