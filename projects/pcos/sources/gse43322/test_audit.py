import csv
import gzip
import hashlib
from collections import Counter
from pathlib import Path
from audit import run,MATRIX_HASH

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def test_source_scope_and_replay():
    result=run()
    assert (result['baseline_PCOS_arrays'],result['placebo_PCOS_arrays'],result['PUFA_PCOS_arrays'],result['control_arrays'])==(8,8,8,7)
    assert result['genes_common']==16910 and result['old_effect_max_abs_replay_delta']<.0001
    cross=list(csv.DictReader(open(HERE/'GSE43322_used_sample_crosswalk.csv')))
    assert Counter(x['source_group'] for x in cross)=={'baseline_PCOS':8,'placebo_PCOS':8,'PUFA_PCOS':8,'healthy_control':7}
    assert len({x['gsm'] for x in cross})==31
    for x in cross:
        assert hashlib.sha256(gzip.open(HERE/f'{x["gsm"]}.soft.txt.gz','rb').read()).hexdigest()==x['source_sha256']
        assert x['matrix_sha256']==MATRIX_HASH
    for m in ['disease_tagged_accession_manifest.csv']:
        rows=list(csv.DictReader(open(ROOT/'results'/m)))
        assert len(rows)==len({x['accession'] for x in rows})==121
        assert {x['gsm'] for x in cross}<={x['accession'] for x in rows}
        assert {'GSE43264','GSE43266'}<={x['accession'] for x in rows}


def test_nested_series_membership_not_new_cohorts():
    import re
    cross=list(csv.DictReader(open(HERE/'GSE43322_used_sample_crosswalk.csv')))
    for series,n in [('GSE43264',15),('GSE43266',16)]:
        text=(HERE/f'{series}.series.soft.txt').read_text()
        ids=re.findall(r'^!Series_sample_id = (GSM\d+)',text,re.M)
        assert len(ids)==n and set(ids)=={x['gsm'] for x in cross if x['subseries']==series}
        assert f'!Series_relation = SubSeries of: GSE43322' in text
