"""Certify finite audit coverage, pinned inputs and exported null witnesses.
Run .venv/bin/python scripts/verify_odd_alp_target_nulls.py.
"""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
    path=ROOT/'data/vector_sector/odd_alp_target_nulls.json'
    data=json.loads(path.read_text());rows=data['results']
    assert data['complete_64_scenario_run'] and len(rows)==448
    assert len({(r['scenario_index'],r['model']) for r in rows})==448
    assert {r['scenario_index'] for r in rows}==set(range(64))
    for name,digest in data['sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    if 'profile_point_path' in data:
        assert hashlib.sha256(Path(data['profile_point_path']).read_bytes()).hexdigest()==data['profile_point_sha256']
    summaries=[];maximum_defect=0.;maximum_relation_defect=0.
    for row in rows:
        assert row['minimum_measured_amplitude']>1e-12
        assert len(row['source_contact_relation_certificates'])==12
        for certificate in row['source_contact_relation_certificates']:
            assert certificate['valid'] and certificate['relation_count']==13
            maximum_relation_defect=max(maximum_relation_defect,certificate['maximum_normalized_relation_defect'])
        for family,result in row['targets'].items():
            if result['null_image_rank']:
                assert result['measured_null_relative_defect']<1e-9
                assert result['target_witness_norm']>0
                maximum_defect=max(maximum_defect,result['measured_null_relative_defect'])
            assert result['stacked_rank']==result['measured_odd_rank']+result['null_image_rank'],(row['model'],family,result)
    for model in sorted({r['model'] for r in rows}):
        group=[r for r in rows if r['model']==model]
        summaries.append(dict(model=model,scenario_count=len(group),targets={family:dict(null_image_rank_range=[min(r['targets'][family]['null_image_rank'] for r in group),max(r['targets'][family]['null_image_rank'] for r in group)],null_sensitive_scenario_count=sum(r['targets'][family]['null_image_rank']>0 for r in group)) for family in ('VV','Vgamma','gammagamma','3P')}))
    output=dict(complete=True,count=len(rows),maximum_measured_null_witness_defect=maximum_defect,
        maximum_source_contact_relation_defect=maximum_relation_defect,model_summaries=summaries,
        scope='finite neutral resonance odd-coefficient audit, fixed nuisances, sampled masses and UV basis; no full physical predictor or band',
        sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest(),
            'scripts/verify_odd_alp_target_nulls.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    target=ROOT/'data/vector_sector/odd_alp_target_null_verification.json'
    target.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');print(target.relative_to(ROOT))

if __name__=='__main__':main()
