"""Merge and independently certify finite profiles. Units: GeV.

Run .venv/bin/python scripts/verify_odd_low_energy_profiles.py --merge-chunks.
No global-minimum or statistical coverage claim follows from finite profiles.
"""
import argparse
import hashlib
import json
import math
import numpy as np
import calibrate_odd_low_energy as C
import profile_odd_low_energy as P
import verify_odd_low_energy_calibration as V
ROOT=C.ROOT

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--merge-chunks',action='store_true')
    args=parser.parse_args()
    path=ROOT/'data/vector_sector/odd_low_energy_profiles.json'
    if args.merge_chunks:
        chunks=[json.loads((ROOT/f'data/vector_sector/odd_low_energy_profiles_{start}_64.json').read_text()) for start in (0,64,128,192)]
        assert all(x['sha256']==chunks[0]['sha256'] for x in chunks)
        output={**chunks[0],'results':[r for x in chunks for r in x['results']]}
        output.update(complete_256_pair_run=True,pair_count=len(output['results']),covered_unsaturated_model_cases=2*len(output['results']))
        path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    output=json.loads(path.read_text())
    for name,digest in output['sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest, name
    rows=output['results'];expected={(i,c) for i in range(64) for c in C.CONTROLS}
    assert output['complete_256_pair_run'] and len(rows)==256
    assert {(r['scenario_index'],r['control']) for r in rows}==expected
    saved=json.loads((ROOT/'data/vector_sector/odd_low_energy_calibration.json').read_text())
    old={(r['scenario_index'],r['control']):r for r in saved['results'] if r['model']=='soft_leading'}
    data=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text())
    pdg=json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text())
    ass=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
    alda=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text())
    tree=C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F)
    fast=C.FastChannels(tree,alda['physical_masses_gev'],math.sqrt(4*math.pi*ass['alpha_em']))
    certificates=[];failed_points=0;valid_points=0;valid_source_points=0;failed_source_points=0;max_constraint=0.;max_chi=0.
    for r in rows:
        likelihood=P.ConditionalLikelihood(old[(r['scenario_index'],r['control'])],data,pdg,ass,fast)
        best=r['best'];branch=V.branch_certificate(best['even_parameters'],r['scenario']['source_hypothesis'],r['scenario']['branch'],ass)
        source=V.source_certificate(best['even_parameters'],r['scenario']['source_hypothesis'],ass)
        assert branch['valid'] and source['valid']
        assert r['maximum_every_evaluation_image_defect']<1e-8
        v=np.asarray(best['coordinates']);res=likelihood.residual(v)
        discrepancy=abs(float(res@res)-best['chi2'])/max(1.,best['chi2'])
        assert discrepancy<1e-8
        if 'one_sided_source_profile' in r:
            descriptor=P.closest_source_quantity(likelihood,v)
            for point in r['one_sided_source_profile']['points']:
                if not point['valid']:
                    failed_source_points+=1;continue
                pv=np.asarray(point['coordinates']);res=likelihood.residual(pv)
                assert abs(float(res@res)-point['chi2'])/max(1.,point['chi2'])<1e-8
                source_value=P.source_quantity(likelihood,pv,descriptor)
                assert abs(source_value-point['signed_source_target_gev'])<1e-9*max(r['one_sided_source_profile']['reference_signed_source_gev'],1e-8)
                likelihood.blocks(pv) # Independent domain and image checks.
                valid_source_points+=1
        for profile in r['profiles']:
            assert profile['null_safety_certificate']<1e-8
            c=np.asarray(profile['w'])/np.asarray(profile['tau']);Gref=np.asarray(profile['reference_physical_amplitudes'])
            for point in profile['points']:
                if not point['valid']:
                    failed_points+=1;continue
                v=np.asarray(point['coordinates']);pred,amp=likelihood.predictions(v);res=likelihood.residual(v)
                constraint=abs(float(c@(amp-Gref))-point['physical_contrast_target'])/max(1.,abs(point['physical_contrast_target']))
                discrepancy=abs(float(res@res)-point['chi2'])/max(1.,point['chi2'])
                assert constraint<1e-7 and discrepancy<1e-8
                p=likelihood.point+likelihood.scales*v[:likelihood.np]
                assert V.branch_certificate(p,r['scenario']['source_hypothesis'],r['scenario']['branch'],ass)['valid']
                assert V.source_certificate(p,r['scenario']['source_hypothesis'],ass)['valid']
                valid_points+=1;max_constraint=max(max_constraint,constraint);max_chi=max(max_chi,discrepancy)
        certificates.append(dict(scenario_index=r['scenario_index'],control=r['control'],best_optimizer_success=best['success'],branch_certified=branch['valid'],source_certified=source['valid'],profiles=len(r['profiles']),profile_reference_resolved=r['profile_reference_resolved']))
    summary=dict(complete=True,pair_count=len(rows),model_case_count=512,
        search_fit_count=sum(r['fit_count'] for r in rows),
        improved_pair_count=sum(r['best']['chi2']<r['old_achieved_chi2']-1e-4 for r in rows),
        continuing_second_batch_improvement_count=sum(r['second_batch_continued_improvement'] for r in rows),
        best_optimizer_failure_count=sum(not r['best']['success'] for r in rows),
        profiled_pair_count=sum(bool(r['profiles']) for r in rows),
        physical_profile_count=sum(len(r['profiles']) for r in rows),
        unresolved_profile_reference_count=sum(not r['profile_reference_resolved'] for r in rows),
        one_sided_source_profile_count=sum('one_sided_source_profile' in r for r in rows),
        valid_source_profile_points=valid_source_points,unresolved_source_profile_points=failed_source_points,
        valid_physical_profile_points=valid_points,unresolved_physical_profile_points=failed_points,
        maximum_constraint_relative_defect=max_constraint,maximum_residual_reproduction_relative_defect=max_chi,
        open_95_level_sides=sum(len(p['open_interval_sides']) for r in rows for p in r['profiles']),
        primary_best_chi2=min(r['best']['chi2'] for r in rows if r['control']=='primary'),
        global_nonlinear_minimum_proven=False,certificates=certificates,
        sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest(),
                'scripts/verify_odd_low_energy_profiles.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    target=ROOT/'data/vector_sector/odd_low_energy_profiles_verification.json'
    target.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(target.relative_to(ROOT))

from pathlib import Path
if __name__=='__main__':main()
