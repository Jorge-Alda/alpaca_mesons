"""Targeted finite restart from source-profile multiplier-one minima.

Run: .venv/bin/python scripts/restart_odd_low_energy_profile_references.py
GeV; original likelihood and source hypotheses remain unchanged.
"""
import hashlib
import json
import math
import numpy as np
import profile_odd_low_energy as P
from reanchor_odd_low_energy_profiles import physical_source_fit
from continue_odd_low_energy_source_profiles import source_point

ROOT=P.ROOT


def main():
    paths=['data/vector_sector/odd_low_energy_profiles.json',
           'data/vector_sector/odd_low_energy_calibration.json',
           'data/odd_parity/odd_low_energy_measurements.json',
           'data/vector_sector/pdg2024.json','data/vector_sector/stage3_results.json',
           'data/vector_sector/alda_arxiv_v1.json',
           'data/vector_sector/odd_low_energy_reanchored_profiles.json',
           'scripts/restart_odd_low_energy_profile_references.py',
           'scripts/reanchor_odd_low_energy_profiles.py',
           'scripts/continue_odd_low_energy_source_profiles.py',
           'scripts/profile_odd_low_energy.py','scripts/calibrate_odd_low_energy.py']
    profiles,calibration,data,pdg,stage3,alda,reprofile=[json.loads((ROOT/p).read_text()) for p in paths[:7]]
    ass=stage3['assumptions']
    tree=P.C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=P.C.J.F,fpi=P.C.J.F)
    fast=P.C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*math.pi*ass['alpha_em']))
    output=dict(scope='targeted multiplier-one finite restarts plus refreshed source grids',
                global_nonlinear_minimum_proven=False,
                sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},results=[])
    path=ROOT/'data/vector_sector/odd_low_energy_targeted_profile_restart.json'
    for previous in reprofile['results']:
        if previous['profile_reference_resolved']:
            continue
        old=profiles['results'][previous['pair_index']]
        base=next(r for r in calibration['results'] if r['model']=='soft_leading' and r['scenario_index']==old['scenario_index'] and r['control']==old['control'])
        likelihood=P.ConditionalLikelihood(base,data,pdg,ass,fast)
        incumbent=min((q for q in previous['one_sided_source_profile']['points'] if q['valid']),key=lambda q:q['chi2'])
        start=np.asarray(incumbent['coordinates'])
        candidate=dict(incumbent);candidate.pop('candidate_fits',None);candidate['method']='accepted profile incumbent'
        fits=[candidate]
        for method in ('original','physical'):
            try:
                if method=='physical':
                    summary=physical_source_fit(likelihood,start)
                else:
                    fit=likelihood.fit(start,max_nfev=1200);summary=likelihood.summary(fit)
                    summary.update(method='original standardized coordinates',valid=bool(fit.success and summary['chi2']<1e20))
                fits.append(summary)
            except (ValueError,np.linalg.LinAlgError,OverflowError) as error:
                fits.append(dict(valid=False,method=method,reason=str(error)))
        best=min((q for q in fits if q['valid']),key=lambda q:q['chi2'])
        reference=np.asarray(best['coordinates'])
        descriptor=P.closest_source_quantity(likelihood,reference)
        if descriptor[0]!='signed_F_s':
            raise ValueError('closest source changed; physical coordinate map must be generalized')
        baseline=P.source_quantity(likelihood,reference,descriptor)
        result=dict(pair_index=previous['pair_index'],scenario_index=old['scenario_index'],control=old['control'],
                    old_reference_chi2=previous['reference']['chi2'],profile_incumbent_chi2=incumbent['chi2'],
                    best=best,fits=fits,additional_unconstrained_improvement=incumbent['chi2']-best['chi2'],
                    source_profile_recomputed=True,source_reference_gev=baseline,source_points=[])
        # Re-anchoring changes the reference and source magnitude; refresh the
        # complete same declared seven-point grid, retaining all outcomes.
        for multiplier in (.05,.25,.5,1.,2.,4.,8.):
            target=max(1e-14,baseline*multiplier)
            seeds=[reference,start]
            if result['source_points'] and result['source_points'][-1]['valid']:
                seeds.append(np.asarray(result['source_points'][-1]['coordinates']))
            point=source_point(likelihood,reference,target,seeds,max_nfev=800)
            point.update(multiplier=multiplier,signed_source_target_gev=target,
                         delta_chi2=point['chi2']-best['chi2'] if point['valid'] else None)
            result['source_points'].append(point)
            print(f'{path.relative_to(ROOT)}: pair {previous["pair_index"]}, target {multiplier} complete',flush=True)
        result['profile_reference_resolved']=not any(q['valid'] and q['delta_chi2']< -1e-4 for q in result['source_points'])
        result['unresolved_source_multipliers']=[q['multiplier'] for q in result['source_points'] if not q['valid']]
        result['lower_reference_points']=[q['multiplier'] for q in result['source_points'] if q['valid'] and q['delta_chi2']< -1e-4]
        result['maximum_every_evaluation_image_defect']=likelihood.maximum_image_defect
        output['results'].append(result);output['complete']=False
        path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    output['complete']=True
    path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
