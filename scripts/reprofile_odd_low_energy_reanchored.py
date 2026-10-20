"""Recompute declared profile grids at re-anchored odd-sector references.

Run: .venv/bin/python scripts/reprofile_odd_low_energy_reanchored.py
Units: GeV. Preserve the frozen likelihood and computational-priority rule.
"""
import hashlib
import json
import math
import numpy as np
import profile_odd_low_energy as P
from continue_odd_low_energy_source_profiles import source_point

ROOT=P.ROOT


def main():
    paths=['data/vector_sector/odd_low_energy_profiles.json',
           'data/vector_sector/odd_low_energy_calibration.json',
           'data/odd_parity/odd_low_energy_measurements.json',
           'data/vector_sector/pdg2024.json','data/vector_sector/stage3_results.json',
           'data/vector_sector/alda_2026nrd_arxiv_v1.json',
           'data/vector_sector/odd_low_energy_profile_reanchoring.json',
           'scripts/reprofile_odd_low_energy_reanchored.py',
           'scripts/continue_odd_low_energy_source_profiles.py',
           'scripts/profile_odd_low_energy.py','scripts/calibrate_odd_low_energy.py']
    profiles,calibration,data,pdg,stage3,alda,reanchoring=[json.loads((ROOT/p).read_text()) for p in paths[:7]]
    ass=stage3['assumptions']
    tree=P.C.J.module('alda_2026nrd_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree',F=P.C.J.F,fpi=P.C.J.F)
    fast=P.C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*math.pi*ass['alpha_em']))
    output=dict(scope='declared one-sided profile grid at re-anchored references',
                physical_profile_priority_policy='interior converged p>=.05 and maxoddpull<=3',
                source_grid=[.05,.25,.5,1.,2.,4.,8.], global_nonlinear_minimum_proven=False,
                sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},results=[])
    path=ROOT/'data/vector_sector/odd_low_energy_reanchored_profiles.json'
    for anchored in reanchoring['results']:
        old=profiles['results'][anchored['pair_index']]
        base=next(r for r in calibration['results'] if r['model']=='soft_leading' and r['scenario_index']==old['scenario_index'] and r['control']==old['control'])
        likelihood=P.ConditionalLikelihood(base,data,pdg,ass,fast)
        reference=likelihood.rebase(anchored['best'])
        descriptor=P.closest_source_quantity(likelihood,reference)
        if descriptor[0]!='signed_F_s':
            raise ValueError('new closest source requires another explicit physical coordinate map')
        baseline=P.source_quantity(likelihood,reference,descriptor)
        maxpull=max(abs(x) for x in anchored['best']['pulls'])
        # Every affected case fails this necessary predeclared criterion.
        # Thus covariance eigendirections are not allocated, regardless of p.
        if maxpull<=3:
            raise ValueError('physical profile eligibility needs full local diagnostics')
        result=dict(pair_index=anchored['pair_index'],scenario_index=old['scenario_index'],control=old['control'],
                    reference=anchored['best'],profiles=[],profile_priority=False,
                    physical_profile_omission_reason='maximum odd-channel absolute pull exceeds predeclared priority threshold3',
                    maximum_absolute_odd_pull=maxpull,
                    one_sided_source_profile=dict(quantity=descriptor[0],reference_signed_source_gev=baseline,
                         positive_numerical_floor_gev=1e-14,one_sided=True,symmetric_gaussian_interval=False,points=[]))
        points=result['one_sided_source_profile']['points']
        starts=[reference]
        competing=[q for q in anchored['fits'] if q['valid'] and q.get('coordinates') and q.get('amplitude_signs')!=anchored['best']['amplitude_signs']]
        if competing:
            starts.append(np.asarray(min(competing,key=lambda q:q['chi2'])['coordinates']))
        for multiplier in output['source_grid']:
            target=max(1e-14,baseline*multiplier)
            seeds=starts.copy()
            if points and points[-1]['valid']:
                seeds.append(np.asarray(points[-1]['coordinates']))
            fitted=source_point(likelihood,reference,target,seeds,max_nfev=800)
            fitted.update(multiplier=multiplier,signed_source_target_gev=target,
                          delta_chi2=fitted['chi2']-anchored['best']['chi2'] if fitted['valid'] else None)
            points.append(fitted)
            print(f'{path.relative_to(ROOT)}: pair {anchored["pair_index"]}, target {multiplier} complete',flush=True)
        result['profile_reference_resolved']=not any(q['valid'] and q['delta_chi2']< -1e-4 for q in points)
        result['lower_reference_points']=[q['multiplier'] for q in points if q['valid'] and q['delta_chi2']< -1e-4]
        result['unresolved_source_multipliers']=[q['multiplier'] for q in points if not q['valid']]
        result['maximum_every_evaluation_image_defect']=likelihood.maximum_image_defect
        output['results'].append(result);output['complete']=False
        path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    output['complete']=True
    path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
