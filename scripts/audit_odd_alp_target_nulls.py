"""Finite conditional measured-odd-null audit; no widths/bands/generic ranks.

Run .venv/bin/python scripts/audit_odd_alp_target_nulls.py only after primary
profile points are ready. Each action model remains independent.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import scipy.linalg as la
import sympy as sp
import calibrate_odd_low_energy as C

ROOT=C.ROOT
ALP=C.J.module('odd_alp_target_null')
TOL=1e-10
MASS_SQUARED=(.01,.16,1.44)
UV_BASIS=np.eye(4)


def null_diagnostic(measured, target, mapping):
    """Finite SVD null witness in declared model coefficient coordinates."""
    A=np.asarray(measured)@mapping; B=np.asarray(target)@mapping
    scales=np.linalg.norm(A,axis=0)
    scales[scales<TOL*max(1.,scales.max(initial=0))]=1.
    balanced=A/scales
    _,s,vh=la.svd(balanced,full_matrices=True)
    rank=int(sum(s>TOL*max(1.,s[0] if len(s) else 0.)))
    N=vh[rank:].T/scales[:,None]
    residual=B@N
    target_rank=C.J.rank(B)[0]
    joint_rank=C.J.rank(np.vstack((A,B)))[0]
    # SVD singular values of the null image use the same coefficient scaling.
    singular=la.svdvals(residual) if residual.size else np.array([])
    reference=max(1.,la.norm(B/scales,2)) if B.size else 1.
    null_rank=int(sum(singular>TOL*reference))
    result=dict(measured_odd_rank=rank,target_rank=target_rank,
        stacked_rank=joint_rank,null_image_rank=null_rank,
        maximum_null_image=float(singular[0] if len(singular) else 0.),
        null_image_reference=float(reference),coefficient_null_dimension=N.shape[1])
    if null_rank:
        _,_,rv=la.svd(residual,full_matrices=False)
        witness=N@rv[0];witness/=la.norm(witness)
        image=B@witness
        defect=la.norm(A@witness)/max(1.,la.norm(A)*la.norm(witness))
        if defect>1e-9:raise ValueError('failed measured-null witness')
        result.update(model_coefficient_witness=witness.tolist(),
            parent_coefficient_witness=(mapping@witness).tolist(),
            measured_null_relative_defect=float(defect),
            target_witness_norm=float(la.norm(image)),
            largest_target_witness_row=int(np.argmax(np.abs(image))))
    return result


def symbolic_method():
    a,b=sp.symbols('a b', nonzero=True)
    A=sp.Matrix([[a,b]]);N=sp.Matrix.hstack(*A.nullspace())
    x,y=sp.symbols('x y');B=sp.Matrix([[x,y]])
    return dict(measured=sp.latex(A),null_basis=sp.latex(N),
                target_null_image=sp.latex(B*N),null_certificate=sp.latex(A*N))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--limit-scenarios',type=int,default=0)
    parser.add_argument('--profile-points',type=Path)
    args=parser.parse_args()
    fitpath=ROOT/'data/vector_sector/odd_low_energy_calibration.json'
    calibration=json.loads(fitpath.read_text())
    if not calibration['complete_64_scenario_run']:raise ValueError('incomplete calibration')
    for name,digest in calibration['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('stale calibration input '+name)
    fits=[r for r in calibration['results'] if r['control']=='primary']
    # Optional explicit replacement protocol: list of {scenario_index,model,
    # even_parameters}; producer identifies stronger sign-search primary points.
    replacements={}
    if args.profile_points:
        replacements={(r['scenario_index'],r['model']):r
                      for r in json.loads(args.profile_points.read_text())}
    jac=json.loads((ROOT/'data/vector_sector/odd_low_energy_jacobians.json').read_text())
    data=json.loads((ROOT/jac['measurement_snapshot']).read_text())
    obs={r['channel']:r for r in data['observations'] if r.get('status','included')!='deferred'}
    assumptions=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
    alda=json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())
    tree=C.J.module('alda_2026nrd_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree')
    chi=tree['pi0']['chi'];e=np.sqrt(4*np.pi*assumptions['alpha_em'])
    fast=C.FastChannels(tree,alda['physical_masses_gev'],e)
    models=C.model_maps();results=[]
    for fit in fits:
        si=fit['scenario_index'];model=fit['model']
        if args.limit_scenarios and si>=args.limit_scenarios:continue
        replacement=replacements.get((si,model),{})
        point=replacement.get('even_parameters',fit['even_parameters'])
        sector,mk=C.J.action_sector(np.array(point),fit['scenario']['source_hypothesis'],assumptions)
        channels=fast(sector,mk)
        # Nonzero measured amplitudes make the amplitude/width row nullspaces
        # identical at fixed nuisance values. Report that domain explicitly.
        coefficients=np.array(replacement.get('gauge_representative_coefficients',fit['gauge_representative_coefficients']))
        amplitudes=np.array([channels[k]['offset']+channels[k]['row']@coefficients
                             for k in fit['measurement_labels']])
        if np.any(np.abs(amplitudes)<1e-12):raise ValueError('zero measured amplitude: width-null conversion invalid')
        measured=np.array([channels[k]['row'] for k in fit['measurement_labels']])
        target={k:[] for k in ('VV','Vgamma','gammagamma','3P')};labels={k:[] for k in target};certificates=[]
        for w in MASS_SQUARED:
            for u,uv in enumerate(UV_BASIS):
                U,Y,_=ALP.sources(chi,w,uv[:3],uv[3],alda['m0_squared_gev2'])
                target_rows=ALP.target_rows(chi,U,Y,w,sector,e)
                certificate=ALP.relation_certificate(target_rows);certificates.append(certificate)
                if not certificate['valid']:raise ValueError('source/contact action-relation gate failed: '+str(certificate))
                for label,row in target_rows.items():
                    family=label.split(':')[0];target[family].append(row)
                    labels[family].append(dict(mass_squared_gev2=w,uv_basis_index=u,channel=label))
        diagnostics={family:null_diagnostic(measured,np.array(rows),models[model])
                     for family,rows in target.items()}
        results.append(dict(scenario_index=si,model=model,scenario=fit['scenario'],
            point_source='profile_primary_replacement' if (si,model) in replacements else 'calibration_primary',
            even_parameters=list(point),minimum_measured_amplitude=float(min(abs(amplitudes))),
            target_labels=labels,targets=diagnostics,source_contact_relation_certificates=certificates))
    inputs=['scripts/audit_odd_alp_target_nulls.py','scripts/implementation/odd_alp_target_null.py',
            'data/vector_sector/odd_low_energy_calibration.json','data/vector_sector/odd_low_energy_jacobians.json',
            'data/vector_sector/odd_cubic_parent_kernels.json','data/vector_sector/alda_2026nrd_arxiv_v1.json',
            'data/vector_sector/stage3_results.json','data/vector_sector/odd_leading_action_count.json',
            'data/vector_sector/odd_three_flavor_action_quotient.json',jac['measurement_snapshot']]
    out=dict(scope='finite conditional neutral bounded resonance odd-coefficient null audit',
        complete_64_scenario_run=not args.limit_scenarios,expected_count=64*7,count=len(results),
        mass_squared_probes_gev2=list(MASS_SQUARED),uv_coordinates=['cAP_u','cAP_d','cAP_s','cGA'],
        fa_gev=1.,rank_tolerance=TOL,symbolic_method=symbolic_method(),
        limits=['fixed fitted even parameters; nuisance nulls not audited',
                'three sampled masses do not certify generic or all-mass ranks',
                'VV and Vgamma independent hadronic anomaly/local coefficients uncomputed',
                'photon contact under inherited normalization has zero odd-coefficient derivative',
                'on-shell narrow daughter projections; no widths or thresholds or bands',
                '3P zero applies to single odd insertion at retained vacuum tree accuracy'],
        results=results,sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    if args.profile_points:
        out['profile_point_path']=str(args.profile_points)
        out['profile_point_sha256']=hashlib.sha256(args.profile_points.read_bytes()).hexdigest()
    path=ROOT/'data/vector_sector'/('odd_alp_target_nulls_debug.json' if args.limit_scenarios else 'odd_alp_target_nulls.json')
    path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(path.relative_to(ROOT))


if __name__=='__main__':main()
