"""Extract a small conditional original-eleven fit report; no new fit scan.

Run: .venv/bin/python scripts/summarize_odd_preliminary_fit.py
Numerical channel-image projection uses scipy least squares/SVD. GeV units;
all c_i,d_i are dimensionless with the original operator M_V factors retained.
"""
from pathlib import Path
import hashlib
import json
import unittest
import sys
import numpy as np
import scipy.linalg as la
sys.path.insert(0,str(Path(__file__).resolve().parent))
import calibrate_odd_low_energy as C

ROOT=C.ROOT
NAMES=[f'c{i}' for i in range(1,8)]+[f'd{i}' for i in range(1,5)]
CENTRAL=dict(rho_input='rho_charged',kstar_input='kstar_charged',phi_input='phi_evaluation',branch='near_ideal')
HYPOTHESES=('universal','single_trace_positive','single_trace_negative','double_trace')

def read(name):
    return json.loads((ROOT/'data/vector_sector'/name).read_text())

def project_original(rows, embedding, parent_coefficients):
    """Minimum Euclidean norm representative, at fixed nuisance parameters."""
    matrix=rows@embedding
    target=rows@parent_coefficients
    coefficients=la.lstsq(matrix,target,cond=1e-10)[0]
    null=la.null_space(matrix,rcond=1e-10)
    defect=float(la.norm(matrix@coefficients-target)/max(1.,la.norm(target)))
    return coefficients,null,defect

def main():
    files=['odd_low_energy_calibration.json','odd_low_energy_profiles.json','odd_low_energy_continuation.json','odd_low_energy_profile_reanchoring.json','odd_low_energy_calibration_verification.json','odd_low_energy_jacobians.json','stage3_results.json','alda_arxiv_v1.json','pdg2024.json']
    archive={name:read(name) for name in files}
    profiles=archive['odd_low_energy_profiles.json']['results']
    overrides={}
    for name in ('odd_low_energy_continuation.json','odd_low_energy_profile_reanchoring.json'):
        for item in archive[name]['results']:
            overrides.setdefault(item['pair_index'],[]).append((name,item['best']))
    calibration=archive['odd_low_energy_calibration.json']['results']
    local={x['fit_index']:x for x in archive['odd_low_energy_calibration_verification.json']['local_prediction_uncertainties']}
    ass=archive['stage3_results.json']['assumptions'];alda=archive['alda_arxiv_v1.json']
    tree=C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F)
    fast=C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*np.pi*ass['alpha_em']))
    embedding=C.model_maps()['original_11']
    selections=[('central_'+hyp,dict(CENTRAL,source_hypothesis=hyp)) for hyp in HYPOTHESES]
    for key,value in [('rho_input','rho_neutral'),('kstar_input','kstar_neutral'),('phi_input','phi_average')]:
        scenario=dict(CENTRAL,source_hypothesis='single_trace_negative');scenario[key]=value
        selections.append(('change_'+key,scenario))
    # Use exact archived phi label, whose name is not a new measurement definition.
    phi_labels=sorted(set(x['scenario']['phi_input'] for x in profiles)-{'phi_evaluation'})
    selections[-1][1]['phi_input']=phi_labels[0]
    selections += [(control,dict(CENTRAL,source_hypothesis='single_trace_negative')) for control in ('without_broad_rho','without_constrained_Kstar')]
    output=[]
    observations=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text())['observations']
    for label,scenario in selections:
        control=label if label.startswith('without_') else 'primary'
        pair_index,pair=next((i,x) for i,x in enumerate(profiles) if x['scenario']==scenario and x['control']==control)
        ci,base=next((i,x) for i,x in enumerate(calibration) if x['scenario']==scenario and x['control']==control and x['model']=='original_11')
        candidates=[('odd_low_energy_calibration.json',base),('odd_low_energy_profiles.json',pair['best']),*overrides.get(pair_index,[])]
        source,best=min(candidates,key=lambda x:x[1]['chi2'])
        sector,mk=C.J.action_sector(np.array(best['even_parameters']),scenario['source_hypothesis'],ass)
        channels=fast(sector,mk);labels=base['measurement_labels'];rows=np.array([channels[k]['row'] for k in labels])
        coef,null,defect=project_original(rows,embedding,np.array(best['gauge_representative_coefficients']))
        # Width and interference invariance follows the documented epsilon adapter.
        canonical=-coef
        if defect>1e-8 or la.norm(rows@embedding@null)>1e-8:
            raise ValueError('original coefficient projection or nullspace failed')
        if not np.allclose(rows@embedding@coef,-rows@embedding@canonical,rtol=1e-12,atol=1e-12):
            raise ValueError('epsilon coefficient adapter failed')
        diagnostic=pair['local_diagnostics'] if source=='odd_low_energy_profiles.json' else local.get(ci) if source=='odd_low_energy_calibration.json' else None
        sigmas=diagnostic['local_prediction_sigmas'] if diagnostic else [None]*len(labels)
        entries=[]
        for k,pred,pull,sigma in zip(labels,best['predictions'],best['pulls'],sigmas):
            obs=next(x for x in observations if x['channel']==k)
            y,error=C.J.data_units(obs)
            entries.append(dict(channel=k,observable=obs['observable'],unit='dimensionless' if obs['observable']=='branching_fraction' else 'GeV',measurement=y,measurement_sigma=error,prediction=pred,pull=pull,local_prediction_sigma=sigma))
        even_names=['m_rho_GeV','m_Kstar_GeV','Delta_V_GeV2','kappa_V','F_V_GeV']
        if scenario['source_hypothesis']!='universal':even_names.append('z_V_GeV_inverse')
        if scenario['source_hypothesis']=='double_trace':even_names.append('y_V_GeV_inverse')
        even_values=best['even_parameters']
        if len(even_names)!=len(even_values):raise ValueError('source parameter naming does not match fitted vector')
        output.append(dict(label=label,scenario=scenario,model='original_11',control=control,archive_source=source,archive_pair_index=pair_index,chi2=best['chi2'],odd_chi2=best['odd_chi2'],ancillary_chi2=best['ancillary_chi2'],success=best['success'],global_minimum_proven=False,canonical_representative_coefficients=dict(zip(NAMES,canonical.tolist())),reference_backend_representative_coefficients=dict(zip(NAMES,coef.tolist())),coefficient_units='dimensionless',representative_policy='minimum Euclidean norm in the eleven original coefficients at fitted nuisance values; no coefficient errors inferred',fixed_nuisance_measured_rank=11-null.shape[1],fixed_nuisance_coefficient_null_dimension=null.shape[1],fixed_nuisance_null_basis_columns=null.tolist(),projection_relative_defect=defect,null_annihilation_norm=float(la.norm(rows@embedding@null)),even_parameters=even_values,even_parameters_by_name=dict(zip(even_names,even_values)),shared_total_widths_gev=best['shared_total_widths_gev'],observations=entries,local_prediction_uncertainty_scope=diagnostic.get('local_covariance_scope') if diagnostic else 'not recomputed for an improved archived point'))
    if len(output)!=9 or len({x['label'] for x in output})!=9:
        raise ValueError('incomplete focused case set')
    result=dict(scope='One fixed experimental-input setup, four distinct source hypotheses, three one-at-a-time input sensitivities and two data-removal checks for the negative single-trace source case, retained as a conditional working example because it describes the primary data; archive extraction without new minimization.',selection_policy='Charged rho and Kstar mass summaries; evaluated phi electronic width; near-ideal inverse mass root; primary twelve odd summaries. These choices define a conditional comparison, not a uniquely preferred source hypothesis.',conventions=r'\epsilon^{0123}=+1; metric (+---); canonical c_i=-c_i^R and d_i=-d_i^R; old archives remain reference backend inputs.',uncertainty_policy='Quoted diagonal experimental summaries, fixed central Alda neutral inputs, tree-level real-photon model. Null directions are unmeasured; representatives are not eleven separately measured coefficients. Local prediction errors are sign-branch conditional; no global nonlinear confidence intervals.',results=output,sha256={str((ROOT/'data/vector_sector'/name).relative_to(ROOT)):hashlib.sha256((ROOT/'data/vector_sector'/name).read_bytes()).hexdigest() for name in files})
    result['sha256']['data/odd_parity/odd_low_energy_measurements.json']=hashlib.sha256((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_bytes()).hexdigest()
    result['sha256']['scripts/summarize_odd_preliminary_fit.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    path=ROOT/'data/vector_sector/odd_preliminary_fit.json';path.write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Conditional preliminary fit of the original eleven operators','',result['scope'],'',result['selection_policy'],'','$'+result['conventions']+'$','',result['uncertainty_policy'],'','| Case | Source hypothesis | Total chi2 | Odd rank | Unmeasured coefficient directions |','|---|---|---:|---:|---:|']
    for item in output:lines.append(f"| {item['label']} | {item['scenario']['source_hypothesis']} | {item['chi2']:.8g} | {item['fixed_nuisance_measured_rank']} | {item['fixed_nuisance_coefficient_null_dimension']} |")
    lines+=['','## Canonical coefficient representatives','','These dimensionless values use a numerical minimum-norm choice for the unmeasured directions. Their zeros or values do not constrain those directions.','', '| Case | '+' | '.join(NAMES)+' |','|---|'+'---:|'*11]
    for item in output:lines.append('| '+item['label']+' | '+' | '.join(f'{item["canonical_representative_coefficients"][k]:.7g}' for k in NAMES)+' |')
    lines+=['','Complete predictions, pulls, available conditional local prediction errors, null bases, nuisance values and input hashes are in `odd_preliminary_fit.json`. No new optimizer search is performed. No individual coefficient uncertainty is assigned.','']
    (path.with_suffix('.md')).write_text('\n'.join(lines))
    print(path.relative_to(ROOT));print(path.with_suffix('.md').relative_to(ROOT))

class ProjectionTests(unittest.TestCase):
    def test_representative_preserves_predictions_and_exposes_ambiguity(self):
        rows=np.array([[1.,2.,0.],[0.,0.,1.]])
        embedding=np.eye(3)
        representative,null,defect=project_original(rows,embedding,np.array([3.,4.,5.]))
        self.assertLess(defect,1e-12)
        self.assertEqual(null.shape,(3,1))
        np.testing.assert_allclose(rows@(representative+7*null[:,0]),[11.,5.],atol=1e-12)
        canonical=-representative
        np.testing.assert_allclose(-rows@canonical,[11.,5.],atol=1e-12)

if __name__=='__main__':main()
