"""Certify conditional calibration completeness, inverse mass roots and gauges.
Run .venv/bin/python scripts/verify_odd_low_energy_calibration.py.
"""
import hashlib,json
import numpy as np
import scipy.linalg as la
import sympy as sp
from scipy.stats import chi2 as chi2_distribution
import calibrate_odd_low_energy as C
import calibrate_vector_stage3 as S
ROOT=C.ROOT

def branch_certificate(parameters,hypothesis,branch,assumptions):
    p=np.asarray(parameters);sector,mk=C.J.action_sector(p,hypothesis,assumptions)
    roots=S.extended_branches([sector.masses[0],mk,*sector.masses[1:]])
    if len(roots)!=2:return {'valid':False,'reason':'inverse mass map does not have two real roots'}
    roots.sort(key=lambda x:abs(x['theta_deg']));selected=roots[0 if branch=='near_ideal' else 1]
    defect=float(max(abs(p[2]-selected['delta_gev2']),abs(p[3]-selected['kappa'])))
    return dict(valid=bool(defect<1e-8),inverse_root_absolute_defect=defect,selected_theta_deg=selected['theta_deg'])

def source_certificate(parameters,hypothesis,assumptions):
    p=np.asarray(parameters);sector,mk=C.J.action_sector(p,hypothesis,assumptions)
    if hypothesis=='universal':margin=float(p[4])
    elif hypothesis in ('single_trace_positive','single_trace_negative'):
        fl=p[4]+2*p[5]*assumptions['xl_gev2'];fs=p[4]+2*p[5]*assumptions['xs_gev2'];margin=float(min(fl,fs if hypothesis=='single_trace_positive' else -fs))
    else:
        reference=C.J.module('vector_sector').source_residues(sector.generators,np.diag([assumptions['xl_gev2']]*2+[assumptions['xs_gev2']]),1.)
        margin=float(min(np.asarray(sector.residues)*np.sign(reference)))
    midpoint=4*(p[1]**2-p[0]**2)/(3*(assumptions['xl_gev2']-assumptions['xs_gev2']))
    physical=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text())['physical_masses_gev']
    open_channels=bool(sector.masses[0]<physical['etap'] and sector.masses[1]<physical['etap'] and sector.masses[0]>physical['eta'] and sector.masses[1]>physical['eta'] and sector.masses[2]>physical['eta'] and mk>max(physical['K_neutral'],physical['K_charged']))
    return dict(valid=bool(margin>0 and open_channels),on_shell_channels_open=open_channels,source_domain_margin_gev=margin,root_midpoint_gap=float(abs(p[3]-midpoint)),near_domain_boundary=bool(margin<1e-7 or abs(p[3]-midpoint)<1e-7))

def local_diagnostics(result,data,ass,fast,mapping):
    obs=[x for x in data['observations'] if x.get('status','included')!='deferred'];labels=result['measurement_labels'];training=[x for x in obs if x['channel'] in labels];held=[x for x in obs if x['channel'] not in labels]
    p=np.asarray(result['even_parameters']);coef=np.asarray(result['gauge_representative_coefficients']);tk=sorted(result['shared_total_widths_gev']);tv=np.array([result['shared_total_widths_gev'][k] for k in tk]);ts=np.array([data['total_widths'][k]['sigma_gev'] for k in tk]);npos=len(p)
    # Parameter scales only balance the linear algebra; they are not priors.
    scales=np.array([.001,.001,.01,.01,.01]+([] if result['scenario']['source_hypothesis']=='universal' else [.1])+([.1] if result['scenario']['source_hypothesis']=='double_trace' else []));hyp=result['scenario']['source_hypothesis']
    pdg=json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text());es=np.array([pdg['masses'][x]['sigma_mev']/1000 for x in (result['scenario']['rho_input'],result['scenario']['kstar_input'],'omega','phi')]+[pdg['leptonic_widths'][x]['sigma_kev']*1e-6 for x in ('rho','omega','phi')]);
    if result['scenario']['phi_input']!='phi_evaluation':es[-1]=pdg['alternatives']['phi_ee_average']['sigma_kev']*1e-6
    errors=np.array([C.J.data_units(x)[1] for x in obs]);selected=np.array([obs.index(x) for x in training]);withheld=np.array([obs.index(x) for x in held],int)
    def evaluate(point,totals):
        sec,mk=C.J.action_sector(point,hyp,ass);channels=fast(sec,mk);td=dict(zip(tk,totals));factors=np.array([channels[x['channel']]['phase_factor']/(td[x['total_width_key']] if x['observable']=='branching_fraction' else 1.) for x in obs]);rows=np.array([channels[x['channel']]['row'] for x in obs]);amp=np.array([channels[x['channel']]['offset'] for x in obs])+rows@coef
        return factors*amp**2,2*factors[:,None]*amp[:,None]*(rows@mapping),C.J.even_observables(sec,mk,ass)
    values,odd,even=evaluate(p,tv);nuisance=[];stepdefect=0.
    for i,scale in enumerate(scales):
        derivatives=[]
        for fraction in (1e-4,5e-5):
            h=scale*fraction;hi=p.copy();lo=p.copy();hi[i]+=h;lo[i]-=h;up=evaluate(hi,tv);dn=evaluate(lo,tv);derivatives.append(np.r_[(up[0]-dn[0])*scale/(2*h)/errors,(up[2]-dn[2])*scale/(2*h)/es,np.zeros(len(tk))])
        stepdefect=max(stepdefect,np.linalg.norm(derivatives[1]-derivatives[0])/max(1.,np.linalg.norm(derivatives[1])));nuisance.append(derivatives[1])
    for i,key in enumerate(tk):
        d=np.zeros(len(obs)+7+len(tk))
        for k,x in enumerate(obs):
            if x.get('total_width_key')==key:d[k]=-values[k]*ts[i]/tv[i]/errors[k]
        d[len(obs)+7+i]=1.;nuisance.append(d)
    allJ=np.column_stack([np.vstack([odd/errors[:,None],np.zeros((7+len(tk),mapping.shape[1]))]),np.array(nuisance).T]);idx=np.r_[selected,np.arange(len(obs),len(obs)+7+len(tk))];A=allJ[idx];norm=np.linalg.norm(A,axis=0);keep=norm>1e-12*max(1.,norm.max(initial=0));B=A[:,keep]/norm[keep];u,sv,vh=la.svd(B,full_matrices=False);rank=int(np.sum(sv>1e-10*max(1.,sv[0])));inverse=(vh[:rank].T/sv[:rank])@u[:,:rank].T
    H=allJ[selected][:,keep]/norm[keep];train_defect=float(np.linalg.norm(H-H@inverse@B)/max(1.,np.linalg.norm(H)));assert train_defect<1e-8,'training observable does not annihilate numerical nullspace';cov=(H@inverse)@(H@inverse).T*errors[selected,None]*errors[None,selected]
    held_results=[]
    for i in withheld:
        h=allJ[i];balanced=h[keep]/norm[keep];reconstructed=balanced@inverse@B;defect=np.linalg.norm(balanced-reconstructed)/max(1.,np.linalg.norm(balanced));excluded_response=np.linalg.norm(h[~keep])/max(1.,np.linalg.norm(h));identifiable=bool(max(defect,excluded_response)<1e-8)
        held_results.append(dict(channel=obs[i]['channel'],joint_nuisance_identifiable=identifiable,normalized_training_null_defect=float(max(defect,excluded_response)),prediction=float(values[i]) if identifiable else None,local_sigma=float(np.linalg.norm(balanced@inverse)*errors[i]) if identifiable else None,prediction_status='local conditional sign-gauge prediction' if identifiable else 'unidentifiable; prediction withheld'))
    return dict(training_observable_null_relative_defect=train_defect,derivative_step_relative_change=float(stepdefect),joint_visible_rank=rank,joint_parameter_column_count=A.shape[1],unmeasured_tangent_dimensions=A.shape[1]-rank,joint_balanced_singular_values=sv.tolist(),local_prediction_covariance=cov.tolist(),local_prediction_sigmas=np.sqrt(np.diag(cov)).tolist(),local_covariance_scope='conditional on fitted sign branch, fixed central Alda and diagonal summaries; physical observables only; no coefficient or global nonlinear uncertainty',held_out=held_results)


def main():
    p=ROOT/'data/vector_sector/odd_low_energy_calibration.json';j=json.loads(p.read_text());ass=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
    for name,digest in j['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise AssertionError('stale calibration '+name)
    expected={(i,m,c) for i in range(64) for m in C.MODELS for c in C.CONTROLS};actual={(r['scenario_index'],r['model'],r['control']) for r in j['results']}
    assert j['complete_64_scenario_run'] and expected==actual and len(j['results'])==1792
    certificates=[];source_certificates=[];groups={};failures=[]
    for index,r in enumerate(j['results']):
        certificate=branch_certificate(r['even_parameters'],r['scenario']['source_hypothesis'],r['scenario']['branch'],ass)
        certificates.append(dict(fit_index=index,**certificate))
        source=source_certificate(r['even_parameters'],r['scenario']['source_hypothesis'],ass);source_certificates.append(dict(fit_index=index,**source))
        if not source['valid']:failures.append(dict(fit_index=index,failure='source sign certificate'))
        if not certificate['valid']:failures.append(dict(fit_index=index,failure='mass branch certificate'))
        if not r['action_image_certified']:failures.append(dict(fit_index=index,failure='fixed gauge containment'))
        if r['visible_start_rank']!=r['visible_fit_rank']:failures.append(dict(fit_index=index,failure='visible rank changed'))
        if not r['success']:failures.append(dict(fit_index=index,failure='optimizer did not converge'))
        groups.setdefault((r['control'],r['model']),[]).append(r)
    summary=[]
    for (control,model),fits in groups.items():
        summary.append(dict(control=control,model=model,count=len(fits),success_count=sum(x['success'] for x in fits),image_certified_count=sum(x['action_image_certified'] for x in fits),odd_chi2_range=[min(x['odd_chi2'] for x in fits),max(x['odd_chi2'] for x in fits)],ancillary_chi2_range=[min(x['ancillary_chi2'] for x in fits),max(x['ancillary_chi2'] for x in fits)],maximum_evaluated_gauge_defect=max(x['maximum_evaluated_gauge_containment_defect'] for x in fits),active_search_bound_fit_count=sum(bool(x['active_search_bounds']) for x in fits),optimizer_failure_count=sum(not x['success'] for x in fits),adequacy_eligible_fit_count=sum(x['success'] and x['action_image_certified'] and not x['active_search_bounds'] and x['visible_fit_rank']<len(x['measurement_labels']) for x in fits)))
    nested=[]
    for i in range(64):
        for control in C.CONTROLS:
            fits={r['model']:r for r in j['results'] if r['scenario_index']==i and r['control']==control};base=np.asarray(fits['soft_leading']['predictions'])
            for model,r in fits.items():
                nested.append(dict(scenario_index=i,control=control,model=model,relative_prediction_changes_from_soft_leading=((np.asarray(r['predictions'])-base)/base).tolist(),status='gauge-independent measured predictions; sign ambiguity prevents unique signed correction decomposition'))
    data=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text());alda=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text());tree=C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F);fast=C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*np.pi*ass['alpha_em']));maps=C.model_maps();local=[]
    for i,r in enumerate(j['results']):
        if r['success'] and r['action_image_certified'] and not r['active_search_bounds'] and certificates[i]['valid'] and source_certificates[i]['valid'] and not source_certificates[i]['near_domain_boundary'] and r['visible_fit_rank']<len(r['measurement_labels']):
            diagnostics=local_diagnostics(r,data,ass,fast,maps[r['model']]);df=len(r['measurement_labels'])+14-diagnostics['joint_visible_rank'];pvalue=float(chi2_distribution.sf(r['chi2'],df)) if df>0 else None
            local.append(dict(fit_index=i,local_covariance_valid=bool(diagnostics['derivative_step_relative_change']<=1e-6),local_effective_degrees_of_freedom=df,descriptive_chi2_survival_probability=pvalue,maximum_odd_channel_pull=max(abs(x) for x in r['pulls']),descriptive_residual_flag=bool((pvalue is not None and pvalue<.05) or max(abs(x) for x in r['pulls'])>3),adequacy_statistics_policy='descriptive conditional local-Gaussian diagnostic; thresholds not a preregistered selection rule; saturated p-values undefined',**diagnostics))
    image_comparisons=[]
    for r in j['results']:
        if r['model']!='soft_leading':continue
        sec,mk=C.J.action_sector(np.asarray(r['even_parameters']),r['scenario']['source_hypothesis'],ass);channels=fast(sec,mk);rows=np.array([channels[k]['row'] for k in r['measurement_labels']]);leading_gauge,_=C.image_gauge(rows,maps['soft_leading']);original_gauge,_=C.image_gauge(rows,maps['original_11']);comparison=dict(scenario_index=r['scenario_index'],control=r['control'],original_in_leading_defect=C.containment(rows,maps['original_11'],leading_gauge),leading_in_original_defect=C.containment(rows,maps['soft_leading'],original_gauge),model_images=[])
        for model,mapping in maps.items():
            gauge,rank=C.image_gauge(rows,mapping);comparison['model_images'].append(dict(model=model,rank=rank,leading_in_model_defect=C.containment(rows,maps['soft_leading'],gauge)))
        image_comparisons.append(comparison)
    s1,s2=sp.symbols('s_1 s_2',positive=True);H=sp.MatrixSymbol('H',2,2);V=sp.MatrixSymbol('V',2,2);theta=sp.symbols('theta',real=True);U=sp.Matrix([[sp.cos(theta),-sp.sin(theta)],[sp.sin(theta),sp.cos(theta)]]);D=sp.diag(1/s1,1/s2);noise_covariance=sp.eye(2);unit_covariance=sp.simplify(U.T*noise_covariance*U);local_formula=sp.latex(H*V*D*unit_covariance*D*V.T*H.T)
    out=dict(complete=True,observable_image_comparisons=image_comparisons,local_prediction_uncertainties=local,local_covariance_derivation=local_formula,fit_count=len(j['results']),branch_certified_count=sum(x['valid'] for x in certificates),source_certified_count=sum(x['valid'] for x in source_certificates),source_domain_certificates=source_certificates,failures=failures,branch_certificates=certificates,group_summary=summary,nested_measured_prediction_changes=nested,sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest(),'scripts/verify_odd_low_energy_calibration.py':hashlib.sha256((ROOT/'scripts/verify_odd_low_energy_calibration.py').read_bytes()).hexdigest()})
    target=ROOT/'data/vector_sector/odd_low_energy_calibration_verification.json';target.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(target.relative_to(ROOT));print('fits',len(j['results']),'branch certificates',out['branch_certified_count'],'failures',len(failures))
    transcript=ROOT/'data/vector_sector/odd_low_energy_calibration.md'
    with transcript.open('w') as h:
        h.write('# Conditional low-energy calibration\n\n')
        for label,expr in j['symbolic_derivation'].items():h.write(label+':\n$$\n'+expr+'\n$$\n\n')
        h.write('Original measurements occur once; fixed diagonal summaries, central Alda source inputs and declared tree kinematics. No OPE/OZI/naturalness priors.\n\n')
        for row in summary:h.write(f"{row['control']}, {row['model']}: odd chi2 range {row['odd_chi2_range']}; {row['success_count']}/64 optimizer successes.\n\n")
        h.write('Saturated interpolation is not a truncation-adequacy test. Invisible action directions are unmeasured. Fixed-start exhaustive signs plus finite nonlinear starts do not establish a global nonlinear minimum. No ALP band or model-selection conclusion is produced.\n')
    print(transcript.relative_to(ROOT))

if __name__=='__main__':main()
