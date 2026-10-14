"""Conditional diagonal-summary calibration; GeV; no OPE/coefficient priors.
Run: .venv/bin/python scripts/calibrate_odd_low_energy.py
A pivot coefficient gauge spans each observed image; omitted coefficients
remain unmeasured. Signs enumerated at the start do not prove a global minimum.
"""
from pathlib import Path
import argparse, hashlib, itertools, json, math, sys, time
import numpy as np
import scipy.linalg as la
from scipy.optimize import least_squares
import sympy as sp
import derive_odd_low_energy_jacobians as J

ROOT=J.ROOT
API=J.module('odd_bounded_amplitudes')
MODELS=('original_11','soft_leading','leading_double_trace','leading_single_trace_corrections','soft_first','promoted_leading','full_bounded')
CONTROLS={'primary':set(),'without_broad_rho':J.BROAD_RHO_CHANNELS,'without_constrained_Kstar':J.CONSTRAINED_KSTAR_CHANNELS,'narrow_nonconstrained_summary_control':J.BROAD_RHO_CHANNELS|J.CONSTRAINED_KSTAR_CHANNELS}

def model_maps():
    n=len(API.parent_ids()); eye=np.eye(n); cat=API.parent_catalog()
    leading=API.retained_columns('all_derivatives_soft',0)
    return {'original_11':API.original_embedding(), 'soft_leading':eye[:,leading],
      'leading_double_trace':eye[:,[i for i,x in enumerate(cat) if x['level']==0]],
      'leading_single_trace_corrections':eye[:,[i for i,x in enumerate(cat) if x['traces']==1 and x['level']<=1]],
      'soft_first':eye[:,API.retained_columns('all_derivatives_soft',1)],
      'promoted_leading':eye[:,API.retained_columns('hard_vector_derivatives_promoted',0)],'full_bounded':eye}

class FastChannels:
    """Exact multilinear trace contractions; no fitted polynomial approximation."""
    def __init__(self,tree,physical,e):
        self.tree,self.physical,self.e=tree,physical,e
        self.basis=J.module('vector_sector').physical_generators(0.)
        self.features={}
        self.labels=[('omega_pi0_gamma',1,'pi0','omega','V_to_P'),('rho0_pi0_gamma',0,'pi0','rho','V_to_P'),('rho0_eta_gamma',0,'eta','rho','V_to_P'),('omega_eta_gamma',1,'eta','omega','V_to_P'),('etap_rho0_gamma',0,'etap','rho','P_to_V'),('etap_omega_gamma',1,'etap','omega','P_to_V'),('phi_eta_gamma',2,'eta','phi','V_to_P')]
        for label,i,p,v,o in self.labels:self._radiative(label,tree[p]['TP'],tree[p]['mp'],None)
        chi=tree['pi0']['chi']
        for label,i,j,p in [('Kstar0_K0_gamma',1,2,'K_neutral'),('Kstarplus_Kplus_gamma',0,2,'K_charged'),('rhoplus_piplus_gamma',0,1,'pi_charged')]:
            self._radiative(label,J.E(i,j),np.sqrt((chi[i,i]+chi[j,j])/2),J.E(j,i))
        for p in ('pi0','eta','etap'):self._tff(p)
    def _factors(self,TP,mp,A,B,item):
        U=-math.sqrt(2)*TP/J.F;chi=self.tree['pi0']['chi'];Y=math.sqrt(2)*(TP@chi+chi@TP)/J.F
        value=item['phase']*item['bose_factor']*API._flavor(item,dict(A=A,B=B,U=Y if item['mode'] in ('Y','Ydd') else U,Y=Y,X=2*chi))*J.OPERATOR_MV**(-item['operator_MV_power'])
        p=item['polynomial']; return value*np.array([p[0]+p[3]*mp**2/2,p[1]-p[3]/2])
    def _radiative(self,label,TP,mp,TV):
        bases=self.basis if TV is None else (TV,)
        direct=np.zeros((len(bases),2,len(API.parent_ids())),complex);exchange=np.zeros((len(bases),3,2,len(API.parent_ids())),complex)
        for i,A in enumerate(bases):
            for k,item in enumerate(API.parent_catalog()):
                if item['mode'] in ('fminus','fmdd'):
                    p=item['polynomial']; fac=item['phase']*math.sqrt(2)*self.e/J.F*J.OPERATOR_MV**(-item['operator_MV_power'])*API._flavor(item,dict(A=A,B=TP@API.Q-API.Q@TP,X=2*self.tree['pi0']['chi']))
                    direct[i,:,k]=fac*np.array([p[0]+p[3]*mp**2/2,p[1]-p[3]/2])
                elif item['kind']=='VJP':direct[i,:,k]=-2*self.e*self._factors(TP,mp,A,API.Q,item)
                else:
                    for j,B in enumerate(self.basis):exchange[i,j,:,k]=math.sqrt(2)*self.e*self._factors(TP,mp,A,B,item)
        self.features[label]=(direct.real,exchange.real)
    def _tff(self,p):
        TP=self.tree[p]['TP'];mp=self.tree[p]['mp'];n=len(API.parent_ids());linear=np.zeros((3,n),complex);quad=np.zeros((3,3,n),complex)
        for k,item in enumerate(API.parent_catalog()):
            if item['mode'] in ('fminus','fmdd'):continue
            for i,A in enumerate(self.basis):
                if item['kind']=='VJP':linear[i,k]=-4*math.sqrt(2)*self._factors(TP,mp,A,API.Q,item)[0]
                else:
                    for j,B in enumerate(self.basis):quad[i,j,k]=2*self._factors(TP,mp,A,B,item)[0]
        offset=float((math.sqrt(2)*3*np.trace(2*API.Q@API.Q@TP)/(8*math.pi**2*J.F)).real)
        self.features[p+'_gamma_gamma']=(linear.real,quad.real,offset)
    def __call__(self,sector,mk):
        # Orthogonal ideal generators have norm 1.
        V=np.array([[np.trace(A@B).real for B in self.basis] for A in sector.generators])
        W=np.array([np.trace(sector.response()@B).real for B in self.basis])
        masses=dict(zip(('rho','omega','phi'),sector.masses));result={}
        for label,i,p,v,o in self.labels:
            m=masses[v];D,E=self.features[label];row=(np.einsum('i,ik->k',V[i],D[:,0])+m*m*np.einsum('i,ik->k',V[i],D[:,1])+np.einsum('i,j,ijk->k',V[i],W,E[:,:,0])+m*m*np.einsum('i,j,ijk->k',V[i],W,E[:,:,1]))/m
            pm=self.physical[p];k=(m*m-pm*pm)**3/(96*math.pi*m**3) if o=='V_to_P' else (pm*pm-m*m)**3/(32*math.pi*pm**3)
            result[label]=dict(row=row,offset=0.,phase_factor=k)
        for label,p,m in [('Kstar0_K0_gamma','K_neutral',mk),('Kstarplus_Kplus_gamma','K_charged',mk),('rhoplus_piplus_gamma','pi_charged',masses['rho'])]:
            D,E=self.features[label];row=(D[0,0]+m*m*D[0,1]+np.einsum('j,jk->k',W,E[0,:,0])+m*m*np.einsum('j,jk->k',W,E[0,:,1]))/m
            result[label]=dict(row=row,offset=0.,phase_factor=(m*m-self.physical[p]**2)**3/(96*math.pi*m**3))
        alpha=self.e**2/(4*math.pi)
        for p in ('pi0','eta','etap'):
            L,Q,offset=self.features[p+'_gamma_gamma'];result[p+'_gamma_gamma']=dict(row=W@L+np.einsum('i,j,ijk->k',W,W,Q),offset=offset,phase_factor=math.pi*alpha**2*self.physical[p]**3/4)
        return result

def image_gauge(rows,mapping):
    A=rows@mapping;norm=np.linalg.norm(A,axis=0);keep=np.flatnonzero(norm>1e-10*max(1.,norm.max(initial=0)))
    if not len(keep):return mapping[:,:0],0
    _,r,p=la.qr(A[:,keep]/norm[keep],mode='economic',pivoting=True)
    rank=int(np.sum(np.abs(np.diag(r))>1e-10));chosen=keep[p[:rank]]
    return mapping[:,chosen]/norm[chosen],rank

def containment(rows,full,gauge):
    A=rows@full;B=rows@gauge
    return float(np.linalg.norm(A-B@la.lstsq(B,A,cond=1e-10)[0])/max(1.,np.linalg.norm(A)))

def signed_starts(A,offset,target,errors,count=4):
    signs=np.array(list(itertools.product((-1.,1.),repeat=len(target))))
    coefficients=la.lstsq(A,(signs*target-offset).T,cond=1e-10)[0].T
    prediction=(coefficients@A.T+offset)**2
    score=np.sum(((prediction-target**2)/errors)**2,axis=1)
    order=np.argsort(score,kind='stable');selected=[]
    for i in order:
        if not any(np.linalg.norm(coefficients[i]-x)<1e-7 for x in selected):selected.append(coefficients[i])
        if len(selected)>=count:break
    return selected,float(score[order[0]]),len(signs)

def symbolic_derivation():
    k,T,g,a,b,sigma,y=sp.symbols('k T g a b sigma y',positive=True)
    G=b+a*g; prediction=k*G**2/T;r=(prediction-y)/sigma
    mass=J.derive_mass();poly=mass['polynomial'];kap=mass['kap'];midpoint=sp.solve(sp.diff(poly,kap),kap)[0]
    return {'inverse_mass_root_midpoint':sp.latex(midpoint),'amplitude':sp.latex(G),'observable':sp.latex(prediction),'residual':sp.latex(r),'odd_derivative':sp.latex(sp.diff(r,g)),'total_derivative':sp.latex(sp.diff(r,T))}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--limit-scenarios',type=int,default=0);parser.add_argument('--sign-starts',type=int,default=4);args=parser.parse_args()
    jacpath=ROOT/'data/vector_sector/odd_low_energy_jacobians.json';artifact=json.loads(jacpath.read_text())
    if not artifact['complete_64_scenario_run'] or len(artifact['scenarios'])!=64:raise ValueError('incomplete measurement Jacobians')
    for p,h in artifact['sha256'].items():
        if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise ValueError('stale Jacobian input '+p)
    data=json.loads((ROOT/artifact['measurement_snapshot']).read_text());cal=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text());pdg=json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text());alda=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text());ass=cal['assumptions']
    obs=[x for x in data['observations'] if x.get('status','included')!='deferred'];tree=J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=J.F,fpi=J.F)
    fast=FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*math.pi*ass['alpha_em']));maps=model_maps();results=[];maxdefect=0.
    for si,s in enumerate(artifact['scenarios'][:args.limit_scenarios or 64]):
        mc=next(x for x in cal['scenarios'] if x['rho_input']==s['rho_input'] and x['kstar_input']==s['kstar_input']);wc=next(x for x in mc['width_tests'] if x['phi_input']==s['phi_input']);br=mc['extended_mass_branches'];bi=sorted(range(len(br)),key=lambda i:abs(br[i]['theta_deg']))[0 if s['branch']=='near_ideal' else 1];bc=br[bi]
        src=next(x for x in wc['fixed_angle_cases'] if x['case']==f'extended_mass_branch_{bi}')[J.WIDTH_KEYS[J.HYPOTHESES.index(s['source_hypothesis'])]]
        point=np.array([pdg['masses'][s['rho_input']]['value_mev']/1000,pdg['masses'][s['kstar_input']]['value_mev']/1000,bc['delta_gev2'],bc['kappa'],src['FV_gev']]+([] if s['source_hypothesis']=='universal' else [src['z_gev_inverse']])+([src['y_gev_inverse']] if s['source_hypothesis']=='double_trace' else []))
        scales=np.array(s['even_parameter_column_scales']);sector,mk=J.action_sector(point,s['source_hypothesis'],ass)
        original=J.channel_rows(sector,mk,tree,alda['physical_masses_gev'],fast.e);new=fast(sector,mk)
        for key in new:
            np.testing.assert_allclose(new[key]['row'],original[key]['row'],rtol=3e-12,atol=3e-12);np.testing.assert_allclose(new[key]['offset'],original[key]['offset'],rtol=3e-12,atol=3e-12)
        even_target=np.array([pdg['masses'][x]['value_mev']/1000 for x in (s['rho_input'],s['kstar_input'],'omega','phi')]+[pdg['leptonic_widths'][x]['value_kev']*1e-6 for x in ('rho','omega','phi')]);even_sigma=np.array([pdg['masses'][x]['sigma_mev']/1000 for x in (s['rho_input'],s['kstar_input'],'omega','phi')]+[pdg['leptonic_widths'][x]['sigma_kev']*1e-6 for x in ('rho','omega','phi')])
        if s['phi_input']!='phi_evaluation':
            alt=pdg['alternatives']['phi_ee_average'];even_target[-1]=alt['value_kev']*1e-6;even_sigma[-1]=alt['sigma_kev']*1e-6
        def getsector(z):
            p=point+scales*z
            sec,kk=J.action_sector(p,s['source_hypothesis'],ass)
            # Quadratic inverse-mass roots are separated by their moving midpoint.
            midpoint=4*(p[1]**2-p[0]**2)/(3*(ass['xl_gev2']-ass['xs_gev2']))
            initial_midpoint=4*(point[1]**2-point[0]**2)/(3*(ass['xl_gev2']-ass['xs_gev2']))
            if (p[3]-midpoint)*(point[3]-initial_midpoint)<=0:raise ValueError('mass root branch domain crossed')
            if s['source_hypothesis']=='universal' and p[4]<=0:raise ValueError('universal source sign crossed')
            if s['source_hypothesis'] in ('single_trace_positive','single_trace_negative'):
                fl=p[4]+2*p[5]*ass['xl_gev2'];fs=p[4]+2*p[5]*ass['xs_gev2']
                if fl<=0 or ((fs>0)!=(s['source_hypothesis']=='single_trace_positive')):raise ValueError('single-trace source sign crossed')
            if s['source_hypothesis']=='double_trace':
                reference=J.module('vector_sector').source_residues(sec.generators,np.diag([ass['xl_gev2']]*2+[ass['xs_gev2']]),1.)
                if np.any(np.asarray(sec.residues)*np.asarray(reference)<=0):raise ValueError('double-trace source sign crossed')
            return sec,kk
        def er(z):
            try:sec,kk=getsector(z);return (J.even_observables(sec,kk,ass)-even_target)/even_sigma
            except ValueError:return np.full(7,1e12)
        # Local branch search bounds are explicit numerical domains, never data priors.
        lower=np.full(len(point),-np.inf);upper=np.full(len(point),np.inf)
        lower[:2]=(np.array([.65,.75])-point[:2])/scales[:2];upper[:2]=(np.array([.9,1.05])-point[:2])/scales[:2]
        lower[2:4]=(-np.array([1.,2.])-point[2:4])/scales[2:4];upper[2:4]=(np.array([1.,2.])-point[2:4])/scales[2:4]
        ef=least_squares(er,np.zeros(len(point)),bounds=(lower,upper),xtol=1e-11,ftol=1e-11,gtol=1e-9,max_nfev=400)
        for control,excluded in CONTROLS.items():
            training=[x for x in obs if x['channel'] not in excluded];held=[x for x in obs if x['channel'] in excluded];labels=[x['channel'] for x in training];y,errors=np.array([J.data_units(x) for x in training]).T
            tk=sorted(set(x['total_width_key'] for x in obs if x['observable']=='branching_fraction'));tv=np.array([data['total_widths'][k]['value_gev'] for k in tk]);ts=np.array([data['total_widths'][k]['sigma_gev'] for k in tk]);npos=len(point)
            cache={}
            def blocks(z):
                key=np.asarray(z[:npos]).tobytes()
                if key not in cache:
                    sec,kk=getsector(z[:npos]);chs=fast(sec,kk);cache.clear();cache[key]=(sec,kk,chs)
                sec,kk,chs=cache[key];tot=tv+ts*z[npos:npos+len(tk)];td=dict(zip(tk,tot));factors=np.array([chs[x['channel']]['phase_factor']/(td[x['total_width_key']] if x['observable']=='branching_fraction' else 1.) for x in training]);rows=np.array([chs[k]['row'] for k in labels]);off=np.array([chs[k]['offset'] for k in labels]);return sec,kk,chs,rows,off,factors
            zstart=np.r_[ef.x,np.zeros(len(tk))];sec,kk,chs,rows,off,factors=blocks(zstart)
            for model,mapping in maps.items():
                gauge,rk=image_gauge(rows,mapping);observed_defects=[];A=np.sqrt(factors)[:,None]*(rows@gauge);b=np.sqrt(factors)*off;starts,bestscan,nscan=signed_starts(A,b,np.sqrt(y),errors,args.sign_starts)
                saturated=rk==len(training);solutions=[]
                lo=np.r_[lower,-tv/ts+1e-8,np.full(rk,-np.inf)];hi=np.r_[upper,np.full(len(tk)+rk,np.inf)]
                def residual(v):
                    try:
                        sec,kk,chs,rr,oo,ff=blocks(v);observed_defects.append(containment(rr,mapping,gauge));amp=oo+rr@gauge@v[npos+len(tk):];odd=(ff*amp**2-y)/errors
                        return np.r_[odd,(J.even_observables(sec,kk,ass)-even_target)/even_sigma,v[npos:npos+len(tk)]]
                    except (ValueError,np.linalg.LinAlgError):return np.full(len(training)+7+len(tk),1e12)
                if saturated:
                    coefficient=la.solve(rows@gauge,np.sqrt(y/factors)*np.sign(A@starts[0]+b)-off)
                    v=np.r_[zstart,coefficient];res=residual(v);solutions=[(float(res@res),v,ef.success,ef.nfev,'exact_odd_interpolation_plus_joint_ancillary_optimum')]
                else:
                    for start in starts:
                        fit=least_squares(residual,np.r_[zstart,start],bounds=(lo,hi),loss='linear',xtol=1e-9,ftol=1e-9,gtol=1e-7,max_nfev=180)
                        solutions.append((float(fit.fun@fit.fun),fit.x,fit.success,fit.nfev,str(fit.message)))
                solutions.sort(key=lambda x:x[0]);chi,v,success,nfev,message=solutions[0];sec,kk,chs,rr,oo,ff=blocks(v);amp=oo+rr@gauge@v[npos+len(tk):];pred=ff*amp**2;res=residual(v);defect=containment(rr,mapping,gauge);maxdefect=max(maxdefect,defect)
                heldresults=[]
                for x in held:
                    h=chs[x['channel']]['row'];full=rr@mapping;u,sv,vh=la.svd(full,full_matrices=True);rank=J.rank(full)[0];null=vh[rank:].T;null_response=np.linalg.norm(h@mapping@null);identifiable=null_response<1e-8*max(1.,np.linalg.norm(h@mapping))
                    heldresults.append(dict(channel=x['channel'],odd_amplitude_identifiable_conditional_on_fixed_nuisances=bool(identifiable),fixed_nuisance_odd_null_response=float(null_response),joint_nuisance_identifiability='not assessed',prediction_status='withheld until shared-width and null-direction propagation' if identifiable else 'unidentifiable; no gauge forecast'))
                physical_coeff=gauge@v[npos+len(tk):];fitted_J=2*ff[:,None]*amp[:,None]*(rr@gauge)/errors[:,None];sv=la.svdvals(fitted_J)
                results.append(dict(scenario_index=si,scenario={k:s[k] for k in ('rho_input','kstar_input','phi_input','branch','source_hypothesis')},control=control,model=model,measurement_labels=labels,visible_start_rank=rk,visible_fit_rank=J.rank(rr@mapping)[0],coefficient_gauge_parent_matrix=gauge.tolist(),gauge_representative_coefficients=physical_coeff.tolist(),gauge_containment_relative_defect=defect,maximum_evaluated_gauge_containment_defect=max([defect,*observed_defects]),action_image_certified=bool(max([defect,*observed_defects])<1e-8),scientific_status='valid_conditional_action_image' if max([defect,*observed_defects])<1e-8 else 'invalid_fixed_gauge',adequacy_inference_allowed=bool(max([defect,*observed_defects])<1e-8 and not saturated),gauge_is_not_null_direction_measurement=True,chi2=chi,odd_chi2=float(res[:len(training)]@res[:len(training)]),ancillary_chi2=float(res[len(training):]@res[len(training):]),nominal_dof=len(training)+7-npos-rk,nominal_dof_caveat='conditional diagonal summaries, local rank; saturation is not adequacy',success=bool(success),nfev=int(nfev),message=message,active_search_bounds=np.flatnonzero((np.abs(v-lo)<1e-7)|(np.abs(v-hi)<1e-7)).tolist(),predictions=pred.tolist(),pulls=((pred-y)/errors).tolist(),physical_amplitudes=amp.tolist(),odd_whitened_singular_values=sv.tolist(),even_parameters=(point+scales*v[:npos]).tolist(),shared_total_widths_gev=dict(zip(tk,(tv+ts*v[npos:npos+len(tk)]).tolist())),enumerated_start_signs=nscan,retained_nonlinear_sign_candidates=len(solutions),sign_minima=[dict(chi2=q[0],success=bool(q[2]),amplitude_signs=np.sign(blocks(q[1])[4]+blocks(q[1])[3]@gauge@q[1][npos+len(tk):]).astype(int).tolist()) for q in solutions],global_nonlinear_minimum_proven=False,held_out=heldresults))
        print(f'scenario {si+1}/64 complete',flush=True)
    inputs=list(artifact['sha256'])+['data/vector_sector/odd_low_energy_jacobians.json','scripts/calibrate_odd_low_energy.py']
    out=dict(scope='conditional low-energy calibration; no model-selection or ALP uncertainty band',complete_64_scenario_run=not args.limit_scenarios,expected_fit_count=64*len(MODELS)*len(CONTROLS),fit_count=len(results),models=list(MODELS),controls=list(CONTROLS),covariance_policy='fixed diagonal PDG summaries; cross covariance unavailable; each original mass/electronic/shared-total observation once',coefficient_policy='fixed pivot representative gauge; all omitted action directions explicitly unmeasured',mass_policy=artifact['vertex_scalar_mass_scheme'],phase_space_policy=artifact['phase_space_scalar_mass_scheme'],search_bounds='m_rho [0.65,0.90], m_Kstar [0.75,1.05] GeV; Delta_V [-1,1] GeV^2; kappa_V [-2,2]; positive total widths; sign(kappa_V-4*(m_Kstar^2-m_rho^2)/(3*(xl-xs))) retained from start; root merger excluded; source Fl>0 and declared Fs sign; doubletrace physical residue signs of positive universal overlaps; numerical domains, not priors',limitations=['local nonlinear optimization, exhaustive fixed-start sign scan only','Alda mixing/quark-spurion covariance fixed','coefficient representatives and topology decomposition gauge dependent','saturated fits cannot establish truncation adequacy','withheld identifiable channels need shared-width and null uncertainty propagation'],symbolic_derivation=symbolic_derivation(),maximum_gauge_containment_defect=maxdefect,results=results,sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    path=ROOT/'data/vector_sector'/('odd_low_energy_calibration_debug.json' if args.limit_scenarios else 'odd_low_energy_calibration.json');path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(path.relative_to(ROOT))

if __name__=='__main__':main()
