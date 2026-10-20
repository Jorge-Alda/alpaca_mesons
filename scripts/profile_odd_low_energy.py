"""Expanded conditional sign searches and nuisance profiles (GeV).

Run .venv/bin/python scripts/profile_odd_low_energy.py.
Every original mass/electronic/shared-width observation occurs once. No new
prior, OPE condition, covariance approximation or ALP prediction is introduced.
"""
from pathlib import Path
import argparse,hashlib,itertools,json,math,time
import numpy as np
import scipy.linalg as la
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_distribution
import sympy as sp
import calibrate_odd_low_energy as C
import verify_odd_low_energy_calibration as V
ROOT=C.ROOT

class ConditionalLikelihood:
    def __init__(self,result,data,pdg,ass,fast):
        self.result,self.ass,self.fast=result,ass,fast
        self.mapping=C.model_maps()['soft_leading'];self.original=C.model_maps()['original_11'];self.point=np.asarray(result['even_parameters']);self.hyp=result['scenario']['source_hypothesis'];self.np=len(self.point)
        self.scales=np.array([.001,.001,.01,.01,.01]+([] if self.hyp=='universal' else [.1])+([.1] if self.hyp=='double_trace' else []))
        self.obs=[x for x in data['observations'] if x.get('status','included')!='deferred' and x['channel'] in result['measurement_labels']];self.labels=[x['channel'] for x in self.obs];self.y,self.errors=np.array([C.J.data_units(x) for x in self.obs]).T
        allobs=[x for x in data['observations'] if x.get('status','included')!='deferred'];self.tk=sorted(set(x['total_width_key'] for x in allobs if x['observable']=='branching_fraction'));self.tv=np.array([data['total_widths'][k]['value_gev'] for k in self.tk]);self.ts=np.array([data['total_widths'][k]['sigma_gev'] for k in self.tk]);self.nn=self.np+len(self.tk)
        s=result['scenario'];self.even_y=np.array([pdg['masses'][k]['value_mev']/1000 for k in (s['rho_input'],s['kstar_input'],'omega','phi')]+[pdg['leptonic_widths'][k]['value_kev']*1e-6 for k in ('rho','omega','phi')]);self.even_errors=np.array([pdg['masses'][k]['sigma_mev']/1000 for k in (s['rho_input'],s['kstar_input'],'omega','phi')]+[pdg['leptonic_widths'][k]['sigma_kev']*1e-6 for k in ('rho','omega','phi')])
        if s['phi_input']!='phi_evaluation':self.even_y[-1]=pdg['alternatives']['phi_ee_average']['value_kev']*1e-6;self.even_errors[-1]=pdg['alternatives']['phi_ee_average']['sigma_kev']*1e-6
        midpoint=4*(self.point[1]**2-self.point[0]**2)/(3*(ass['xl_gev2']-ass['xs_gev2']));self.root_sign=np.sign(self.point[3]-midpoint);self.cache={};self.maximum_image_defect=0.;self.domain_rejections=0
        self.gauge=None;self.old_coefficients=np.asarray(result['gauge_representative_coefficients']);temp=np.r_[np.zeros(self.np),[(result['shared_total_widths_gev'][k]-self.tv[i])/self.ts[i] for i,k in enumerate(self.tk)]]
        self.lower=np.r_[np.full(self.np,-np.inf),-self.tv/self.ts+1e-8];self.upper=np.full(self.nn,np.inf);self.lower[:2]=(np.array([.65,.75])-self.point[:2])/self.scales[:2];self.upper[:2]=(np.array([.9,1.05])-self.point[:2])/self.scales[:2];self.lower[2:4]=(-np.array([1.,2.])-self.point[2:4])/self.scales[2:4];self.upper[2:4]=(np.array([1.,2.])-self.point[2:4])/self.scales[2:4]
        sec,mk,channels,rows,off,factors=self.blocks(temp);self.gauge,self.rank=C.image_gauge(rows,self.mapping);self.lower=np.r_[self.lower,np.full(self.rank,-np.inf)];self.upper=np.r_[self.upper,np.full(self.rank,np.inf)]
        odd=la.lstsq(rows@self.gauge,rows@self.old_coefficients,cond=1e-10)[0];self.start=np.r_[temp,odd]
    def sector(self,z):
        p=self.point+self.scales*z[:self.np];sec,mk=C.J.action_sector(p,self.hyp,self.ass);midpoint=4*(p[1]**2-p[0]**2)/(3*(self.ass['xl_gev2']-self.ass['xs_gev2']))
        if self.root_sign*(p[3]-midpoint)<=0:raise ValueError('inverse mass branch crossed')
        if self.hyp=='universal' and p[4]<=0:raise ValueError('universal FV sign crossed')
        if self.hyp in ('single_trace_positive','single_trace_negative'):
            fl=p[4]+2*p[5]*self.ass['xl_gev2'];fs=p[4]+2*p[5]*self.ass['xs_gev2']
            if fl<=0 or (fs>0)!=(self.hyp=='single_trace_positive'):raise ValueError('single-trace source signs crossed')
        # The double-trace model adds an independent spurion source term.
        # Its signed residues are predictions of (F_V,z_V,y_V), not a
        # positivity constraint relative to the universal-source residues.
        return sec,mk
    def blocks(self,v):
        if np.any(v[:self.nn]<self.lower[:self.nn]) or np.any(v[:self.nn]>self.upper[:self.nn]):raise ValueError('numerical domain bounds crossed')
        key=np.asarray(v[:self.np]).tobytes()
        if key not in self.cache:
            sec,mk=self.sector(v);channels=self.fast(sec,mk);rows=np.array([channels[k]['row'] for k in self.labels]);original_gauge,_=C.image_gauge(rows,self.original);leading_gauge,_=C.image_gauge(rows,self.mapping);defects=[C.containment(rows,self.original,leading_gauge),C.containment(rows,self.mapping,original_gauge)]
            if self.gauge is not None:defects.append(C.containment(rows,self.mapping,self.gauge))
            self.maximum_image_defect=max(self.maximum_image_defect,*defects)
            if max(defects)>1e-8:raise ValueError('training image certificate failed')
            self.cache.clear();self.cache[key]=(sec,mk,channels,rows)
        sec,mk,channels,rows=self.cache[key];td=dict(zip(self.tk,self.tv+self.ts*v[self.np:self.nn]));off=np.array([channels[k]['offset'] for k in self.labels]);factors=np.array([channels[x['channel']]['phase_factor']/(td[x['total_width_key']] if x['observable']=='branching_fraction' else 1.) for x in self.obs])
        if np.any(factors<=0):raise ValueError('closed channel or nonpositive total width')
        return sec,mk,channels,rows,off,factors
    def residual(self,v):
        try:
            sec,mk,channels,rows,off,factors=self.blocks(v);amp=off+rows@self.gauge@v[self.nn:]
            return np.r_[(factors*amp**2-self.y)/self.errors,(C.J.even_observables(sec,mk,self.ass)-self.even_y)/self.even_errors,v[self.np:self.nn]]
        except (ValueError,np.linalg.LinAlgError):self.domain_rejections+=1;return np.full(len(self.obs)+7+len(self.tk),1e12)
    def predictions(self,v):
        sec,mk,channels,rows,off,factors=self.blocks(v);amp=off+rows@self.gauge@v[self.nn:];return factors*amp**2,amp
    def rebase(self,result):
        p=np.asarray(result['even_parameters']);nuis=np.r_[(p-self.point)/self.scales,[(result['shared_total_widths_gev'][k]-self.tv[i])/self.ts[i] for i,k in enumerate(self.tk)]];rows=self.blocks(nuis)[3];coef=np.asarray(result['gauge_representative_coefficients']);odd=la.lstsq(rows@self.gauge,rows@coef,cond=1e-10)[0];return np.r_[nuis,odd]
    def fit(self,start,max_nfev=220):
        fit=least_squares(self.residual,start,bounds=(self.lower,self.upper),loss='linear',xtol=1e-10,ftol=1e-10,gtol=1e-7,max_nfev=max_nfev)
        return fit
    def sign_starts(self,count):
        sec,mk,channels,rows,off,factors=self.blocks(self.start);A=np.sqrt(factors)[:,None]*(rows@self.gauge);b=np.sqrt(factors)*off
        starts,score,n=C.signed_starts(A,b,np.sqrt(self.y),self.errors,count)
        return [np.r_[self.start[:self.nn],x] for x in starts],score,n
    def summary(self,fit):
        predictions,amp=self.predictions(fit.x);n=len(self.obs);res=self.residual(fit.x);active=np.flatnonzero((np.abs(fit.x-self.lower)<1e-7)|(np.abs(fit.x-self.upper)<1e-7)).tolist();params=(self.point+self.scales*fit.x[:self.np]).tolist();source=V.source_certificate(params,self.hyp,self.ass);branch=V.branch_certificate(params,self.hyp,self.result['scenario']['branch'],self.ass)
        return dict(chi2=float(res@res),odd_chi2=float(res[:n]@res[:n]),ancillary_chi2=float(res[n:]@res[n:]),success=bool(fit.success),nfev=int(fit.nfev),message=str(fit.message),active_search_bounds=active,source_certificate=source,branch_certificate=branch,predictions=predictions.tolist(),physical_amplitudes=amp.tolist(),amplitude_signs=np.sign(amp).astype(int).tolist(),pulls=((predictions-self.y)/self.errors).tolist(),even_parameters=params,shared_total_widths_gev=dict(zip(self.tk,(self.tv+self.ts*fit.x[self.np:self.nn]).tolist())),gauge_representative_coefficients=(self.gauge@fit.x[self.nn:]).tolist(),coordinates=fit.x.tolist())

def deterministic_eigenvectors(covariance,count=2,tolerance=1e-8):
    """Largest eigenspaces; canonical projected axes resolve degeneracy."""
    values,vectors=la.eigh(covariance);order=np.argsort(values)[::-1];values=values[order];vectors=vectors[:,order];chosen=[];i=0
    while len(chosen)<count and i<len(values):
        group=np.flatnonzero(np.abs(values-values[i])<=tolerance*max(1.,abs(values[i])));group=group[group>=i];basis=vectors[:,group];projector=basis@basis.T;local=[]
        for axis in range(len(values)):
            v=projector[:,axis].copy()
            for old in local:v-=old*(old@v)
            norm=np.linalg.norm(v)
            if norm>1e-8:
                v/=norm;first=np.flatnonzero(np.abs(v)>1e-10)[0]
                if v[first]<0:v=-v
                local.append(v);chosen.append((float(values[i]),v,len(group)))
                if len(local)==len(group) or len(chosen)==count:break
        i=int(group[-1])+1
    return chosen

def identifiable_covariance(J,H,tolerance=1e-10):
    norms=np.linalg.norm(J,axis=0);keep=norms>1e-12*max(1.,norms.max(initial=0));B=J[:,keep]/norms[keep];T=H[:,keep]/norms[keep];u,sv,vh=la.svd(B,full_matrices=False);rank=int(np.sum(sv>tolerance*max(1.,sv[0])));inverse=(vh[:rank].T/sv[:rank])@u[:,:rank].T;defect=float(max(np.linalg.norm(T-T@inverse@B)/max(1.,np.linalg.norm(T)),np.linalg.norm(H[:,~keep])/max(1.,np.linalg.norm(H))))
    return (T@inverse)@(T@inverse).T,rank,defect,sv

def amplitude_diagnostics(likelihood,v,mapping=None):
    mapping=likelihood.mapping if mapping is None else mapping;np_=likelihood.np;nn=likelihood.nn;n=len(likelihood.obs);sec,mk,channels,rows,off,factors=likelihood.blocks(v);coef=likelihood.gauge@v[nn:];G=off+rows@coef;odd=2*factors[:,None]*G[:,None]*(rows@mapping)/likelihood.errors[:,None];Jodd=np.vstack([odd,np.zeros((14,mapping.shape[1]))]);Hodd=rows@mapping;nuisance=[];Hnuisance=[];stepdefect=0.
    def values(point):
        sec,mk,chs,rr,oo,ff=likelihood.blocks(point);amp=oo+rr@coef;return np.r_[ff*amp**2/likelihood.errors,C.J.even_observables(sec,mk,likelihood.ass)/likelihood.even_errors,point[np_:nn]],amp
    for i in range(nn):
        derivatives=[]
        for step in (1e-4,5e-5):
            hi=v.copy();lo=v.copy();hi[i]+=step;lo[i]-=step
            if i>=np_:
                k=i-np_;jr=np.zeros(n+14)
                for j,item in enumerate(likelihood.obs):
                    if item.get('total_width_key')==likelihood.tk[k]:jr[j]=-factors[j]*G[j]**2*likelihood.ts[k]/(likelihood.tv[k]+likelihood.ts[k]*v[i])/likelihood.errors[j]
                jr[n+7+k]=1.;derivatives.append((jr,np.zeros(n)))
            else:
                up=values(hi);down=values(lo);derivatives.append(((up[0]-down[0])/(2*step),(up[1]-down[1])/(2*step)))
        stepdefect=max(stepdefect,np.linalg.norm(derivatives[0][0]-derivatives[1][0])/max(1.,np.linalg.norm(derivatives[1][0])));nuisance.append(derivatives[1][0]);Hnuisance.append(derivatives[1][1])
    joint=np.column_stack([Jodd,np.array(nuisance).T]);gradient=np.column_stack([Hodd,np.array(Hnuisance).T]);cov,rank,defect,sv=identifiable_covariance(joint,gradient)
    return dict(physical_amplitude_covariance=cov.tolist(),physical_amplitude_null_relative_defect=defect,joint_rank=rank,derivative_step_relative_change=stepdefect,balanced_singular_values=sv.tolist(),physical_amplitude_covariance_valid=bool(defect<1e-8 and stepdefect<1e-6),phase_space_and_shared_width_correlations_retained=True)


def constraint_coordinates(likelihood,v,c,Gref,target):
    rows,off=likelihood.blocks(v)[3:5];h=c@(rows@likelihood.gauge);shift=c@(off-Gref);norm=h@h
    if norm<1e-20:raise ValueError('physical contrast lost odd-coordinate sensitivity')
    Q=la.qr(h[:,None],mode='full')[0];null=Q[:,1:];part=h*(target-shift)/norm
    return part,null

def fit_contrast(likelihood,c,Gref,target,starts):
    candidates=[]
    for start in starts:
        try:part,null=constraint_coordinates(likelihood,start,c,Gref,target);z=np.r_[start[:likelihood.nn],null.T@(start[likelihood.nn:]-part)]
        except ValueError:continue
        def expand(z):
            v=np.r_[z[:likelihood.nn],np.zeros(likelihood.rank)];part,null=constraint_coordinates(likelihood,v,c,Gref,target);return np.r_[v[:likelihood.nn],part+null@z[likelihood.nn:]]
        def residual(z):
            try:return likelihood.residual(expand(z))
            except ValueError:return np.full(len(likelihood.obs)+14,1e12)
        lower=np.r_[likelihood.lower[:likelihood.nn],np.full(likelihood.rank-1,-np.inf)];upper=np.r_[likelihood.upper[:likelihood.nn],np.full(likelihood.rank-1,np.inf)]
        fit=least_squares(residual,z,bounds=(lower,upper),loss='linear',xtol=1e-9,ftol=1e-9,gtol=1e-7,max_nfev=220)
        try:
            v=expand(fit.x);pred,amp=likelihood.predictions(v);defect=float(abs(c@(amp-Gref)-target));chi=float(likelihood.residual(v)@likelihood.residual(v));valid=bool(fit.success and defect<1e-7*max(1.,abs(target)) and chi<1e20)
            candidates.append(dict(chi2=chi,success=bool(fit.success),valid=valid,nfev=int(fit.nfev),constraint_absolute_defect=defect,coordinates=v.tolist(),amplitude_signs=np.sign(amp).astype(int).tolist()))
        except ValueError:continue
    candidates.sort(key=lambda x:(not x['valid'],x['chi2']))
    return candidates[0] if candidates else dict(valid=False,success=False,chi2=None,reason='no feasible profile start')

def closest_source_quantity(likelihood,reference):
    p=likelihood.point+likelihood.scales*reference[:likelihood.np];hyp=likelihood.hyp
    if hyp=='universal':return ('F_V',4,1.)
    if hyp in ('single_trace_positive','single_trace_negative'):
        fl=p[4]+2*p[5]*likelihood.ass['xl_gev2'];sign=1. if hyp=='single_trace_positive' else -1.;fs=sign*(p[4]+2*p[5]*likelihood.ass['xs_gev2']);return ('F_l' if fl<=fs else 'signed_F_s',4,sign if fs<fl else 1.)
    sec,mk=likelihood.sector(reference);universal=C.J.module('vector_sector').source_residues(sec.generators,np.diag([likelihood.ass['xl_gev2']]*2+[likelihood.ass['xs_gev2']]),1.);signs=np.sign(universal);i=int(np.argmin(np.asarray(sec.residues)*signs));weights=[]
    for j in range(4,likelihood.np):
        q=np.zeros(likelihood.np);q[:4]=p[:4];q[j]=1.;tmp,_=C.J.action_sector(q,hyp,likelihood.ass);weights.append(abs(tmp.residues[i]*likelihood.scales[j]))
    return ('signed_residue_'+str(i),4+int(np.argmax(weights)),float(signs[i]))

def source_quantity(likelihood,v,descriptor):
    name,eliminated,sign=descriptor;p=likelihood.point+likelihood.scales*v[:likelihood.np]
    if name=='F_V':return float(p[4])
    if name=='F_l':return float(p[4]+2*p[5]*likelihood.ass['xl_gev2'])
    if name=='signed_F_s':return float(sign*(p[4]+2*p[5]*likelihood.ass['xs_gev2']))
    sec,mk=C.J.action_sector(p,likelihood.hyp,likelihood.ass);return float(sign*sec.residues[int(name.split('_')[-1])])

def boundary_profiles(likelihood,reference,competing):
    descriptor=closest_source_quantity(likelihood,reference);index=descriptor[1];mask=np.arange(len(reference))!=index;baseline=source_quantity(likelihood,reference,descriptor);points=[]
    for multiplier in (.05,.25,.5,1.,2.,4.,8.):
        target=max(1e-14,baseline*multiplier);candidates=[]
        def expand(z):
            v=np.zeros(len(reference));v[mask]=z;low=source_quantity(likelihood,v,descriptor);v[index]=1.;slope=source_quantity(likelihood,v,descriptor)-low
            if abs(slope)<1e-14:raise ValueError('source constraint lost sensitivity')
            v[index]=(target-low)/slope;return v
        def residual(z):
            try:return likelihood.residual(expand(z))
            except ValueError:return np.full(len(likelihood.obs)+14,1e12)
        for start in (reference,competing):
            fit=least_squares(residual,start[mask],bounds=(likelihood.lower[mask],likelihood.upper[mask]),loss='linear',max_nfev=220,xtol=1e-9,ftol=1e-9,gtol=1e-7)
            try:
                v=expand(fit.x);chi=float(likelihood.residual(v)@likelihood.residual(v));defect=abs(source_quantity(likelihood,v,descriptor)-target);valid=bool(fit.success and chi<1e20 and defect<1e-9*max(baseline,1e-8));candidates.append(dict(valid=valid,success=bool(fit.success),chi2=chi if valid else None,coordinates=v.tolist(),constraint_absolute_defect=float(defect),nfev=int(fit.nfev)))
            except ValueError:continue
        candidates.sort(key=lambda x:(not x['valid'],x['chi2'] if x['chi2'] is not None else float('inf')));point=candidates[0] if candidates else dict(valid=False,chi2=None,reason='no feasible source-constraint point');point.update(multiplier=multiplier,signed_source_target_gev=target);points.append(point)
    return dict(quantity=descriptor[0],reference_signed_source_gev=baseline,positive_numerical_floor_gev=1e-14,one_sided=True,symmetric_gaussian_interval=False,points=points)


def symbolic_profile_derivation():
    q,h1,h2,shift,a1,a2=sp.symbols('q h_1 h_2 b a_1 a_2',real=True);h=sp.Matrix([[h1,h2]]);a=sp.Matrix([a1,a2]);contrast=(h*a)[0]+shift;solution=sp.solve(sp.Eq(contrast,q),a1)[0]
    s=sp.symbols('s',positive=True);sigma,G,f=sp.symbols('sigma G f',positive=True);tau=sp.solve(sp.Eq(2*f*G*s,sigma),s)[0]
    return dict(physical_contrast=sp.latex(contrast),exact_eliminated_coordinate=sp.latex(solution),amplitude_error_scale=sp.latex(tau),linear_gaussian_profile=sp.latex(sp.simplify((q/sp.Symbol('sigma_q',positive=True))**2)))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--limit-pairs',type=int,default=0);parser.add_argument('--search-only',action='store_true');parser.add_argument('--start-pair',type=int,default=0);args=parser.parse_args()
    suffix=f'_{args.start_pair}_{args.limit_pairs}' if args.limit_pairs or args.start_pair else ''
    source=ROOT/'data/vector_sector/odd_low_energy_calibration.json';saved=json.loads(source.read_text());verification=json.loads((ROOT/'data/vector_sector/odd_low_energy_calibration_verification.json').read_text())
    for path,digest in saved['sha256'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:raise ValueError('stale calibration input '+path)
    data=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text());pdg=json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text());ass=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions'];alda=json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text());tree=C.J.module('alda_2026nrd_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F);fast=C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*math.pi*ass['alpha_em']));maps=C.model_maps();results=[]
    pairs=[x for x in saved['results'] if x['model']=='soft_leading'];pairs=pairs[args.start_pair:args.start_pair+(args.limit_pairs or len(pairs))]
    for index,old in enumerate(pairs):
        original=next(x for x in saved['results'] if x['scenario_index']==old['scenario_index'] and x['control']==old['control'] and x['model']=='original_11');likelihood=ConditionalLikelihood(old,data,pdg,ass,fast);starts,score,nscan=likelihood.sign_starts(32);budget_starts=starts[:16]+[likelihood.start,likelihood.rebase(original)];fits=[likelihood.fit(start) for start in budget_starts];fit_summaries=[likelihood.summary(f) for f in fits];oldbest=min(old['chi2'],original['chi2']);firstbest=min(x['chi2'] for x in fit_summaries)
        adaptive=firstbest<oldbest-1e-4
        if adaptive:
            newfits=[likelihood.fit(start) for start in starts[16:32]];fits+=newfits;fit_summaries+=[likelihood.summary(f) for f in newfits]
        ordered=sorted(range(len(fits)),key=lambda i:(not fit_summaries[i]['success'],fit_summaries[i]['chi2']));bestfit=fits[ordered[0]];best=fit_summaries[ordered[0]];merged={**old,**best};local=V.local_diagnostics(merged,data,ass,fast,maps['soft_leading']);df=len(old['measurement_labels'])+14-local['joint_visible_rank'];pvalue=float(chi2_distribution.sf(best['chi2'],df)) if df>0 else None;priority=bool(best['success'] and not best['active_search_bounds'] and not best['source_certificate']['near_domain_boundary'] and pvalue is not None and pvalue>=.05 and max(abs(x) for x in best['pulls'])<=3)
        record=dict(scenario_index=old['scenario_index'],control=old['control'],scenario=old['scenario'],paired_models=['original_11','soft_leading'],training_only_equivalence=True,old_achieved_chi2=oldbest,best=best,fit_count=len(fits),fixed_nuisance_enumerated_sign_count=nscan,adaptive_second_batch=adaptive,second_batch_continued_improvement=bool(adaptive and best['chi2']<firstbest-1e-4),maximum_every_evaluation_image_defect=likelihood.maximum_image_defect,domain_rejections=likelihood.domain_rejections,search_minima=[{k:v for k,v in x.items() if k not in ('coordinates','gauge_representative_coefficients')} for x in fit_summaries],local_diagnostics=local,local_residual_dof=df,descriptive_pvalue=pvalue,profile_priority=priority,global_minimum_proven=False,profiles=[])
        old_boundary=V.source_certificate(old['even_parameters'],old['scenario']['source_hypothesis'],ass)['near_domain_boundary'] or V.source_certificate(original['even_parameters'],original['scenario']['source_hypothesis'],ass)['near_domain_boundary']
        if old_boundary and not args.search_only:record['one_sided_source_profile']=boundary_profiles(likelihood,bestfit.x,fits[ordered[1]].x)
        if priority and not args.search_only:
            sec,mk,channels,rows,off,factors=likelihood.blocks(bestfit.x);Gref=np.asarray(best['physical_amplitudes']);tau=likelihood.errors/(2*factors*np.abs(Gref));ampdiag=amplitude_diagnostics(likelihood,bestfit.x);record['physical_amplitude_diagnostics']=ampdiag;assert ampdiag['physical_amplitude_covariance_valid'],'physical amplitude covariance failed';ampcov=np.asarray(ampdiag['physical_amplitude_covariance']);normalized=ampcov/tau[:,None]/tau[None,:];contrasts=deterministic_eigenvectors(normalized,2)
            competing=next((fits[k].x for k in ordered[1:] if fit_summaries[k]['amplitude_signs']!=best['amplitude_signs']),bestfit.x)
            for eigenvalue,w,multiplicity in contrasts:
                c=w/tau;points=[];previous=bestfit.x;reference_chi=best['chi2']
                for t in np.arange(-4.,4.01,.5):
                    point=fit_contrast(likelihood,c,Gref,float(t*np.sqrt(max(eigenvalue,0))),[previous,bestfit.x,competing]);point['standardized_target']=float(t);point['physical_contrast_target']=float(t*np.sqrt(max(eigenvalue,0)));point['delta_chi2']=point['chi2']-reference_chi if point['valid'] else None;points.append(point)
                    if point['valid']:previous=np.asarray(point['coordinates'])
                for sign in (-1.,1.):
                    if not any(x['valid'] and sign*x['standardized_target']>0 and x['delta_chi2']>=3.84 for x in points):
                        for extension in (6.,8.):
                            t=sign*extension;point=fit_contrast(likelihood,c,Gref,float(t*np.sqrt(max(eigenvalue,0))),[previous,bestfit.x,competing]);point.update(standardized_target=t,physical_contrast_target=float(t*np.sqrt(max(eigenvalue,0))),delta_chi2=point['chi2']-reference_chi if point['valid'] else None);points.append(point)
                            if point['valid'] and point['delta_chi2']>=3.84:break
                record['profiles'].append(dict(eigenvalue=eigenvalue,degenerate_eigenspace_dimension=multiplicity,w=w.tolist(),tau=tau.tolist(),reference_physical_amplitudes=Gref.tolist(),quantity='fixed physical training-amplitude contrast with regenerated nuisance-dependent rows',null_safety_certificate=ampdiag['physical_amplitude_null_relative_defect'],points=points,lower_reference_found=any(x['valid'] and x['delta_chi2']< -1e-4 for x in points),open_interval_sides=[sign for sign in (-1,1) if not any(x['valid'] and sign*x['standardized_target']>0 and x['delta_chi2']>=3.84 for x in points)]))
        # A discovered lower reference remains unresolved until re-anchored;
        # never report a likelihood interval relative to an obsolete minimum.
        record['profile_reference_resolved']=not any(p['lower_reference_found'] for p in record['profiles'])
        for profile in record['profiles']:
            profile['likelihood_level_brackets']={}
            points=sorted(profile['points'],key=lambda p:p['standardized_target'])
            for level in (1.,3.84):
                brackets=[]
                if record['profile_reference_resolved']:
                    for left,right in zip(points,points[1:]):
                        if left['valid'] and right['valid'] and (left['delta_chi2']-level)*(right['delta_chi2']-level)<=0:
                            brackets.append([left['standardized_target'],right['standardized_target']])
                profile['likelihood_level_brackets'][str(level)]=brackets
        results.append(record)
        checkpoint=ROOT/'data/vector_sector'/f'odd_low_energy_profiles_checkpoint{suffix}.json'
        checkpoint.write_text(json.dumps(dict(complete=False,start_pair=args.start_pair,pair_count=len(results),results=results),indent=2,allow_nan=False)+'\n')
        print(f'pair {index+1}/{len(pairs)} complete',flush=True)
    inputs=['data/vector_sector/odd_low_energy_calibration.json','data/vector_sector/odd_low_energy_calibration_verification.json','scripts/profile_odd_low_energy.py','scripts/calibrate_odd_low_energy.py','scripts/verify_odd_low_energy_calibration.py','data/odd_parity/odd_low_energy_measurements.json']
    output=dict(scope='finite expanded sign search and physical training-contrast profiles; conditional diagonal summaries',complete_256_pair_run=not args.limit_pairs and args.start_pair==0,pair_count=len(results),covered_unsaturated_model_cases=2*len(results),sign_start_policy='16 distinct fixed-point targets plus both old model minima; next16 if firstbatch improves oldchi2 by >1e-4',profile_priority_policy='postsearch interior converged p>=.05 and maxoddpull<=3; computational priority only, not model-selection rule',profile_grid_policy='two largest standardized amplitude-covariance eigendirections; t=-4:0.5:4; missing brackets extended6,8sigma',global_nonlinear_minimum_proven=False,symbolic_derivation=symbolic_profile_derivation(),results=results,sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    path=ROOT/'data/vector_sector'/f'odd_low_energy_profiles{suffix}.json';path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');print(path.relative_to(ROOT))

if __name__=='__main__':main()
