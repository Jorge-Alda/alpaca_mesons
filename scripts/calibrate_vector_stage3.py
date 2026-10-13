"""Conditional stage-3 mass/residue diagnostics, not a precision global fit.

Run: .venv/bin/python scripts/calibrate_vector_stage3.py
GeV internally, keV for electronic widths; independent summary errors assumed.
Charged/neutral mass scenarios and phi width alternatives remain separate.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp
from scipy.optimize import least_squares
from numpy.polynomial import Polynomial
from derive_vector_mass_extension import derive as derive_mass

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'data/vector_sector/pdg2024.json'
OUTPUT = ROOT / 'data/vector_sector/stage3_results.json'
XL, XS = 4*(6568+11740)/2*1e-6, 4*235921*1e-6
ME = 0.00051099895  # fixed electron-mass assumption, GeV
SYMBOLIC = derive_mass()
MASS_KERNEL = sp.lambdify(tuple(SYMBOLIC[k] for k in ('r','k','Delta','kap','xl','xs')),
                          SYMBOLIC['calibrated'], 'numpy')
DETERMINANT_COEFFICIENTS = sp.lambdify(tuple(SYMBOLIC[k] for k in ('r','k','W','P','xl','xs')),
    list(reversed(SYMBOLIC['polynomial'].all_coeffs())), 'numpy')


def mass_matrix(mrho, mkstar, delta, kappa=0.):
    return np.asarray(MASS_KERNEL(mrho**2,mkstar**2,delta,kappa,XL,XS),dtype=float)


def angle(matrix):
    return .5*np.arctan2(2*matrix[0,1],matrix[1,1]-matrix[0,0])


def predicted_masses(parameters):
    mr, mk, delta = parameters
    eigen = np.linalg.eigvalsh(mass_matrix(mr,mk,delta))
    if np.any(eigen <= 0):
        return np.full(4,10.)
    return np.r_[mr,mk,np.sqrt(eigen)]


def mass_diagnostic(observed, errors):
    starts = [np.r_[observed[:2],d] for d in (-.03,.02,.08)]
    fits = [least_squares(lambda p:(predicted_masses(p)-observed)/errors,
                         p,bounds=([.5,.7,-.25],[.95,1.,.25]),
                         xtol=1e-13,ftol=1e-13,gtol=1e-13) for p in starts]
    fit = min(fits,key=lambda f:np.dot(f.fun,f.fun))
    assert fit.success
    mr,mk,mo,mp=observed
    # Fix the trace to both physical isoscalar masses; compare determinants.
    delta_trace=mo**2+mp**2-2*mk**2
    matrix=mass_matrix(mr,mk,delta_trace)
    return dict(parameters=fit.x.tolist(),predicted_mev=(1000*predicted_masses(fit.x)).tolist(),
                residual_mev=(1000*(predicted_masses(fit.x)-observed)).tolist(),
                conditional_chi2=float(fit.fun@fit.fun),nominal_dof=1,
                theta_deg=float(np.degrees(angle(mass_matrix(*fit.x)))),
                trace_delta_gev2=float(delta_trace),
                trace_anchored_predicted_mev=(1000*np.sqrt(np.linalg.eigvalsh(matrix))).tolist(),
                determinant_residual_gev4=float(np.linalg.det(matrix)-mo**2*mp**2))


def extended_branches(observed):
    mr,mk,mo,mp=observed
    t=mo**2+mp**2-2*mk**2
    delta=Polynomial([t,-(2*XL+XS)])
    roots=Polynomial(DETERMINANT_COEFFICIENTS(mr**2,mk**2,mo**2,mp**2,XL,XS)).roots()
    output=[]
    for k in sorted(roots,key=lambda r:float(np.real(r))):
        if abs(np.imag(k))>1e-9:
            continue
        k=float(np.real(k)); shift=float(delta(k))
        matrix=mass_matrix(mr,mk,shift,k)
        np.testing.assert_allclose(np.linalg.eigvalsh(matrix),[mo**2,mp**2],rtol=1e-10)
        output.append(dict(kappa=k,delta_gev2=shift,theta_deg=float(np.degrees(angle(matrix))),
                           exact_mass_interpolation=True))
    return output


def propagated_errors(function, point, errors):
    """Local delta method; input covariance assumed diagonal, no model errors."""
    columns=[]
    for i,error in enumerate(errors):
        step=max(float(error)*.01,1e-9)
        high,low=point.copy(),point.copy()
        high[i]+=step;low[i]-=step
        columns.append((np.asarray(function(high))-np.asarray(function(low)))/(2*step))
    jacobian=np.column_stack(columns)
    covariance=(jacobian*errors)@(jacobian*errors).T
    return np.sqrt(np.diag(covariance)).tolist()


def overlaps(theta):
    return np.array([1/np.sqrt(2), np.cos(theta)/(3*np.sqrt(2))+np.sin(theta)/3,
                     np.sin(theta)/(3*np.sqrt(2))-np.cos(theta)/3])


def width(residues,masses,alpha):
    ratio=(ME/masses)**2
    return 1e6*8*np.pi*alpha**2*residues**2/(3*masses)*(1+2*ratio)*np.sqrt(1-4*ratio)


def fixed_angle_residues(theta,masses,values,errors,alpha):
    coefficient=width(overlaps(theta),masses,alpha)
    fv2=np.sum(coefficient*values/errors**2)/np.sum(coefficient**2/errors**2)
    prediction=coefficient*fv2
    universal=dict(FV_gev=float(np.sqrt(fv2)),prediction_kev=prediction.tolist(),
                   conditional_chi2=float(np.sum(((prediction-values)/errors)**2)),nominal_dof=2,
                   separate_effective_FV_gev=np.sqrt(values/coefficient).tolist())
    def predict_source(p):
        fl,fs=p
        residues=np.array([fl/np.sqrt(2),fl*np.cos(theta)/(3*np.sqrt(2))+fs*np.sin(theta)/3,
                           fl*np.sin(theta)/(3*np.sqrt(2))-fs*np.cos(theta)/3])
        return width(residues,masses,alpha)
    # Positive light/strange residue branch continuously connected to FV>0.
    fit=least_squares(lambda p:(predict_source(p)-values)/errors,[.16,.16],
                      bounds=([.01,.01],[.5,.5]),xtol=1e-13,ftol=1e-13,gtol=1e-13)
    assert fit.success
    fl,fs=fit.x; z=(fs-fl)/(2*(XS-XL)); fv=fl-2*z*XL
    single=dict(Fl_gev=float(fl),Fs_gev=float(fs),FV_gev=float(fv),z_gev_inverse=float(z),
                prediction_kev=predict_source(fit.x).tolist(),conditional_chi2=float(fit.fun@fit.fun),nominal_dof=1)
    alternative=least_squares(lambda p:(predict_source(p)-values)/errors,[.16,-.16],
        bounds=([.01,-.5],[.5,-.01]),xtol=1e-13,ftol=1e-13,gtol=1e-13)
    assert alternative.success
    fl_alt,fs_alt=alternative.x;z_alt=(fs_alt-fl_alt)/(2*(XS-XL))
    negative=dict(Fl_gev=float(fl_alt),Fs_gev=float(fs_alt),FV_gev=float(fl_alt-2*z_alt*XL),
        z_gev_inverse=float(z_alt),prediction_kev=predict_source(alternative.x).tolist(),
        conditional_chi2=float(alternative.fun@alternative.fun),nominal_dof=1,
        sign_assumption='Fl positive, Fs negative; disconnected from small breaking around universal positive FV')
    # Three independent current residues interpolate three positive widths.
    magnitudes=np.sqrt(values/width(np.ones(3),masses,alpha))
    signed=magnitudes*np.sign(overlaps(theta))
    rotation=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    ideal=rotation@signed[1:]
    fl=np.sqrt(2)*signed[0]
    y=(ideal[0]-fl/(3*np.sqrt(2)))/(np.sqrt(2)*(XL-XS)/3)
    fs=-3*(ideal[1]-y*(XL-XS)/3)
    z=(fs-fl)/(2*(XS-XL));fv=fl-2*z*XL
    restored=np.r_[fl/np.sqrt(2),rotation.T@np.array([fl/(3*np.sqrt(2))+np.sqrt(2)*y*(XL-XS)/3,
                                                  -fs/3+y*(XL-XS)/3])]
    np.testing.assert_allclose(width(restored,masses,alpha),values,rtol=1e-12)
    complete=dict(FV_gev=float(fv),z_gev_inverse=float(z),y_gev_inverse=float(y),
                  exact_width_interpolation=True,nominal_dof=0,
                  sign_assumption='same physical-residue signs as universal source at this angle')
    return dict(universal=universal,single_trace_spurion=single,
                single_trace_spurion_sign_alternative=negative,plus_double_trace_spurion=complete)


def main():
    data=json.loads(INPUT.read_text());alpha=data['assumptions']['alpha_em']['value']
    for source in data['sources'].values():
        assert hashlib.sha256((ROOT/source['local_path']).read_bytes()).hexdigest()==source['sha256']
    outcomes=[]
    for rk in ('rho_neutral','rho_charged'):
        for kk in ('kstar_charged','kstar_neutral'):
            rows=[data['masses'][key] for key in (rk,kk,'omega','phi')]
            observed=np.array([r['value_mev'] for r in rows])/1000
            errors=np.array([r['sigma_mev'] for r in rows])/1000
            minimal=mass_diagnostic(observed,errors);branches=extended_branches(observed)
            for i,branch in enumerate(branches):
                def branch_parameters(point):
                    b=extended_branches(point)[i]
                    return [b['kappa'],b['delta_gev2'],b['theta_deg']]
                branch['conditional_input_sigma']=dict(zip(('kappa','delta_gev2','theta_deg'),
                    propagated_errors(branch_parameters,observed,errors)))
            theta_cases=[('minimal_mass_fit',minimal['theta_deg'])]+[(f'extended_mass_branch_{i}',b['theta_deg']) for i,b in enumerate(branches)]
            widths=[]
            for choice in ('phi_evaluation','phi_average'):
                wr=[data['leptonic_widths'][k] for k in ('rho','omega','phi')]
                if choice=='phi_average':wr[-1]=data['alternatives']['phi_ee_average']
                y=np.array([r['value_kev'] for r in wr]);sigma=np.array([r['sigma_kev'] for r in wr])
                # Use physical neutral rho mass for electronic decays in every mass scenario.
                masses=np.array([data['masses'][k]['value_mev'] for k in ('rho_neutral','omega','phi')])/1000
                def free_residual(p):return (width(p[0]*overlaps(p[1]),masses,alpha)-y)/sigma
                free=least_squares(free_residual,[.156,-.06],bounds=([.01,-np.pi/9],[.5,np.pi/9]),
                                   xtol=1e-13,ftol=1e-13,gtol=1e-13)
                assert free.success
                conditional_sigma=np.sqrt(np.diag(np.linalg.inv(free.jac.T@free.jac)))
                widths.append(dict(phi_input=choice,near_ideal_universal_width_only=dict(FV_gev=float(free.x[0]),
                    theta_deg=float(np.degrees(free.x[1])),conditional_chi2=float(free.fun@free.fun),nominal_dof=1,
                    conditional_input_sigma=dict(FV_gev=float(conditional_sigma[0]),theta_deg=float(np.degrees(conditional_sigma[1]))),
                    angle_domain_deg=[-20,20]),
                    fixed_angle_cases=[dict(case=name,theta_deg=theta,**fixed_angle_residues(np.radians(theta),masses,y,sigma,alpha))
                                       for name,theta in theta_cases]))
            outcomes.append(dict(rho_input=rk,kstar_input=kk,minimal_mass=minimal,extended_mass_branches=branches,width_tests=widths))
    result=dict(input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        assumptions=dict(covariance='diagonal summary errors, conditional diagnostics only',
            masses='PDG listing fitted masses used as tree masses; no EM or pole conversion',
            xl_gev2=XL,xs_gev2=XS,spurion='fixed Alda arXiv-v1 central isospin-average inputs',
            model_errors='not included in diagnostic chi2; charged/neutral scenarios and width alternatives remain separate',
            electron_mass_gev=ME,alpha_em=alpha,fit_scope='no odd-parity or three-pion inputs'),scenarios=outcomes)
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(OUTPUT.relative_to(ROOT))
    for o in outcomes:
        print(o['rho_input'],o['kstar_input'],'minimal mass chi2=',round(o['minimal_mass']['conditional_chi2'],2),
              'theta=',round(o['minimal_mass']['theta_deg'],3),
              'extended angles=',[round(b['theta_deg'],3) for b in o['extended_mass_branches']])


if __name__=='__main__':main()
