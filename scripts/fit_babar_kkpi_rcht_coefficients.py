"""Fit RChT coefficient coordinates to BaBar:2007ceh Table 7.

Run: .venv/bin/python scripts/fit_babar_kkpi_rcht_coefficients.py --output /tmp/kkpi_rcht.json
Masses are in GeV; c_i,d_i and delta are dimensionless. Coefficients use the
canonical epsilon convention. The K* is an external stable vector in the
source's quasi-two-body component definition.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares, brentq
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 0.0072973525692838015
GEV_MINUS_TWO_TO_NB = 389379.3656
F = 0.093
OPERATOR_MASS = 0.770
KAON_MASS = 0.493677
PION_MASS = 0.13957039
CHI_PLUS_LIGHT = 0.036616
CHI_PLUS_STRANGE = 0.943684


def load_inputs():
    """Load the existing canonical preliminary representative and BaBar bins."""
    parameters = json.loads((ROOT/'data/vector_sector/odd_preliminary_fit.json').read_text())['parameters']
    table = np.genfromtxt(ROOT/'data/vector_sector/babar_2007_isovector_kstar_k_tau_region.csv', delimiter=',', names=True)
    return parameters, np.column_stack((table['energy_low_GeV'], table['energy_high_GeV'])), table['cross_section_nb'], table['statistical_error_nb']


def kallen(squared_energy, first_squared_mass, second_squared_mass):
    """Kallen function in GeV^4."""
    return (squared_energy-first_squared_mass-second_squared_mass)**2-4*first_squared_mass*second_squared_mass


def rho_inverse_propagator(squared_energy, parameters):
    """Upper-rim rho propagator; P-wave pi-pi running width is an assumption."""
    mass = parameters['m_rho_GeV']
    phase = np.maximum(1-4*PION_MASS**2/np.maximum(np.asarray(squared_energy), 1e-30), 0)**1.5
    pole_phase = (1-4*PION_MASS**2/mass**2)**1.5
    # M Gamma(s)=Gamma(M²)*s/M*[beta(s)/beta(M²)]³.
    imaginary_part = parameters['Gamma_total_rho0_GeV']*np.asarray(squared_energy)/mass*phase/pole_phase
    return mass**2-np.asarray(squared_energy)-1j*imaginary_part


def rcht_current_form_factor(squared_energy, parameters, delta=0., c4=None):
    """F_VP=G/e for Q_I=diag(1/2,-1/2,0), charged K* K channel.

    delta changes c1,c5 by delta/2 and preserves the preliminary real-photon
    amplitudes. Spurions and the pseudoscalar mass inside vertices follow
    the preliminary Alda:2026nrd tree-input prescription, independently of
    physical phase-space masses. F_l is the even-sector isovector source.
    """
    squared_energy = np.asarray(squared_energy)
    vector_mass = parameters['m_Kstar_GeV']
    vertex_pseudoscalar_squared_mass = (CHI_PLUS_LIGHT+CHI_PLUS_STRANGE)/4
    strange_light_difference = (CHI_PLUS_STRANGE-CHI_PLUS_LIGHT)/4
    c1, c5 = parameters['c1']+delta/2, parameters['c5']+delta/2
    c2, c3, c6 = (parameters[key] for key in ('c2','c3','c6'))
    c4 = parameters['c4'] if c4 is None else c4
    direct_numerator = .5*((c1+c2+8*c3-c5)*vertex_pseudoscalar_squared_mass
        +(-c1+c2+c5-2*c6)*vector_mass**2
        +(c1-c2+c5)*squared_energy+8*c4*strange_light_difference)
    exchange_numerator = .5*((parameters['d1']+8*parameters['d2']-parameters['d3'])*vertex_pseudoscalar_squared_mass
        +parameters['d3']*(vector_mass**2+squared_energy))
    light_source = parameters['F_V_GeV']+2*parameters['z_V_GeV_inverse']*CHI_PLUS_LIGHT
    direct = 2*np.sqrt(2)/(F*vector_mass*OPERATOR_MASS)*direct_numerator
    exchange = -4*light_source/(F*vector_mass)*exchange_numerator/rho_inverse_propagator(squared_energy, parameters)
    return direct+exchange


def excited_propagator(squared_energy, mass, width):
    """BaBar rho' running-width ansatz with B_4pi=.64 and power 3/2."""
    squared_energy = np.asarray(squared_energy)
    four_pion_phase = np.maximum(squared_energy-16*PION_MASS**2, 0)**1.5/np.maximum(squared_energy, 1e-30)
    pole_phase = (mass**2-16*PION_MASS**2)**1.5/mass**2
    running_width = width*(.64*four_pion_phase/pole_phase+.36)
    return 1/(mass**2-squared_energy-1j*np.sqrt(np.maximum(squared_energy, 0))*running_width)


def cross_section(energies, parameters, delta=0., c4=None, pole_parameters=None, subtraction=True):
    """Isovector sigma_1 in nb, including the two charge conjugate channels.

    H^mu=i epsilon^{mu nu alpha beta} eps*_nu pV_alpha pK_beta F_VP.
    Single-charge spin/orientation integration gives pi alpha² lambda³/²
    |F_VP|²/(6 s³). Multiplying by two for c.c. yields BaBar Eq.(13)'s
    isovector component sigma_1 (charged and neutral reduced amplitudes
    have equal absolute values and opposite isovector signs).
    """
    squared_energy = np.asarray(energies)**2
    form_factor = rcht_current_form_factor(squared_energy, parameters, delta, c4)
    if pole_parameters is not None:
        residue_real, residue_imaginary, mass, width = pole_parameters
        response = excited_propagator(squared_energy, mass, width)
        if subtraction:
            response -= 1/mass**2
        form_factor = form_factor+(residue_real+1j*residue_imaginary)*response
    discriminant = np.maximum(kallen(squared_energy, parameters['m_Kstar_GeV']**2, KAON_MASS**2), 0)
    return GEV_MINUS_TWO_TO_NB*np.pi*ALPHA**2*discriminant**1.5/(3*squared_energy**3)*np.abs(form_factor)**2


def prepare_bin_prediction(edges, parameters, order=16):
    """Prepare uniform dE quadrature for the published bin estimates."""
    nodes, weights = leggauss(order)
    lower, upper = edges.T
    energies = (lower[:, None]+upper[:, None])/2+(upper-lower)[:, None]*nodes/2
    def predict(delta=0., c4=None, pole_parameters=None, subtraction=True):
        return cross_section(energies, parameters, delta, c4, pole_parameters, subtraction)@weights/2
    return predict


def statistical_residual(prediction, observed, errors, relative_systematics=(), diagonal=False):
    """Profile explicit common nuisances, or use the diagonal sensitivity case."""
    if diagonal:
        errors = np.sqrt(errors**2+sum(value**2 for value in relative_systematics)*prediction**2)
        return (prediction-observed)/errors
    residual = (prediction-observed)/errors
    if not relative_systematics:
        return residual
    columns = prediction[:, None]*np.asarray(relative_systematics)[None, :]/errors[:, None]
    nuisances = -np.linalg.solve(np.eye(len(relative_systematics))+columns.T@columns, columns.T@residual)
    return np.r_[residual+columns@nuisances, nuisances]


def describe_fit(result, names, observed, prediction, lower, upper):
    """Store the achieved coefficient fit and its local numerical diagnostics."""
    singular_values = np.linalg.svd(result.jac*np.maximum(np.abs(result.x), .01), compute_uv=False)
    rank = int(sum(singular_values > singular_values[0]*1e-7))
    active_bounds = [name for name, value, low, high in zip(names, result.x, lower, upper) if min(value-low, high-value) < 1e-5*(high-low)]
    statistic = float(result.fun@result.fun)
    errors = np.sqrt(np.diag(np.linalg.pinv(result.jac.T@result.jac, rcond=1e-13)))
    return dict(parameters=dict(zip(names, map(float,result.x))), chi_squared=statistic,
        nominal_degrees_of_freedom=len(observed)-len(result.x), local_rank=rank,
        probability=float(chi2.sf(statistic,len(observed)-len(result.x))),
        local_errors=None if rank<len(result.x) or active_bounds else dict(zip(names,map(float,errors))),
        scaled_jacobian_singular_values=singular_values.tolist(), parameters_at_bounds=active_bounds,
        predictions_nb=prediction.tolist(), success=bool(result.success), optimizer_message=result.message,
        global_nonlinear_minimum_proven=False,
        probability_qualification='Nominal count; boundary and unresolved-direction distributions are not calibrated')


def profile_crossings(grid, chi_squared, minimum):
    """Interpolate all delta-chi² crossings, retaining disconnected intervals."""
    output = {}
    for level in (1., 3.841458820694124):
        values = np.asarray(chi_squared)-minimum-level
        crossings = []
        for left in range(len(grid)-1):
            if values[left]*values[left+1] < 0:
                crossings.append(float(grid[left]+(grid[left+1]-grid[left])*(-values[left])/(values[left+1]-values[left])))
        output[str(level)] = crossings
    return output


def fit_c4(parameters, edges, observed, errors, relative_systematics=(), diagonal=False):
    """Fit the actual c4 coefficient in the source-inspired six-bin control."""
    predict = prepare_bin_prediction(edges[:6], parameters)
    residual = lambda values: statistical_residual(predict(c4=values[0]),observed[:6],errors[:6],relative_systematics,diagonal)
    starts = [least_squares(residual,[value],bounds=([-.5],[.5]),xtol=1e-12,ftol=1e-12,gtol=1e-12) for value in (-.2,-.05,0.,.05,.2)]
    result = min(starts,key=lambda item:item.fun@item.fun)
    output = describe_fit(result,['c4'],observed[:6],predict(c4=result.x[0]),[-.5],[.5])
    distinct_minima = []
    for candidate in sorted(starts, key=lambda item: item.fun@item.fun):
        if all(abs(candidate.x[0]-entry['c4']) > 1e-5 for entry in distinct_minima):
            distinct_minima.append(dict(c4=float(candidate.x[0]), chi_squared=float(candidate.fun@candidate.fun)))
    output['distinct_minima']=distinct_minima
    grid = np.unique(np.r_[np.linspace(-.25,.25,1001), result.x[0]])
    statistic = [float(residual([value])@residual([value])) for value in grid]
    crossings={}
    for level in (1.,3.841458820694124):
        values=np.asarray(statistic)-output['chi_squared']-level
        crossings[str(level)]=[float(brentq(lambda coefficient: float(residual([coefficient])@residual([coefficient]))-output['chi_squared']-level,grid[left],grid[left+1],xtol=1e-10)) for left in range(len(grid)-1) if values[left]*values[left+1]<0]
    output.update(profile=dict(c4=grid.tolist(),chi_squared=statistic,delta_chi_squared_crossings=profile_crossings(grid,statistic,output['chi_squared'])),
        fitted_coefficients={**{key:parameters[key] for key in [*(f'c{i}' for i in range(1,8)),*(f'd{i}' for i in range(1,5))]},'c4':float(result.x[0])},
        interpretation='BaBar c4 extraction conditional on the remaining preliminary coefficients, stable K*, first six bins, and no excited pole')
    output['profile']['delta_chi_squared_crossings']=crossings
    output['profile']['scan_endpoints_are_confidence_limits']=False
    output['profile']['minimum_chi_squared_minus_reported_minimum']=min(statistic)-output['chi_squared']
    output['systematics_relative']=list(relative_systematics)
    output['systematic_correlation']='diagonal assumption' if diagonal else 'separate fully correlated Gaussian nuisances'
    output['fit_bin_edges_GeV']=edges[:6].tolist()
    return output


def fit_delta(parameters, edges, observed, errors, *, relative_systematics=(), diagonal=False, subtraction=True, real_residue=False, profile=True):
    """Fit c1-c2+c5 through its exact real-photon null coordinate."""
    predict = prepare_bin_prediction(edges,parameters)
    names = ['delta_c1_minus_c2_plus_c5','rho_prime_residue_real_GeV','rho_prime_residue_imaginary_GeV','rho_prime_mass_GeV','rho_prime_width_GeV']
    lower = np.array([-.3,-10.,-10.,1.39,.03]); upper=np.array([.3,10.,10.,1.9,1.5])
    active = [0,1,3,4] if real_residue else list(range(5))
    def expand(values, fixed_delta=None):
        full = np.array([0.,0.,0.,1.5,.3])
        indices = [index for index in active if fixed_delta is None or index != 0]
        full[indices]=values
        if fixed_delta is not None: full[0]=fixed_delta
        return full
    def residual(values,fixed_delta=None):
        full=expand(values,fixed_delta)
        prediction=predict(delta=full[0],pole_parameters=full[1:],subtraction=subtraction)
        return statistical_residual(prediction,observed,errors,relative_systematics,diagonal)
    starts=[]
    for delta in (-.04,0.,.04):
        for mass,width in ((1.46,.2),(1.55,.4),(1.7,.8)):
            for sign in (-1,1):
                initial=np.array([delta,sign*.5,sign*.2,mass,width])
                result=least_squares(residual,initial[active],bounds=(lower[active],upper[active]),x_scale='jac',max_nfev=2000,ftol=1e-10,xtol=1e-10,gtol=1e-10)
                starts.append(result)
    result=min(starts,key=lambda item:item.fun@item.fun)
    best=expand(result.x)
    output=describe_fit(result,[names[index] for index in active],observed,predict(delta=best[0],pole_parameters=best[1:],subtraction=subtraction),lower[active],upper[active])
    output['fixed_rho_prime_imaginary_residue']=0. if real_residue else None
    output['subtracted_excited_pole_at_zero']=subtraction
    output['systematics_relative']=list(relative_systematics)
    output['systematic_correlation']='diagonal assumption' if diagonal else 'separate fully correlated Gaussian nuisances'
    output['fitted_c1_minus_c2_plus_c5']=float(parameters['c1']-parameters['c2']+parameters['c5']+best[0])
    output['fitted_coefficients']={key:float(parameters[key]+(best[0]/2 if key in ('c1','c5') else 0.)) for key in [*(f'c{i}' for i in range(1,8)),*(f'd{i}' for i in range(1,5))]}
    output['number_of_deterministic_starts']=len(starts)
    output['fit_bin_edges_GeV']=edges.tolist()
    if profile:
        grid=np.unique(np.r_[np.linspace(-.12,.12,81),np.linspace(best[0]-.02,best[0]+.02,41),best[0]])
        profile_indices=[index for index in active if index!=0]
        profile_results=[]
        profile_nuisances=[]
        # Independent starts plus nearby optimum avoid a single continuation branch.
        continuation=best[profile_indices]
        def fit_fixed_delta(delta, continuation, additional_starts=()):
            candidates=[]
            alternative_starts = [np.array([sign*.4,0.,1.5,.25]) if not real_residue else np.array([sign*.4,1.5,.25]) for sign in (-1,1)]
            for initial in [continuation,best[profile_indices],*alternative_starts,*additional_starts]:
                item=least_squares(lambda values:residual(values,delta),initial,bounds=(lower[profile_indices],upper[profile_indices]),x_scale='jac',max_nfev=1200,ftol=1e-9,xtol=1e-9,gtol=1e-9)
                candidates.append(item)
            return min(candidates,key=lambda item:item.fun@item.fun)
        for delta in grid:
            chosen=fit_fixed_delta(delta,continuation)
            continuation=chosen.x
            profile_results.append(float(chosen.fun@chosen.fun))
            profile_nuisances.append(chosen.x.copy())
        refined_crossings={}
        changed_brackets=[]
        for level in (1.,3.841458820694124):
            values=np.asarray(profile_results)-output['chi_squared']-level
            crossings=[]
            for left in range(len(grid)-1):
                if values[left]*values[left+1] < 0:
                    def level_residual(delta):
                        fitted=fit_fixed_delta(delta,profile_nuisances[left],(profile_nuisances[left+1],))
                        return float(fitted.fun@fitted.fun)-output['chi_squared']-level
                    if level_residual(grid[left])*level_residual(grid[left+1]) < 0:
                        crossings.append(float(brentq(level_residual,grid[left],grid[left+1],xtol=1e-6)))
                    else:
                        changed_brackets.append(dict(level=level,left=float(grid[left]),right=float(grid[left+1]),reason='Additional endpoint starts found a lower profile branch; coarse crossing removed'))
            refined_crossings[str(level)]=crossings
        output['profile']=dict(delta=grid.tolist(),c1_minus_c2_plus_c5=(grid+parameters['c1']-parameters['c2']+parameters['c5']).tolist(),chi_squared=profile_results,
            delta_chi_squared_crossings=refined_crossings,
            changed_crossing_brackets=changed_brackets,
            minimum_chi_squared_minus_reported_minimum=min(profile_results)-output['chi_squared'],
            scan_endpoints_are_confidence_limits=False)
    return output


def fit_joint_coefficients(parameters, edges, observed, errors):
    """Fit c4 and the new VJP combination simultaneously using BaBar alone.

    The preliminary values fix the other coefficient combinations. Varying
    c4 changes the old radiative fit, whose likelihood is not included here.
    This explicit BaBar-only comparison reveals the coefficient correlation.
    """
    predict=prepare_bin_prediction(edges,parameters)
    names=['delta_c1_minus_c2_plus_c5','c4','rho_prime_residue_real_GeV','rho_prime_residue_imaginary_GeV','rho_prime_mass_GeV','rho_prime_width_GeV']
    lower=np.array([-.3,-.3,-10.,-10.,1.39,.03])
    upper=np.array([.3,.3,10.,10.,1.9,1.5])
    def residual(values):
        return (predict(delta=values[0],c4=values[1],pole_parameters=values[2:])-observed)/errors
    candidates=[]
    for delta in (-.05,0.,.05):
        for c4 in (-.05,.02,.08):
            for sign in (-1,1):
                for mass,width in ((1.48,.2),(1.65,.5)):
                    initial=[delta,c4,.4*sign,.15*sign,mass,width]
                    candidates.append(least_squares(residual,initial,bounds=(lower,upper),x_scale='jac',max_nfev=3000,ftol=1e-10,xtol=1e-10,gtol=1e-10))
    result=min(candidates,key=lambda item:item.fun@item.fun)
    output=describe_fit(result,names,observed,predict(delta=result.x[0],c4=result.x[1],pole_parameters=result.x[2:]),lower,upper)
    covariance=np.linalg.pinv(result.jac.T@result.jac,rcond=1e-13)
    diagonal=np.sqrt(np.diag(covariance))
    correlation=covariance/np.outer(diagonal,diagonal)
    output.update(coefficient_covariance=covariance[:2,:2].tolist(),coefficient_correlation=correlation[:2,:2].tolist(),
        full_covariance=covariance.tolist(),number_of_deterministic_starts=len(candidates),
        fitted_c1_minus_c2_plus_c5=float(parameters['c1']-parameters['c2']+parameters['c5']+result.x[0]),
        fixed_input_qualification='BaBar-only joint c4/delta fit; the remaining preliminary coefficient combinations are held fixed, and their original data likelihood is not included',
        covariance_qualification='Local linear diagnostic only; unreliable if bounds are active or eigenvalues are unresolved')
    output['fit_bin_edges_GeV']=edges.tolist()
    output['fitted_coefficients']={key:float(result.x[1] if key=='c4' else parameters[key]+(result.x[0]/2 if key in ('c1','c5') else 0.)) for key in [*(f'c{i}' for i in range(1,8)),*(f'd{i}' for i in range(1,5))]}
    distinct=[]
    for candidate in sorted(candidates,key=lambda item:item.fun@item.fun):
        if all(np.linalg.norm(candidate.x[:2]-np.array([entry['delta_c1_minus_c2_plus_c5'],entry['c4']])) > 1e-4 for entry in distinct):
            distinct.append(dict(chi_squared=float(candidate.fun@candidate.fun),**dict(zip(names,map(float,candidate.x)))))
    output['distinct_minima']=distinct
    return output


def run_fits():
    """Run coefficient fits and state every fixed physical input and assumption."""
    parameters,edges,observed,errors=load_inputs()
    results={'c4_first_six_statistical':fit_c4(parameters,edges,observed,errors)}
    results['c4_first_six_common_systematics']=fit_c4(parameters,edges,observed,errors,(.15,.051))
    results['c4_first_six_diagonal_systematics']=fit_c4(parameters,edges,observed,errors,(.15,.051),True)
    for name,systematics,diagonal in (('statistical',(),False),('common_systematics',(.15,.051),False),('diagonal_systematics',(.15,.051),True)):
        results['delta_complex_rho_prime_'+name]=fit_delta(parameters,edges,observed,errors,relative_systematics=systematics,diagonal=diagonal)
    results['delta_real_rho_prime_statistical']=fit_delta(parameters,edges,observed,errors,real_residue=True)
    results['delta_unsubtracted_rho_prime_statistical']=fit_delta(parameters,edges,observed,errors,subtraction=False)
    results['joint_c4_delta_complex_rho_prime_statistical']=fit_joint_coefficients(parameters,edges,observed,errors)
    extended_path=ROOT/'data/vector_sector/babar_2007_isovector_kstar_k_below_2gev.csv'
    if extended_path.exists():
        extended=np.genfromtxt(extended_path,delimiter=',',names=True)
        extended_edges=np.column_stack((extended['energy_low_GeV'],extended['energy_high_GeV']))
        results['delta_complex_rho_prime_below_2gev_statistical']=fit_delta(parameters,extended_edges,extended['cross_section_nb'],extended['statistical_error_nb'])
        results['joint_c4_delta_complex_rho_prime_below_2gev_statistical']=fit_joint_coefficients(parameters,extended_edges,extended['cross_section_nb'],extended['statistical_error_nb'])
    return dict(bibtex_keys=['BaBar:2007ceh','Alda:2026nrd','Ruiz-Femenia:2003jdx','GomezDumm:2010ra','ParticleDataGroup:2024cfk'],
        scope='Actual RChT coefficient fits to the BaBar isovector K*(892)K component',
        conventions='epsilon^0123=+1; canonical c_i,d_i=-reference c_i,d_i; current scalar per canonical coefficient follows the existing linear gradients',
        normalization_derivation=dict(current='H^mu=i epsilon^{mu nu alpha beta} eps*_nu pV_alpha pK_beta F_VP',
            spin_sum='sum H_i H_j*=s |F_VP|² (p² delta_ij-p_i p_j)',
            single_charge_cross_section='pi alpha² lambda(s,MKstar²,MK²)^(3/2) |F_VP|²/(6s³)',
            charge_conjugate_factor=2,
            isospin_component='BaBar Eq.13: sigma_charged_pair=|sqrt(sigma0)-sqrt(sigma1) exp(i phase)|², so pure isovector charged pair plus c.c. equals sigma1',
            source_traces='FQ(QI)=1/2; FQ(T_rho)=1/sqrt2; R_rho=F_l/sqrt2, F_l=F_V+2z_V chi_plus_l'),
        inputs=dict(preliminary_parameters=parameters,F_GeV=F,operator_mass_GeV=OPERATOR_MASS,kaon_phase_space_mass_GeV=KAON_MASS,pion_mass_GeV=PION_MASS,chi_plus_light_GeV2=CHI_PLUS_LIGHT,chi_plus_strange_GeV2=CHI_PLUS_STRANGE),
        assumptions=['External stable K* and quasi-two-body VP phase space; no inclusive KKpi substitution',
            'BaBar sigma1 is one isovector charge pair plus c.c.; cross section pi alpha² lambda³/² |FVP|²/(3s³)',
            'Remaining original coefficients and even-source inputs fixed at preliminary central representative; no OPE coefficient constraints',
            'delta c1=delta c5=delta/2 preserves real-photon VJP numerators exactly',
            'Uniform dE bin averaging; diagonal statistical errors because covariance is unpublished',
            'rho P-wave pi-pi width and rho-prime BaBar 4pi width with fixed B4pi=.64',
            'Excited pole subtraction Dprime(s)-Dprime(0) matches its contribution to zero at the real-photon point; explicit matching assumption',
            'Unsubtracted variant does not preserve preliminary radiative observations and is a model sensitivity only',
            '15% Dalitz and 5.1% Table1 systematics correlations are explicit common/diagonal assumptions; diagonal case is Pearson chi2 with predicted multiplicative variance',
            'Numerical search bounds are not theory priors; profile scans retain their finite endpoint qualification'],
        exact_invisible_directions=['c7','d4','common c2,c5,c6 shift','all g1,...,g5','other directions not fixed by the single isovector VP form factor remain unmeasured'],
        coefficient_coordinate=dict(initial_c1_minus_c2_plus_c5=parameters['c1']-parameters['c2']+parameters['c5'],changes={'c1':'delta/2','c5':'delta/2'}),
        bin_edges_GeV=edges.tolist(),observed_nb=observed.tolist(),statistical_errors_nb=errors.tolist(),fits=results)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    arguments=parser.parse_args()
    arguments.output.write_text(json.dumps(run_fits(),indent=2,allow_nan=False)+'\n')
    print(arguments.output)
