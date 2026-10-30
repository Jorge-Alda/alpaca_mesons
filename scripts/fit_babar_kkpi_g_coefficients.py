"""Fit local VPPP coefficients to inclusive BaBar:2007ceh KKpi spectra.

Run: .venv/bin/python scripts/fit_babar_kkpi_g_coefficients.py --output data/vector_sector/babar_kkpi_g_fit.json
The vector pieces of Qin:2024ulb Eqs. (kpkmpi0),(kskpi), Appendix C are
used with no short-distance elimination. Natural units, GeV, canonical
epsilon^0123=+1; canonical odd coefficients are minus source coefficients.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2

try:
    from implementation.kkpi_vector_current import energy_bin_matrices
except ModuleNotFoundError:
    from scripts.implementation.kkpi_vector_current import energy_bin_matrices

ROOT = Path(__file__).resolve().parents[1]
F = .093
MV = .770
ALPHA = 1/137.035999084
MK = .495644
MPI = .1372745
COEFFICIENT_NAMES = ('g1_plus_2g2_minus_g3', 'g2', 'g4', '2g4_plus_g5')
CHANNELS = ('kkpi0', 'kskpi')


def load_inputs():
    parameters = json.loads((ROOT/'data/vector_sector/odd_preliminary_fit.json').read_text())['parameters']
    data = {}
    for channel in CHANNELS:
        path = ROOT/f'data/vector_sector/babar_2007_{channel}_inclusive_below_1p78.csv'
        table = np.genfromtxt(path, delimiter=',', names=True)
        data[channel] = dict(edges=np.column_stack((table['energy_low_GeV'], table['energy_high_GeV'])),
            observed=table['cross_section_nb'], errors=table['statistical_error_nb'],
            sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return parameters, data


def triangle(s, first, second):
    return (s-first-second)**2-4*first*second


def propagator(name, invariant, parameters, beta=0.):
    """Ground state plus one excited multiplet; independent PDG pole inputs.

    The same beta modifies all vector species in a declared channel, as in
    Qin:2024ulb Eq.(bw). Ground rho/K* widths use the source Appendix D.
    Excited masses/widths are the independent PDG figures printed there,
    not the paper's fitted pole values.
    """
    invariant = np.asarray(invariant)
    beta_first,beta_second = beta if np.ndim(beta) else (beta,0.)
    if name == 'rho':
        mass = parameters['m_rho_GeV']
        pion_phase = np.maximum(1-4*MPI**2/np.maximum(invariant, 1e-20), 0)**1.5
        kaon_phase = np.maximum(1-4*MK**2/np.maximum(invariant, 1e-20), 0)**1.5
        width = mass*invariant/(96*np.pi*F**2)*(pion_phase+.5*kaon_phase)
        excited_mass = parameters.get('rho_prime_mass_GeV',1.465)
        excited_width = parameters.get('rho_prime_width_GeV',.400)
        ratio = pion_phase/(1-4*MPI**2/excited_mass**2)**1.5
        excited_running_width = excited_width*np.sqrt(np.maximum(invariant,0))/excited_mass*ratio
        second_mass,second_width=1.720,.250
        second_running_width=second_width*np.sqrt(np.maximum(invariant,0))/second_mass*pion_phase/(1-4*MPI**2/second_mass**2)**1.5
    elif name == 'kstar':
        mass = parameters['m_Kstar_GeV']
        def phase(second_mass):
            return np.where(invariant > (MK+second_mass)**2,
                np.maximum(triangle(invariant,MK**2,second_mass**2),0)**1.5, 0)
        width = mass*(phase(MPI)+phase(.547862))/(128*np.pi*F**2*np.maximum(invariant,1e-20)**2)
        excited_mass, excited_width = 1.414, .232
        ratio = np.maximum(triangle(invariant,MK**2,MPI**2),0)/triangle(excited_mass**2,MK**2,MPI**2)
        excited_running_width = np.where(invariant > (MK+MPI)**2,
            excited_width*excited_mass**4/np.maximum(invariant,1e-20)**2*ratio**1.5,0)
        second_mass,second_width=1.718,.322
        second_ratio=np.maximum(triangle(invariant,MK**2,MPI**2),0)/triangle(second_mass**2,MK**2,MPI**2)
        second_running_width=np.where(invariant>(MK+MPI)**2,
            second_width*second_mass**4/np.maximum(invariant,1e-20)**2*second_ratio**1.5,0)
    else:
        mass = {'omega':.78266, 'phi':1.019461}[name]
        width = parameters[f'Gamma_total_{name}_GeV']
        excited_mass, excited_width = {'omega':(1.410,.290),'phi':(1.680,.150)}[name]
        if name=='phi':
            excited_mass=parameters.get('phi_prime_mass_GeV',excited_mass)
            excited_width=parameters.get('phi_prime_width_GeV',excited_width)
        else:
            excited_mass=parameters.get('omega_prime_mass_GeV',excited_mass)
            excited_width=parameters.get('omega_prime_width_GeV',excited_width)
        excited_running_width = excited_width
        second_mass,second_width={'omega':(1.670,.315),'phi':(2.162,.100)}[name]
        second_running_width=second_width
    return (1/(mass**2-invariant-1j*mass*width)
            +beta_first/(excited_mass**2-invariant-1j*excited_mass*excited_running_width)
            +beta_second/(second_mass**2-invariant-1j*second_mass*second_running_width))


def form_factor(q_squared, pair_s, pair_t, parameters, channel, g_coordinates, beta=0., source_breaking=True):
    """Full coherent vector form factor; tensor terms explicitly omitted.

    Variables follow Qin:2024ulb: s=(pK+pKbar)^2, t=(pKbar+ppi)^2.
    The source's DR(q²,mK²,t) in Fd1 is read as DR(q²,t,mK²),
    consistently with its definition and the t/u pair exchange.
    """
    pair_s, pair_t = np.broadcast_arrays(pair_s, pair_t)
    pair_u = q_squared+2*MK**2+MPI**2-pair_s-pair_t
    derivative_combination, g2, g4, mass_combination = g_coordinates
    vector_source = parameters['F_V_GeV']
    vector_pion = parameters['G_V_GeV']
    c = [0]+[parameters[f'c{i}'] for i in range(1,8)]
    d = [0]+[parameters[f'd{i}'] for i in range(1,5)]
    def cr1(pair, mass_squared):
        return ((c[1]-c[2]+c[5])*q_squared
            -(c[1]-c[2]-c[5]+2*c[6])*pair
            +(c[1]+c[2]+8*c[3]-c[5])*mass_squared)
    def cr2(pair):
        return cr1(pair,MK**2)+24*c[4]*(MK**2-MPI**2)
    def dr(pair, mass_squared):
        return (d[1]+8*d[2]-d[3])*mass_squared+d[3]*(q_squared+pair)
    gr1 = (derivative_combination*(q_squared+pair_s-4*MK**2-MPI**2)
        +4*mass_combination*MK**2-4*g2*(q_squared-2*MK**2-MPI**2)
        -4*g4*(MK**2-MPI**2))
    gr2 = (derivative_combination*(q_squared-pair_s-MPI**2)
        +2*mass_combination*MPI**2-2*g2*(q_squared-2*MK**2-MPI**2)
        +4*g4*(MK**2-MPI**2))
    rho_q, omega_q, phi_q = (propagator(name,q_squared,parameters,beta) for name in ('rho','omega','phi'))
    rho_s = propagator('rho',pair_s,parameters,beta)
    omega_s = propagator('omega',pair_s,parameters,beta)
    kstar_t = propagator('kstar',pair_t,parameters,beta)
    kstar_u = propagator('kstar',pair_u,parameters,beta)
    light_source = vector_source+2*parameters['z_V_GeV_inverse']*.036616 if source_breaking else vector_source
    strange_source = vector_source+2*parameters['z_V_GeV_inverse']*.943684 if source_breaking else vector_source
    rho_source,omega_source,phi_source = light_source*rho_q,light_source*omega_q,strange_source*phi_q
    contact_prefactor = 2*np.sqrt(2)/(F**3*MV)
    direct_prefactor = 2*np.sqrt(2)*vector_pion/(3*F**3*MV)
    exchange_prefactor = 2*vector_pion/F**3
    if channel == 'kkpi0':
        # The canonical WZW offset reverses the reference scalar sign.
        anomaly = np.full(pair_s.shape,3/(12*np.pi**2*F**3),complex)
        local = contact_prefactor*((rho_source+omega_source/3)*gr1-2*phi_source*gr2/3)
        direct = -direct_prefactor*((rho_s+3*omega_s)*cr1(pair_s,MPI**2)
            +cr2(pair_t)*kstar_t+cr2(pair_u)*kstar_u)
        exchange = exchange_prefactor*(2*dr(pair_s,MPI**2)*(rho_source*omega_s+omega_source*rho_s/3)
            +(rho_source+omega_source/3-2*phi_source/3)*(dr(pair_t,MK**2)*kstar_t+dr(pair_u,MK**2)*kstar_u))
    elif channel == 'kskpi':
        anomaly = np.zeros(pair_s.shape,complex)
        local = contact_prefactor*(derivative_combination*(pair_t-pair_u)*rho_source-omega_source*gr1/3+2*phi_source*gr2/3)
        direct = direct_prefactor*(rho_s*cr1(pair_s,MPI**2)
            -2*cr1(pair_t,MK**2)*kstar_t+cr2(pair_u)*kstar_u)
        exchange = -exchange_prefactor*(2*omega_source*rho_s*dr(pair_s,MPI**2)/3
            -(rho_source-omega_source/3+2*phi_source/3)*dr(pair_t,MK**2)*kstar_t
            +(rho_source+omega_source/3-2*phi_source/3)*dr(pair_u,MK**2)*kstar_u)
    else:
        raise ValueError(channel)
    return anomaly+local+direct+exchange


def polynomial_basis(q_squared, pair_s, pair_t, parameters, channel, source_breaking=True):
    """Eleven columns represent beta² and the four g×beta products exactly."""
    zero = np.zeros(4)
    fixed_zero = form_factor(q_squared,pair_s,pair_t,parameters,channel,zero,0,source_breaking)
    fixed_plus = form_factor(q_squared,pair_s,pair_t,parameters,channel,zero,1,source_breaking)
    fixed_minus = form_factor(q_squared,pair_s,pair_t,parameters,channel,zero,-1,source_breaking)
    fixed_linear = (fixed_plus-fixed_minus)/2
    fixed_quadratic = (fixed_plus+fixed_minus)/2-fixed_zero
    columns = [fixed_zero,fixed_linear,fixed_quadratic]
    zero_gradients, beta_gradients = [], []
    for index in range(4):
        unit = np.eye(4)[index]
        at_zero = form_factor(q_squared,pair_s,pair_t,parameters,channel,unit,0,source_breaking)-fixed_zero
        at_plus = form_factor(q_squared,pair_s,pair_t,parameters,channel,unit,1,source_breaking)-fixed_plus
        zero_gradients.append(at_zero)
        beta_gradients.append(at_plus-at_zero)
    return np.stack(columns+zero_gradients+beta_gradients,axis=-1)


def prediction_coefficients(g_coordinates,beta):
    if np.ndim(beta):
        first,second=beta
        return np.r_[1,first,second,first**2,first*second,second**2,
            g_coordinates,first*np.asarray(g_coordinates),second*np.asarray(g_coordinates)]
    return np.r_[1,beta,beta**2,g_coordinates,beta*np.asarray(g_coordinates)]


def control_parameters(parameters):
    """Qin:2024ulb constraints converted into canonical coefficient signs."""
    result=dict(parameters)
    result['G_V_GeV']=F**2/result['F_V_GeV']
    for index in range(1,8):
        result[f'c{index}']=0.
    for index in range(1,5):
        result[f'd{index}']=0.
    result['c6']=3*MV/(192*np.sqrt(2)*np.pi**2*result['G_V_GeV'])
    result['d3']=3*MV**2/(192*np.pi**2*result['F_V_GeV']*result['G_V_GeV'])
    result['d1']=result['d3']-F**2/(8*result['F_V_GeV']**2)
    result['z_V_GeV_inverse']=0.
    return result


def control_polynomial_basis(q_squared,pair_s,pair_t,parameters,channel,multiplets=1,a2=False):
    fixed_g2=-3*MV/(192*np.sqrt(2)*np.pi**2*parameters['F_V_GeV'])
    coefficient_count=5 if a2 else 3
    def evaluate(coordinates,beta):
        varied=dict(parameters,c4=coordinates[2])
        result=form_factor(q_squared,pair_s,pair_t,varied,channel,
            [0.,fixed_g2,coordinates[0],coordinates[1]],beta,False)
        if a2 and channel=='kskpi':
            pair_u=q_squared+2*MK**2+MPI**2-pair_s-pair_t
            tensor_response=1/(1.3182**2-pair_s-1j*1.3182*.107)
            g_tensor=.0233
            result += (-2*g_tensor/F**3*coordinates[3]*(pair_t-pair_u)*tensor_response
                +np.sqrt(2)*parameters['F_V_GeV']*g_tensor/F**3*coordinates[4]
                *(pair_t-pair_u)*tensor_response*propagator('rho',q_squared,parameters,beta))
        return result
    if multiplets==2:
        fixed_zero=evaluate(np.zeros(coefficient_count),(0,0))
        first_plus=evaluate(np.zeros(coefficient_count),(1,0));first_minus=evaluate(np.zeros(coefficient_count),(-1,0))
        second_plus=evaluate(np.zeros(coefficient_count),(0,1));second_minus=evaluate(np.zeros(coefficient_count),(0,-1))
        linear_first=(first_plus-first_minus)/2
        linear_second=(second_plus-second_minus)/2
        square_first=(first_plus+first_minus)/2-fixed_zero
        square_second=(second_plus+second_minus)/2-fixed_zero
        mixed=evaluate(np.zeros(coefficient_count),(1,1))-fixed_zero-linear_first-linear_second-square_first-square_second
        columns=[fixed_zero,linear_first,linear_second,square_first,mixed,square_second]
        gradients_zero=[];gradients_first=[];gradients_second=[]
        for unit in np.eye(coefficient_count):
            zero=evaluate(unit,(0,0))-fixed_zero
            gradients_zero.append(zero)
            gradients_first.append(evaluate(unit,(1,0))-first_plus-zero)
            gradients_second.append(evaluate(unit,(0,1))-second_plus-zero)
        return np.stack(columns+gradients_zero+gradients_first+gradients_second,axis=-1)
    fixed_zero=evaluate(np.zeros(3),0)
    fixed_plus=evaluate(np.zeros(3),1)
    fixed_minus=evaluate(np.zeros(3),-1)
    columns=[fixed_zero,(fixed_plus-fixed_minus)/2,(fixed_plus+fixed_minus)/2-fixed_zero]
    zero_gradients=[];beta_gradients=[]
    for unit in np.eye(3):
        zero=evaluate(unit,0)-fixed_zero
        plus=evaluate(unit,1)-fixed_plus
        zero_gradients.append(zero);beta_gradients.append(plus-zero)
    return np.stack(columns+zero_gradients+beta_gradients,axis=-1)


def prepare_matrices(parameters,data,energy_order=6,pair_order=48,angle_order=48,source_breaking=True,control=False,multiplets=1,cutoff=None,a2=False):
    matrices = {}
    for channel in CHANNELS:
        masses = ((.493677,.493677,.1349768) if channel == 'kkpi0'
            else (.497611,.493677,.13957039))
        def basis(q,s,pair_s13):
            # The integration helper supplies s13; the source calls s23 t.
            pair_s23=q+sum(mass**2 for mass in masses)-s-pair_s13
            if control:
                scalar=control_polynomial_basis(q,s,pair_s23,parameters,channel,multiplets,a2)
                return scalar*np.exp(-q**3/cutoff**6) if cutoff else scalar
            return polynomial_basis(q,s,pair_s23,parameters,channel,source_breaking)
        matrices[channel] = energy_bin_matrices(data[channel]['edges'],masses,basis,ALPHA,
            energy_order=energy_order,pair_order=pair_order,angle_order=angle_order)
        if channel == 'kskpi':
            matrices[channel] *= 2
    return matrices


def predict_channel(values,matrices,channel,excited,coefficient_count=4,multiplets=1):
    channel_index=CHANNELS.index(channel)
    beta = values[coefficient_count+channel_index] if excited else 0.
    if multiplets==2:
        beta=(beta,values[coefficient_count+2+channel_index] if excited else 0.)
    coefficients = prediction_coefficients(values[:coefficient_count],beta)
    return np.einsum('i,bij,j->b',coefficients,matrices[channel],coefficients)


def residual_vector(values,matrices,data,excited,systematics=True,coefficient_count=4,multiplets=1):
    residuals = []
    for channel,relative in zip(CHANNELS,(.065,.051)):
        prediction = predict_channel(values,matrices,channel,excited,coefficient_count,multiplets)
        residual = (prediction-data[channel]['observed'])/data[channel]['errors']
        if systematics:
            column = relative*prediction/data[channel]['errors']
            nuisance = -column@residual/(1+column@column)
            residuals.extend(np.r_[residual+column*nuisance,nuisance])
        else:
            residuals.extend(residual)
    return np.asarray(residuals)


def fit_coefficients(matrices,data,excited,systematics=True,profiles=True,coefficient_names=COEFFICIENT_NAMES,multiplets=1,random_starts=None,extra_initial=None):
    coefficient_count=len(coefficient_names)
    number = coefficient_count+2*multiplets if excited else coefficient_count
    lower = np.r_[np.full(coefficient_count,-3.),np.full(number-coefficient_count,-2.)]
    upper = -lower
    random = np.random.default_rng(4818)
    initial_points = [np.zeros(number)]
    if extra_initial is not None:
        initial_points.append(np.asarray(extra_initial))
    if coefficient_names[0]=='g4':
        beta_start=[-.2685,-.4271]+([-.0289,.0377] if multiplets==2 else [])
        coefficient_start=[.047,.0294,.0004]+([-.0342,-.7022] if coefficient_count==5 else [])
        initial_points.append(np.r_[coefficient_start,beta_start if excited else []])
    if coefficient_names[0]=='g4' and multiplets==2 and excited:
        initial_points.append(np.r_[[.17695,.37691,.0068]+([0.,0.] if coefficient_count==5 else []),
            [.10736,-.33721,-.03463,-.27871]])
    for index in range((23 if multiplets==2 else 11) if random_starts is None else random_starts):
        initial = random.normal(0,.08,number)
        if coefficient_names[0]=='g4':
            initial[2] *= .05
            if coefficient_count==5:
                initial[4]=random.uniform(-1.,1.)
        if excited:
            initial[coefficient_count:] = random.uniform(-.8,.2,2*multiplets)
        initial_points.append(initial)
    function = lambda values:residual_vector(values,matrices,data,excited,systematics,coefficient_count,multiplets)
    candidates = [least_squares(function,initial,bounds=(lower,upper),x_scale='jac',
        max_nfev=1500,ftol=1e-10,xtol=1e-10,gtol=1e-9) for initial in initial_points]
    best = min(candidates,key=lambda result:result.fun@result.fun)
    statistic = float(best.fun@best.fun)
    covariance = np.linalg.pinv(best.jac.T@best.jac,rcond=1e-13)
    jacobian_g = best.jac[:,:coefficient_count]
    jacobian_nuisance = best.jac[:,coefficient_count:]
    if jacobian_nuisance.shape[1]:
        jacobian_g = jacobian_g-jacobian_nuisance@np.linalg.pinv(jacobian_nuisance)@jacobian_g
    singular_values = np.linalg.svd(jacobian_g,compute_uv=False)
    coefficient_rank = int(sum(singular_values > singular_values[0]*1e-7))
    beta_names=[f'beta_excited_{channel}' for channel in CHANNELS]
    if multiplets==2:
        beta_names += [f'beta_second_excited_{channel}' for channel in CHANNELS]
    names = list(coefficient_names)+(beta_names if excited else [])
    output = dict(parameters=dict(zip(names,map(float,best.x))),chi_squared=statistic,
        degrees_of_freedom=sum(len(data[c]['observed']) for c in CHANNELS)-number,
        coefficient_rank_after_nuisance_projection=coefficient_rank,
        coefficient_singular_values=singular_values.tolist(),local_covariance=covariance.tolist(),
        local_errors=dict(zip(names,map(float,np.sqrt(np.diag(covariance))))),
        optimizer_success=bool(best.success),optimizer_message=best.message,
        active_bounds=[names[i] for i in range(number) if min(best.x[i]-lower[i],upper[i]-best.x[i])<1e-4],
        systematics='separate common 6.5% and 5.1% Gaussian normalization nuisances' if systematics else 'statistical only',
        number_of_starts=len(candidates),global_minimum_proven=False)
    output['probability']=float(chi2.sf(statistic,output['degrees_of_freedom']))
    if coefficient_names[0]=='g4':
        local_vppp_information=np.linalg.pinv(covariance[:2,:2],rcond=1e-13)
        vppp_values=np.sqrt(np.maximum(np.linalg.eigvalsh(local_vppp_information),0))[::-1]
        output['vppp_rank_after_c4_and_beta_projection']=int(sum(vppp_values>vppp_values[0]*1e-7))
        output['vppp_rank_after_all_other_parameters']=int(sum(vppp_values>vppp_values[0]*1e-7))
        output['vppp_profiled_nuisances']=names[2:]
        output['vppp_information_singular_values']=vppp_values.tolist()
    if coefficient_count==4:
        output['canonical_representative'] = dict(g1=float(best.x[0]-2*best.x[1]),g2=float(best.x[1]),
            g3=0.,g4=float(best.x[2]),g5=float(best.x[3]-2*best.x[2]))
    output['predictions_nb']={channel:predict_channel(best.x,matrices,channel,excited,coefficient_count,multiplets).tolist() for channel in CHANNELS}
    output['channel_statistical_chi_squared']={channel:float(np.sum(((np.asarray(output['predictions_nb'][channel])-data[channel]['observed'])/data[channel]['errors'])**2)) for channel in CHANNELS}
    output['normalization_gaussian_pulls']={}
    output['normalization_multipliers']={}
    output['profiled_predictions_nb']={}
    for channel,relative in zip(CHANNELS,(.065,.051)):
        prediction=np.asarray(output['predictions_nb'][channel])
        column=relative*prediction/data[channel]['errors']
        residual=(prediction-data[channel]['observed'])/data[channel]['errors']
        nuisance=float(-column@residual/(1+column@column)) if systematics else 0.
        output['normalization_gaussian_pulls'][channel]=nuisance
        output['normalization_multipliers'][channel]=1+relative*nuisance
        output['profiled_predictions_nb'][channel]=(prediction*(1+relative*nuisance)).tolist()
    output['local_errors_qualification']='Local linear diagnostic; no physical calibration if residuals are inadequate, bounds active, or minimization incomplete'
    if profiles:
        output['profiles']={}
        profile_names=coefficient_names[:2] if coefficient_count==5 else coefficient_names
        for index,name in enumerate(profile_names):
            error = np.sqrt(covariance[index,index])
            half_span = min(2.,max(.002 if name=='c4' else .03,5*error))
            grid = np.unique(np.r_[np.linspace(max(lower[index],best.x[index]-half_span),
                min(upper[index],best.x[index]+half_span),13),best.x[index]])
            active = [i for i in range(number) if i!=index]
            values=[];convergence=[]
            continuation=best.x[active]
            for fixed in grid:
                def profile_residual(free):
                    full=best.x.copy();full[active]=free;full[index]=fixed
                    return function(full)
                fits=[least_squares(profile_residual,initial,bounds=(lower[active],upper[active]),
                    x_scale='jac',max_nfev=350,ftol=1e-8,xtol=1e-8,gtol=1e-8)
                    for initial in (best.x[active],)]
                result=min(fits,key=lambda item:item.fun@item.fun)
                continuation=result.x
                values.append(float(result.fun@result.fun))
                convergence.append(bool(result.success))
            crossings={}
            for level in (1.,3.841458820694124):
                difference=np.asarray(values)-statistic-level
                crossings[str(level)]=[float(grid[j]+(grid[j+1]-grid[j])*(-difference[j])/(difference[j+1]-difference[j]))
                    for j in range(len(grid)-1) if difference[j]*difference[j+1]<0]
            output['profiles'][name]=dict(grid=grid.tolist(),chi_squared=values,crossings=crossings,
                optimizer_success=convergence,endpoints_are_confidence_limits=False,
                search_qualification='One multistart optimum initializes each bounded profile; branch completeness is not proven')
        profile_minimum=min(min(item['chi_squared']) for item in output['profiles'].values())
        output['minimum_profile_chi_squared_minus_reported_minimum']=float(profile_minimum-statistic)
        output['profile_intervals_are_final']=bool(profile_minimum>=statistic-1e-4
            and all(all(item['optimizer_success']) for item in output['profiles'].values())
            and output['probability']>.01)
    return output


def run_analysis(source_breaking=True):
    parameters,data=load_inputs()
    matrices=prepare_matrices(parameters,data,source_breaking=source_breaking)
    ground=fit_coefficients(matrices,data,False,profiles=False)
    excited=fit_coefficients(matrices,data,True,profiles=False)
    statistical=fit_coefficients(matrices,data,True,systematics=False,profiles=False)
    constrained=control_parameters(parameters)
    control_matrices=prepare_matrices(constrained,data,source_breaking=False,control=True)
    names=('g4','2g4_plus_g5','c4')
    control_ground=fit_coefficients(control_matrices,data,False,profiles=False,coefficient_names=names)
    control_excited=fit_coefficients(control_matrices,data,True,coefficient_names=names)
    control_statistical=fit_coefficients(control_matrices,data,True,systematics=False,profiles=False,coefficient_names=names)
    refined=prepare_matrices(parameters,data,energy_order=8,pair_order=72,angle_order=72,source_breaking=source_breaking)
    best=np.array(list(excited['parameters'].values()))
    changes={channel:float(np.max(np.abs(predict_channel(best,refined,channel,True)
        -predict_channel(best,matrices,channel,True)))) for channel in CHANNELS}
    control_refined=prepare_matrices(constrained,data,energy_order=8,pair_order=72,angle_order=72,source_breaking=False,control=True)
    control_best=np.array(list(control_excited['parameters'].values()))
    control_changes={channel:float(np.max(np.abs(predict_channel(control_best,control_refined,channel,True,3)
        -predict_channel(control_best,control_matrices,channel,True,3)))) for channel in CHANNELS}
    return dict(bibtex_keys=['BaBar:2007ceh','Qin:2024ulb','Alda:2026nrd','Ruiz-Femenia:2003jdx'],
        scope='Full inclusive electromagnetic KKpi spectra with local VPPP coefficients fitted',
        conventions='canonical epsilon^0123=+1; c_i,d_i,g_i=-reference; WZW scalar reversed and linear gradients retained',
        inputs=dict(preliminary_parameters=parameters,F_GeV=F,operator_mass_GeV=MV,vertex_kaon_mass_GeV=MK,vertex_pion_mass_GeV=MPI,
            source_breaking=source_breaking,source_light_GeV=parameters['F_V_GeV']+2*parameters['z_V_GeV_inverse']*.036616 if source_breaking else parameters['F_V_GeV'],
            source_strange_GeV=parameters['F_V_GeV']+2*parameters['z_V_GeV_inverse']*.943684 if source_breaking else parameters['F_V_GeV']),
        data={channel:{key:value.tolist() if isinstance(value,np.ndarray) else value for key,value in block.items()} for channel,block in data.items()},
        assumptions=['Vector pieces of Qin:2024ulb only; tensor contributions omitted in the below-1.78-GeV window',
            'Ideal neutral-vector mixing; all virtual-current rho/omega poles use F_l and phi poles use F_s from the preliminary source point' if source_breaking else 'Ideal mixing and universal FV source comparison',
            'All c_i,d_i,FV,GV fixed to the existing preliminary central point; no g OPE constraint',
            'Physical phase-space masses; averaged kaon/pion masses inside source form factors',
            'First excited vector multiplet with PDG masses/widths fixed, common beta across vector species within each channel',
            'Both charge conjugates included for KS K±pi∓; inclusive Table2/3 only, no Table7 double counting',
            'Uniform dE bin averages; statistical covariance unpublished; separate common normalization systematics'],
        exact_null='g1 and g3 may be shifted together; the displayed g3=0 is only a representative',
        fits=dict(ground_vector_only=ground,excited_vector_profiled=excited,excited_vector_statistical=statistical),
        source_constrained_control=dict(parameters=constrained,
            assumptions=['Qin:2024ulb high-energy relations are fixed model assumptions, not measured likelihood rows',
                'g1+2g2-g3=0; canonical g2=-Nc MV/(192sqrt2 pi² FV)',
                'c1-c2+c5=0; c1+c2+8c3-c5=0; canonical c1-c2-c5+2c6=+Nc MV/(96sqrt2 pi² GV)',
                'canonical d3=+Nc MV²/(192pi² FV GV); canonical d1+8d2-d3=-F²/(8FV²)',
                'FV GV=F²; FV fixed at the preliminary central value; universal source as in the source model'],
            fits=dict(ground_vector_only=control_ground,excited_vector_profiled=control_excited,excited_vector_statistical=control_statistical),
            quadrature_refinement_max_absolute_change_nb=control_changes),
        quadrature_refinement_max_absolute_change_nb=changes)


def run_two_multiplet_analysis(a2=False):
    preliminary,data=load_inputs()
    parameters=control_parameters(preliminary)
    fits={};matrices_by_cutoff={}
    for cutoff in (2.5,2.3):
        matrices=prepare_matrices(parameters,data,source_breaking=False,control=True,multiplets=2,cutoff=cutoff,a2=a2)
        matrices_by_cutoff[cutoff]=matrices
        names=('g4','2g4_plus_g5','c4')+(('C_T_GeV_inverse','D_T') if a2 else ())
        fits[str(cutoff)]=fit_coefficients(matrices,data,True,
            coefficient_names=names,multiplets=2,profiles=cutoff==2.5 and not a2)
        print('Two multiplets, cutoff',cutoff,'chi2',fits[str(cutoff)]['chi_squared'],
            fits[str(cutoff)]['parameters'],flush=True)
    refined=prepare_matrices(parameters,data,energy_order=8,pair_order=72,angle_order=72,
        source_breaking=False,control=True,multiplets=2,cutoff=2.5,a2=a2)
    best=np.array(list(fits['2.5']['parameters'].values()))
    coefficient_count=5 if a2 else 3
    changes={channel:float(np.max(np.abs(predict_channel(best,refined,channel,True,coefficient_count,2)
        -predict_channel(best,matrices_by_cutoff[2.5],channel,True,coefficient_count,2)))) for channel in CHANNELS}
    fixed_g2=-3*MV/(192*np.sqrt(2)*np.pi**2*parameters['F_V_GeV'])
    for fit in fits.values():
        fit['canonical_g_representative']=dict(g1=-2*fixed_g2,g2=fixed_g2,g3=0.,
            g4=fit['parameters']['g4'],g5=fit['parameters']['2g4_plus_g5']-2*fit['parameters']['g4'])
        if a2:
            fit['tensor_coefficient_products']=dict(gT_times_CT_dimensionless=.0233*fit['parameters']['C_T_GeV_inverse'],
                gT_times_DT_GeV=.0233*fit['parameters']['D_T'])
    return dict(bibtex_keys=['BaBar:2007ceh','Qin:2024ulb','Alda:2026nrd'],
        scope='Inclusive BaBar source-constrained VPPP fit with two excited vector multiplets'+(' and a2 exchange' if a2 else ''),
        conventions='Canonical odd coefficients=-reference; canonical WZW scalar is positive',
        inputs=dict(preliminary_FV_GeV=parameters['F_V_GeV'],constrained_parameters=parameters,
            F_GeV=F,operator_mass_GeV=MV,fixed_canonical_g2=fixed_g2,
            gT_GeV=.0233 if a2 else None,a2_mass_GeV=1.3182 if a2 else None,a2_width_GeV=.107 if a2 else None),
        constraints_as_model_assumptions=['g1+2g2-g3=0','FV GV=F²',
            'g2=-Nc MV/(192sqrt2 pi² FV)',
            'c1-c2+c5=0','c1+c2+8c3-c5=0',
            'c1-c2-c5+2c6=+Nc MV/(96sqrt2 pi² GV)',
            'd3=+Nc MV²/(192pi² FV GV)','d1+8d2-d3=-F²/(8FV²)'],
        assumptions=['Qin WZW+VPPP+VJP+VVP and source a2 terms in KS channel; K2* tails omitted' if a2 else 'Qin vector terms WZW+VPPP+VJP+VVP; tensor terms omitted below1.78GeV',
            'First and second excited multiplet residues are independent real nuisance parameters in each channel',
            'Independent fixed PDG resonance masses/widths; no literature fitted g values are likelihood inputs',
            'Universal source and ideal mixing as the source model; fixed FV from preliminary central point',
            'CT has GeV^-1 by operator dimension; the source parameter table unit label GeV is not used',
            'Regulator exp[-(q²)^3/Lambda^6] applied to full amplitude',
            'Inclusive Tables2/3 only; 6.5%/5.1% common normalization nuisances; no Dalitz15% added',
            'Uniform energy-bin average, physical phase-space masses, averaged vertex masses'],
        data={channel:{key:value.tolist() if isinstance(value,np.ndarray) else value for key,value in block.items()} for channel,block in data.items()},
        fits=fits,quadrature_refinement_max_absolute_change_nb=changes,
        exact_null='delta g1=delta g3 remains invisible; constraints are assumptions, not measurements')


def run_phi_prime_scan():
    preliminary,data=load_inputs()
    parameters=control_parameters(preliminary)
    fits=[];continuation=None
    names=('g4','2g4_plus_g5','c4','C_T_GeV_inverse','D_T')
    for mass in (1.680,1.640,1.710):
        for width in (.150,.240,.325):
            varied=dict(parameters,phi_prime_mass_GeV=mass,phi_prime_width_GeV=width)
            matrices=prepare_matrices(varied,data,source_breaking=False,control=True,
                multiplets=2,cutoff=2.5,a2=True)
            fit=fit_coefficients(matrices,data,True,coefficient_names=names,
                multiplets=2,profiles=False,random_starts=5,extra_initial=continuation)
            continuation=np.array(list(fit['parameters'].values()))
            fits.append(dict(phi_prime_mass_GeV=mass,phi_prime_width_GeV=width,fit=fit))
            print('Phi-prime grid',mass,width,'chi2',fit['chi_squared'],flush=True)
    selected=min(fits,key=lambda item:item['fit']['chi_squared'])
    varied=dict(parameters,phi_prime_mass_GeV=selected['phi_prime_mass_GeV'],
        phi_prime_width_GeV=selected['phi_prime_width_GeV'])
    matrices=prepare_matrices(varied,data,source_breaking=False,control=True,multiplets=2,cutoff=2.5,a2=True)
    preferred=fit_coefficients(matrices,data,True,coefficient_names=names,multiplets=2,
        profiles=True,random_starts=5,extra_initial=np.array(list(selected['fit']['parameters'].values())))
    refined=prepare_matrices(varied,data,energy_order=8,pair_order=72,angle_order=72,
        source_breaking=False,control=True,multiplets=2,cutoff=2.5,a2=True)
    best=np.array(list(preferred['parameters'].values()))
    changes={channel:float(np.max(np.abs(predict_channel(best,refined,channel,True,5,2)
        -predict_channel(best,matrices,channel,True,5,2)))) for channel in CHANNELS}
    preferred['tensor_coefficient_products']=dict(
        gT_times_CT_dimensionless=.0233*preferred['parameters']['C_T_GeV_inverse'],
        gT_times_DT_GeV=.0233*preferred['parameters']['D_T'])
    return dict(bibtex_keys=['BaBar:2007ceh','Qin:2024ulb','Alda:2026nrd'],
        scope='Source-constrained inclusive BaBar VPPP fit with two excited vector multiplets, a2, and a discrete phi-prime shape scan',
        conventions='Canonical epsilon^0123=+1; all single-epsilon odd coefficients are minus source values',
        inputs=dict(constrained_parameters=parameters,F_GeV=F,operator_mass_GeV=MV,
            gT_GeV=.0233,a2_mass_GeV=1.3182,a2_width_GeV=.107,cutoff_GeV=2.5),
        assumptions=['Qin:2024ulb high-energy relations are fixed model assumptions, not likelihood data',
            'g1+2g2-g3=0; canonical g2=-Nc MV/(192sqrt2 pi² FV); FV GV=F²',
            'Canonical c/d constraints as in source_constrained_control of babar_kkpi_g_fit.json',
            'Tensor a2 terms in KS channel included; K2* tails omitted; CT units GeV^-1',
            'All first/second excited-vector beta residues are real and fitted separately in each charge channel',
            'Phi-prime mass/width scanned on an explicit bounded grid; no published same-data fit is used as an independent prior',
            'The grid minimum is a discrete model comparison, not a continuous pole confidence interval',
            'Inclusive Table2/3 only; separate 5.1%/6.5% common normalization nuisances; no Dalitz15%'],
        data={channel:{key:value.tolist() if isinstance(value,np.ndarray) else value for key,value in block.items()} for channel,block in data.items()},
        grid_fits=fits,selected_phi_prime_mass_GeV=selected['phi_prime_mass_GeV'],
        selected_phi_prime_width_GeV=selected['phi_prime_width_GeV'],
        preferred_grid_point_fit=preferred,quadrature_refinement_max_absolute_change_nb=changes,
        exact_null='delta g1=delta g3; imposed gA/g2 values are model assumptions')


def run_pole_grid():
    """Continue the inclusive fit with bounded rho'/omega' shape variations."""
    preliminary,data=load_inputs()
    parameters=control_parameters(preliminary)
    parameters.update(phi_prime_mass_GeV=1.640,phi_prime_width_GeV=.240)
    prior_path=ROOT/'data/vector_sector/babar_kkpi_g_fit_phi_scan.json'
    prior=json.loads(prior_path.read_text())
    continuation=np.array(list(prior['preferred_grid_point_fit']['parameters'].values()))
    names=('g4','2g4_plus_g5','c4','C_T_GeV_inverse','D_T')
    def fit_at_poles(varied,initial):
        matrices=prepare_matrices(varied,data,source_breaking=False,control=True,
            multiplets=2,cutoff=2.5,a2=True)
        return fit_coefficients(matrices,data,True,coefficient_names=names,
            multiplets=2,profiles=False,random_starts=5,extra_initial=initial)
    rho_grid=[]
    for mass in (1.465,1.510,1.550):
        for width in (.220,.300,.400):
            varied=dict(parameters,rho_prime_mass_GeV=mass,rho_prime_width_GeV=width)
            fit=fit_at_poles(varied,continuation)
            continuation=np.array(list(fit['parameters'].values()))
            rho_grid.append(dict(rho_prime_mass_GeV=mass,rho_prime_width_GeV=width,fit=fit))
            print('Rho-prime grid',mass,width,'chi2',fit['chi_squared'],flush=True)
    selected_rho=min(rho_grid,key=lambda item:item['fit']['chi_squared'])
    fixed_rho=dict(parameters,rho_prime_mass_GeV=selected_rho['rho_prime_mass_GeV'],
        rho_prime_width_GeV=selected_rho['rho_prime_width_GeV'])
    selected_parameters=dict(fixed_rho,omega_prime_mass_GeV=1.410,omega_prime_width_GeV=.290)
    selected_fit=selected_rho['fit']
    omega_grid=[]
    # This threshold only decides whether to execute the requested model
    # sensitivity. It is not a calibration acceptance theorem.
    if selected_fit['probability']<.01:
        continuation=np.array(list(selected_fit['parameters'].values()))
        for mass in (1.410,1.4223,1.430):
            for width in (.180,.240,.290):
                varied=dict(fixed_rho,omega_prime_mass_GeV=mass,omega_prime_width_GeV=width)
                fit=fit_at_poles(varied,continuation)
                continuation=np.array(list(fit['parameters'].values()))
                omega_grid.append(dict(omega_prime_mass_GeV=mass,omega_prime_width_GeV=width,fit=fit))
                print('Omega-prime grid',mass,width,'chi2',fit['chi_squared'],flush=True)
        selected_omega=min(omega_grid,key=lambda item:item['fit']['chi_squared'])
        if selected_omega['fit']['chi_squared']<selected_fit['chi_squared']:
            selected_parameters=dict(fixed_rho,omega_prime_mass_GeV=selected_omega['omega_prime_mass_GeV'],
                omega_prime_width_GeV=selected_omega['omega_prime_width_GeV'])
            selected_fit=selected_omega['fit']
    matrices=prepare_matrices(selected_parameters,data,source_breaking=False,control=True,
        multiplets=2,cutoff=2.5,a2=True)
    preferred=fit_coefficients(matrices,data,True,coefficient_names=names,multiplets=2,
        profiles=True,random_starts=5,extra_initial=np.array(list(selected_fit['parameters'].values())))
    refined=prepare_matrices(selected_parameters,data,energy_order=8,pair_order=72,angle_order=72,
        source_breaking=False,control=True,multiplets=2,cutoff=2.5,a2=True)
    best=np.array(list(preferred['parameters'].values()))
    changes={channel:float(np.max(np.abs(predict_channel(best,refined,channel,True,5,2)
        -predict_channel(best,matrices,channel,True,5,2)))) for channel in CHANNELS}
    preferred['tensor_coefficient_products']=dict(
        gT_times_CT_dimensionless=.0233*preferred['parameters']['C_T_GeV_inverse'],
        gT_times_DT_GeV=.0233*preferred['parameters']['D_T'])
    return dict(bibtex_keys=['BaBar:2007ceh','Qin:2024ulb','Alda:2026nrd'],
        scope='Inclusive VPPP fit with source constraints, a2, and bounded rho-prime/omega-prime pole grids',
        conventions='Canonical epsilon^0123=+1; odd coefficients are minus source convention',
        inputs=dict(constrained_parameters=parameters,F_GeV=F,operator_mass_GeV=MV,
            gT_GeV=.0233,a2_mass_GeV=1.3182,a2_width_GeV=.107,cutoff_GeV=2.5,
            initialization_artifact=str(prior_path.relative_to(ROOT)),
            initialization_sha256=hashlib.sha256(prior_path.read_bytes()).hexdigest()),
        model_assumptions=prior['assumptions']+[
            'rho-prime mass grid1.465,1.510,1.550GeV and width grid.220,.300,.400GeV; omega-prime grid1.410,1.4223,1.430GeV and .180,.240,.290GeV',
            'Phi-prime fixed1.640,.240GeV at preceding grid minimum',
            'Pole grids compare bounded model shapes; no published same-data fitted pole is used as a likelihood prior',
            'Omega-prime grid executed if best rho-grid nominal chi-square probability is below.01'],
        data={channel:{key:value.tolist() if isinstance(value,np.ndarray) else value for key,value in block.items()} for channel,block in data.items()},
        rho_prime_grid=rho_grid,omega_prime_grid=omega_grid,
        selected_pole_parameters={key:value for key,value in selected_parameters.items() if 'prime_' in key},
        preferred_grid_point_fit=preferred,quadrature_refinement_max_absolute_change_nb=changes,
        exact_null='delta g1=delta g3; gA=0 and g2 are imposed source constraints, not measurements')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--universal-source',action='store_true')
    parser.add_argument('--two-multiplets',action='store_true')
    parser.add_argument('--a2',action='store_true')
    parser.add_argument('--phi-prime-scan',action='store_true')
    parser.add_argument('--pole-grid',action='store_true')
    arguments=parser.parse_args()
    content=(run_pole_grid() if arguments.pole_grid else run_phi_prime_scan() if arguments.phi_prime_scan
        else run_two_multiplet_analysis(arguments.a2) if arguments.two_multiplets else run_analysis(not arguments.universal_source))
    arguments.output.write_text(json.dumps(content,indent=2,allow_nan=False)+'\n')
    print(arguments.output)
