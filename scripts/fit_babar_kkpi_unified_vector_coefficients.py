"""One inclusive BaBar likelihood for all ten visible odd-vector coordinates.

Run: .venv/bin/python scripts/fit_babar_kkpi_unified_vector_coefficients.py --output data/vector_sector/babar_kkpi_unified_vector_coefficients_fit.json
GeV, canonical epsilon^0123=+1. No OPE relation or odd-coefficient prior is
used. Source equations are Qin:2024ulb's full vector current, with the
preliminary light/strange photon-source residues inserted consistently.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares
from scipy.stats import chi2

try:
    from fit_babar_kkpi_g_coefficients import (
        load_inputs, form_factor, propagator, F, MV, MK, MPI, ALPHA,
        CHANNELS, triangle)
except ModuleNotFoundError:
    from scripts.fit_babar_kkpi_g_coefficients import (
        load_inputs, form_factor, propagator, F, MV, MK, MPI, ALPHA,
        CHANNELS, triangle)

ROOT=Path(__file__).resolve().parents[1]
CONVERSION=389379.3656
COORDINATE_NAMES=(
    'c1_minus_c2_plus_c5','c1_minus_c2_minus_c5_plus_2c6',
    'c1_plus_c2_plus_8c3_minus_c5','c4','d1_plus_8d2_minus_d3','d3',
    'g1_plus_2g2_minus_g3','g2','g4','2g4_plus_g5')
NUISANCE_NAMES=('gT_times_CT_dimensionless','gT_times_DT_GeV',
    'beta_first_kkpi0','beta_first_kskpi','beta_second_kkpi0','beta_second_kskpi',
    'normalization_pull_kkpi0','normalization_pull_kskpi')
PARAMETER_NAMES=COORDINATE_NAMES+NUISANCE_NAMES
RELATIVE_SYSTEMATICS={'kkpi0':.065,'kskpi':.051}
POWERS_LINEAR=((0,0),(1,0),(0,1))
POWERS_QUADRATIC=POWERS_LINEAR+((2,0),(1,1),(0,2))
BASIS_TERMS=[(-1,0,0)]
for coordinate in range(12):
    powers=(POWERS_QUADRATIC if coordinate in (4,5)
        else ((0,0),) if coordinate==10 else POWERS_LINEAR)
    BASIS_TERMS.extend((coordinate,first,second) for first,second in powers)


def original_coefficients(coordinates):
    """A representative of the six-dimensional exact coordinate-map null."""
    c_a,c_b,c_mass,c4,d_mass,d3,g_a,g2,g4,g_mass=coordinates[:10]
    return dict(c1=c_a,c2=0.,c3=(c_mass-c_a)/8,c4=c4,c5=0.,
        c6=(c_b-c_a)/2,c7=0.,d1=d_mass+d3,d2=0.,d3=d3,d4=0.,
        g1=g_a-2*g2,g2=g2,g3=0.,g4=g4,g5=g_mass-2*g4)


def preliminary_coordinates(parameters):
    return np.array([
        parameters['c1']-parameters['c2']+parameters['c5'],
        parameters['c1']-parameters['c2']-parameters['c5']+2*parameters['c6'],
        parameters['c1']+parameters['c2']+8*parameters['c3']-parameters['c5'],
        parameters['c4'],parameters['d1']+8*parameters['d2']-parameters['d3'],
        parameters['d3'],0.,0.,0.,0.])


def coherent_current(q_squared,pair_s,pair_t,even_parameters,channel,coordinates,beta):
    """All c,d,g coordinates enter this same electromagnetic scalar."""
    coefficients=original_coefficients(coordinates)
    varied=dict(even_parameters,**{key:value for key,value in coefficients.items() if not key.startswith('g')})
    result=form_factor(q_squared,pair_s,pair_t,varied,channel,
        [coordinates[6],coordinates[7],coordinates[8],coordinates[9]],beta,True)
    if channel=='kskpi':
        pair_u=q_squared+2*MK**2+MPI**2-pair_s-pair_t
        tensor_response=1/(1.3182**2-pair_s-1j*1.3182*.107)
        light_source=even_parameters['F_V_GeV']+2*even_parameters['z_V_GeV_inverse']*.036616
        result += ((-2/F**3)*coordinates[10]*(pair_t-pair_u)*tensor_response
            +(np.sqrt(2)*light_source/F**3)*coordinates[11]*(pair_t-pair_u)
            *tensor_response*propagator('rho',q_squared,even_parameters,beta))
    return result


def amplitude_basis(q_squared,pair_s,pair_t,parameters,channel):
    """Exact polynomial in coefficients and the two excited-multiplet residues."""
    zero=np.zeros(12)
    offset=coherent_current(q_squared,pair_s,pair_t,parameters,channel,zero,(0.,0.))
    columns=[offset]
    for coordinate in range(12):
        unit=np.eye(12)[coordinate]
        fixed=coherent_current(q_squared,pair_s,pair_t,parameters,channel,unit,(0.,0.))-offset
        if coordinate==10:
            columns.append(fixed)
            continue
        first_plus=coherent_current(q_squared,pair_s,pair_t,parameters,channel,unit,(1.,0.))-offset
        second_plus=coherent_current(q_squared,pair_s,pair_t,parameters,channel,unit,(0.,1.))-offset
        if coordinate in (4,5):
            first_minus=coherent_current(q_squared,pair_s,pair_t,parameters,channel,unit,(-1.,0.))-offset
            second_minus=coherent_current(q_squared,pair_s,pair_t,parameters,channel,unit,(0.,-1.))-offset
            first_linear=(first_plus-first_minus)/2
            second_linear=(second_plus-second_minus)/2
            first_square=(first_plus+first_minus)/2-fixed
            second_square=(second_plus+second_minus)/2-fixed
            mixed=(coherent_current(q_squared,pair_s,pair_t,parameters,channel,unit,(1.,1.))
                -offset-fixed-first_linear-second_linear-first_square-second_square)
            columns.extend((fixed,first_linear,second_linear,first_square,mixed,second_square))
        else:
            columns.extend((fixed,first_plus-fixed,second_plus-fixed))
    return np.stack(columns,axis=-1)


def coefficients_and_derivatives(values,channel):
    """Amplitude coefficient vector and its exact derivative in fitted coordinates."""
    channel_index=CHANNELS.index(channel)
    first_index,second_index=12+channel_index,14+channel_index
    beta_first,beta_second=values[first_index],values[second_index]
    vector=np.zeros(len(BASIS_TERMS))
    derivative=np.zeros((len(BASIS_TERMS),len(PARAMETER_NAMES)))
    for row,(coordinate,first_power,second_power) in enumerate(BASIS_TERMS):
        if coordinate<0:
            vector[row]=1.
            continue
        monomial=beta_first**first_power*beta_second**second_power
        vector[row]=values[coordinate]*monomial
        derivative[row,coordinate]=monomial
        if first_power:
            derivative[row,first_index]=values[coordinate]*first_power*beta_first**(first_power-1)*beta_second**second_power
        if second_power:
            derivative[row,second_index]=values[coordinate]*second_power*beta_first**first_power*beta_second**(second_power-1)
    return vector,derivative


def prepare_matrices(parameters,data,energy_order=6,pair_order=48,angle_order=48,cutoff=2.5):
    """Integrate all coherent products with the existing absolute normalization."""
    energy_nodes,energy_weights=leggauss(energy_order)
    pair_nodes,pair_weights=leggauss(pair_order)
    angle_nodes,angle_weights=leggauss(angle_order)
    matrices={}
    for channel in CHANNELS:
        masses=(.493677,.493677,.1349768) if channel=='kkpi0' else (.497611,.493677,.13957039)
        first_mass,second_mass,third_mass=masses
        output=[]
        for lower,upper in data[channel]['edges']:
            matrix=np.zeros((len(BASIS_TERMS),len(BASIS_TERMS)))
            for energy_node,energy_weight in zip(energy_nodes,energy_weights):
                energy=(lower+upper)/2+(upper-lower)*energy_node/2
                q_squared=energy**2
                pair_lower=(first_mass+second_mass)**2
                pair_upper=(energy-third_mass)**2
                pair_s=(pair_lower+pair_upper)/2+(pair_upper-pair_lower)*pair_nodes/2
                outer=np.maximum(triangle(q_squared,pair_s,third_mass**2),0)
                inner=np.maximum(triangle(pair_s,first_mass**2,second_mass**2),0)
                span=np.sqrt(outer*inner)/pair_s
                pair_s13=(first_mass**2+third_mass**2+(pair_s+first_mass**2-second_mass**2)
                    *(q_squared-pair_s-third_mass**2)/(2*pair_s))[:,None]+span[:,None]*angle_nodes/2
                pair_s23=q_squared+sum(mass**2 for mass in masses)-pair_s[:,None]-pair_s13
                basis=amplitude_basis(q_squared,pair_s[:,None],pair_s23,parameters,channel)
                basis *= np.exp(-q_squared**3/cutoff**6)
                current_squared=(outer*inner/(16*pair_s))[:,None]*(1-angle_nodes**2)
                measure=((pair_upper-pair_lower)/2*pair_weights[:,None]
                    *angle_weights[None,:]*span[:,None]/2)
                weight=CONVERSION*ALPHA**2/(48*np.pi*q_squared**3)*measure*current_squared*energy_weight/2
                flat=basis.reshape(-1,len(BASIS_TERMS));flat_weight=weight.reshape(-1)
                matrix += (flat.conj().T@(flat_weight[:,None]*flat)).real
            output.append(matrix*(2 if channel=='kskpi' else 1))
        matrices[channel]=np.asarray(output)
    return matrices


def predictions_and_derivatives(values,matrices,channel):
    vector,derivative=coefficients_and_derivatives(values,channel)
    contracted=np.einsum('bij,j->bi',matrices[channel],vector)
    return contracted@vector,2*contracted@derivative


def likelihood(values,matrices,data):
    residuals=[];jacobians=[]
    for channel_index,channel in enumerate(CHANNELS):
        prediction,derivative=predictions_and_derivatives(values,matrices,channel)
        normalization_index=16+channel_index
        scale=1+RELATIVE_SYSTEMATICS[channel]*values[normalization_index]
        residuals.extend((scale*prediction-data[channel]['observed'])/data[channel]['errors'])
        jacobian=scale*derivative/data[channel]['errors'][:,None]
        jacobian[:,normalization_index] += RELATIVE_SYSTEMATICS[channel]*prediction/data[channel]['errors']
        jacobians.extend(jacobian)
    for index in (16,17):
        residuals.append(values[index])
        row=np.zeros(len(values));row[index]=1.;jacobians.append(row)
    return np.asarray(residuals),np.asarray(jacobians)


def coordinate_nulls():
    c=sp.symbols('c1:8');d=sp.symbols('d1:5');g=sp.symbols('g1:6')
    coordinates=sp.Matrix([c[0]-c[1]+c[4],c[0]-c[1]-c[4]+2*c[5],
        c[0]+c[1]+8*c[2]-c[4],c[3],d[0]+8*d[1]-d[2],d[2],
        g[0]+2*g[1]-g[2],g[1],g[3],2*g[3]+g[4]])
    matrix=coordinates.jacobian(c+d+g)
    original_order=list(map(str,c+d+g))
    named=({'c7':1},{'d4':1},{'g1':1,'g3':1},{'d1':-8,'d2':1},
        {'c2':1,'c5':1,'c6':1},{'c1':1,'c2':1,'c3':sp.Rational(-1,4)})
    named_records=[]
    for changes in named:
        vector=sp.Matrix([changes.get(name,0) for name in original_order])
        named_records.append(dict(changes={key:str(value) for key,value in changes.items()},
            exact_coordinate_residual_zero=matrix*vector==sp.zeros(10,1)))
    return dict(original_order=original_order,coordinate_map=[[str(x) for x in row] for row in matrix.tolist()],
        exact_rank=int(matrix.rank()),named_exact_nulls=named_records,
        exact_null_basis=[[str(value) for value in vector] for vector in matrix.nullspace()])


def zero_width_amplitude_certificate(parameters):
    """Exact rational pole map before using any finite-width rate information.

    Columns c/g are divided by sqrt(2) to work entirely over QQ; the
    column rescaling is inverted in the returned physical null vectors.
    The excited-pole comparison isolates excitation-model lifting from
    width lifting. A zero-width three-body rate integral is not defined
    across an on-shell K* pole, so this certificate is at amplitude level.
    """
    rational=lambda value:sp.Rational(str(value))
    f,mv,mk2,mpi2=map(rational,(F,MV,MK**2,MPI**2))
    fv,gv,z=map(rational,(parameters['F_V_GeV'],parameters['G_V_GeV'],parameters['z_V_GeV_inverse']))
    light=fv+2*z*rational(.036616);strange=fv+2*z*rational(.943684)
    masses={'rho':rational(parameters['m_rho_GeV']),'omega':rational(.78266),
        'phi':rational(1.019461),'kstar':rational(parameters['m_Kstar_GeV'])}
    first={'rho':rational(1.465),'omega':rational(1.410),'phi':rational(1.680),'kstar':rational(1.414)}
    second={'rho':rational(1.720),'omega':rational(1.670),'phi':rational(2.162),'kstar':rational(1.718)}
    def row(q,s,t,channel,excited):
        u=q+2*mk2+mpi2-s-t
        def pole(name,x):
            value=1/(masses[name]**2-x)
            if excited:
                value += -sp.Rational(1,3)/(first[name]**2-x)+sp.Rational(1,20)/(second[name]**2-x)
            return value
        def cr1(pair,mass):
            result=sp.zeros(1,10);result[0]=q;result[1]=-pair;result[2]=mass
            return result
        def cr2(pair):
            result=cr1(pair,mk2);result[3]=24*(mk2-mpi2)
            return result
        def dr(pair,mass):
            result=sp.zeros(1,10);result[4]=mass;result[5]=q+pair
            return result
        gr1=sp.zeros(1,10);gr2=sp.zeros(1,10)
        gr1[6]=q+s-4*mk2-mpi2;gr1[7]=-4*(q-2*mk2-mpi2)
        gr1[8]=-4*(mk2-mpi2);gr1[9]=4*mk2
        gr2[6]=q-s-mpi2;gr2[7]=-2*(q-2*mk2-mpi2)
        gr2[8]=4*(mk2-mpi2);gr2[9]=2*mpi2
        rq,wq,pq=light*pole('rho',q),light*pole('omega',q),strange*pole('phi',q)
        rs,ws,kt,ku=pole('rho',s),pole('omega',s),pole('kstar',t),pole('kstar',u)
        if channel=='kkpi0':
            local=2/(f**3*mv)*((rq+wq/3)*gr1-2*pq*gr2/3)
            direct=-2*gv/(3*f**3*mv)*((rs+3*ws)*cr1(s,mpi2)+kt*cr2(t)+ku*cr2(u))
            exchange=2*gv/f**3*(2*(rq*ws+wq*rs/3)*dr(s,mpi2)
                +(rq+wq/3-2*pq/3)*(kt*dr(t,mk2)+ku*dr(u,mk2)))
        else:
            antisymmetric=sp.zeros(1,10);antisymmetric[6]=(t-u)*rq
            local=2/(f**3*mv)*(antisymmetric-wq*gr1/3+2*pq*gr2/3)
            direct=2*gv/(3*f**3*mv)*(rs*cr1(s,mpi2)-2*kt*cr1(t,mk2)+ku*cr2(u))
            exchange=-2*gv/f**3*(2*wq*rs*dr(s,mpi2)/3
                -(rq-wq/3+2*pq/3)*kt*dr(t,mk2)
                +(rq+wq/3-2*pq/3)*ku*dr(u,mk2))
        # direct/local already have sqrt(2) removed; exchange did not.
        return local+direct+exchange
    points=[(q,s,t,channel) for q in (sp.Rational(49,25),sp.Rational(64,25),sp.Rational(81,25))
        for s in (sp.Rational(21,20),sp.Rational(6,5))
        for t in (sp.Rational(11,20),sp.Rational(7,10)) for channel in CHANNELS]
    output={}
    for excited in (False,True):
        matrix=sp.Matrix.vstack(*(row(*point,excited) for point in points))
        domain=sp.polys.matrices.DomainMatrix.from_Matrix(matrix)
        rank=int(domain.rank())
        _,pivot_rows=domain.transpose().rref()
        minor_nonzero=None
        if rank==10:
            minor=domain.extract(list(pivot_rows[:10]),list(range(10)))
            minor_nonzero=minor.det()!=domain.domain.zero
        null_rows=domain.nullspace().to_Matrix()
        physical_nulls=[]
        for vector in null_rows.tolist():
            # QQ columns correspond to A/sqrt2 for c,g and A for d.
            physical=[value/(sp.sqrt(2) if index not in (4,5) else 1)
                for index,value in enumerate(vector)]
            scale=next(value for value in physical if value!=0)
            physical_nulls.append([sp.simplify(value/scale) for value in physical])
        beta=(-1/3,1/20) if excited else (0.,0.)
        finite_rows=[]
        for q,s,t,channel in points:
            q,s,t=map(float,(q,s,t))
            zero=np.zeros(12)
            offset=coherent_current(q,np.array([[s]]),np.array([[t]]),parameters,channel,zero,beta)[0,0]
            gradient=[]
            for index in range(10):
                unit=np.zeros(12);unit[index]=1
                gradient.append(coherent_current(q,np.array([[s]]),np.array([[t]]),parameters,channel,unit,beta)[0,0]-offset)
            finite_rows.extend((np.real(gradient),np.imag(gradient)))
        finite_matrix=np.array(finite_rows)
        finite_singular=np.linalg.svd(finite_matrix,compute_uv=False)
        width_lifts=[]
        for null in physical_nulls:
            vector=np.array([float(sp.N(value)) for value in null])
            width_lifts.append(float(np.linalg.norm(finite_matrix@vector)
                /(np.linalg.norm(finite_matrix)*np.linalg.norm(vector))))
        output['two_excited_multiplets' if excited else 'ground_poles']=dict(
            exact_rational_amplitude_rank=rank,number_of_rows=matrix.rows,
            exact_pivot_row_indices=list(map(int,pivot_rows)),
            exact_full_rank_minor_nonzero=bool(minor_nonzero) if minor_nonzero is not None else None,
            evaluation_points=[[str(value) for value in point] for point in points],
            rational_input_policy='Every declared decimal mass/source input is interpreted as its exact decimal rational; c/g columns divided by sqrt2',
            exact_physical_coefficient_nulls=[[str(value) for value in vector] for vector in physical_nulls],
            finite_width_amplitude_singular_values=finite_singular.tolist(),
            finite_width_amplitude_rank=int(sum(finite_singular>finite_singular[0]*1e-8)),
            relative_finite_width_response_along_zero_width_nulls=width_lifts)
    return output


def fit_all_coefficients(parameters,data,matrices,starting_points=None):
    initial=np.r_[preliminary_coordinates(parameters),0.,0.,-.3,-.4,0.,.04,0.,0.]
    lower=np.r_[np.full(10,-3.),-1.,-1.,np.full(4,-2.),-8.,-8.]
    upper=-lower
    scales=np.r_[np.full(4,.02),np.full(6,.1),.01,.02,np.full(4,.3),1.,1.]
    random=np.random.default_rng(7511)
    starts=[initial] if starting_points is None else list(starting_points)
    if starting_points is None:
        for index in range(23):
            candidate=initial.copy()
            candidate[:10] += random.normal(0,scales[:10])
            candidate[10:12]=random.normal(0,[.01,.03])
            candidate[12:16]=random.uniform(-.7,.2,4)
            starts.append(candidate)
    fits=[]
    for initial_point in starts:
        result=least_squares(lambda values:likelihood(values,matrices,data)[0],initial_point,
            jac=lambda values:likelihood(values,matrices,data)[1],bounds=(lower,upper),
            x_scale=scales,max_nfev=4000,ftol=1e-10,xtol=1e-10,gtol=1e-9)
        fits.append(result)
    result=min(fits,key=lambda value:value.fun@value.fun)
    pre_polish_statistic=float(result.fun@result.fun)
    polishing=[]
    for iteration in range(2):
        improved=least_squares(lambda values:likelihood(values,matrices,data)[0],result.x,
            jac=lambda values:likelihood(values,matrices,data)[1],bounds=(lower,upper),
            x_scale='jac',max_nfev=20000,ftol=1e-9,xtol=1e-9,gtol=1e-8)
        polishing.append(dict(iteration=iteration+1,chi_squared=float(improved.fun@improved.fun),
            nfev=improved.nfev,success=bool(improved.success),message=improved.message,
            optimality=float(improved.optimality)))
        fits.append(improved)
        if improved.fun@improved.fun<=result.fun@result.fun:
            result=improved
        if improved.success:
            break
    statistic=float(result.fun@result.fun)
    jacobian=result.jac
    nuisance=jacobian[:,10:]
    profiled=jacobian[:,:10]-nuisance@np.linalg.pinv(nuisance,rcond=1e-12)@jacobian[:,:10]
    singular_values=np.linalg.svd(profiled,compute_uv=False)
    rank_by_tolerance={str(tolerance):int(sum(singular_values>singular_values[0]*tolerance))
        for tolerance in (1e-6,1e-7,1e-8,1e-9,1e-10)}
    covariance=np.linalg.pinv(jacobian.T@jacobian,rcond=1e-12)
    predictions={};blocks={}
    for channel_index,channel in enumerate(CHANNELS):
        raw,_=predictions_and_derivatives(result.x,matrices,channel)
        scale=1+RELATIVE_SYSTEMATICS[channel]*result.x[16+channel_index]
        predictions[channel]=(raw*scale).tolist()
        blocks[channel]=float(np.sum(((raw*scale-data[channel]['observed'])/data[channel]['errors'])**2))
    minimum_branches=[]
    for item in sorted(fits,key=lambda value:value.fun@value.fun):
        value=float(item.fun@item.fun)
        if not minimum_branches or all(abs(value-entry['chi_squared'])>1e-3 for entry in minimum_branches):
            minimum_branches.append(dict(chi_squared=value,parameters=dict(zip(PARAMETER_NAMES,map(float,item.x))),
                optimizer_success=bool(item.success)))
    return dict(parameters=dict(zip(PARAMETER_NAMES,map(float,result.x))),
        fitted_original_representative=original_coefficients(result.x),chi_squared=statistic,
        nominal_degrees_of_freedom=45-16,nominal_probability=float(chi2.sf(statistic,45-16)),
        integrated_rate_profiled_singular_values=singular_values.tolist(),
        integrated_rate_rank_by_relative_tolerance=rank_by_tolerance,
        local_covariance=covariance.tolist(),local_errors=dict(zip(PARAMETER_NAMES,map(float,np.sqrt(np.diag(covariance))))),
        active_numerical_bounds=[PARAMETER_NAMES[i] for i in range(len(result.x))
            if min(result.x[i]-lower[i],upper[i]-result.x[i])<1e-5*(upper[i]-lower[i])],
        optimizer_success=bool(result.success),optimizer_message=result.message,
        optimizer_nfev=int(result.nfev),optimizer_optimality=float(result.optimality),
        unscaled_gradient_infinity_norm=float(np.max(np.abs(result.jac.T@result.fun))),
        polishing=polishing,
        multistart_minima=minimum_branches,number_of_starts=len(starts),global_minimum_proven=False,
        executed_initial_multistart_search=dict(number_of_starts=len(starts),
            best_pre_polish_chi_squared=pre_polish_statistic,
            qualification='Initial deterministic starts followed by analytic-Jacobian polishing; global minimality is not proven'),
        predictions_nb=predictions,channel_chi_squared=blocks,
        coefficient_error_qualification='Local linear diagnostics; exact and width-scheme nulls must be assessed before interpreting individual errors',
        profile_likelihoods_completed=False),result.x


def profile_g_coordinates(values,matrices,data,minimum,coordinates=(8,9,6,7)):
    """Profile every other fitted coefficient and nuisance for four VPPP coordinates."""
    lower=np.r_[np.full(10,-3.),-1.,-1.,np.full(4,-2.),-8.,-8.]
    upper=-lower
    profiles={}
    for coordinate in coordinates:
        active=[index for index in range(len(values)) if index!=coordinate]
        grid=np.unique(np.r_[[-3.,-1.5,0.,1.5,3.],values[coordinate]])
        order=sorted(range(len(grid)),key=lambda index:abs(grid[index]-values[coordinate]))
        solved={}
        for point_index in order:
            fixed=grid[point_index]
            def expanded(free):
                full=values.copy();full[active]=free;full[coordinate]=fixed
                return full
            def residual(free):
                return likelihood(expanded(free),matrices,data)[0]
            def jacobian(free):
                return likelihood(expanded(free),matrices,data)[1][:,active]
            initial_points=[values[active]]
            if solved:
                nearest=min(solved,key=lambda index:abs(grid[index]-fixed))
                initial_points.append(solved[nearest][1].x)
            fits=[least_squares(residual,initial,jac=jacobian,bounds=(lower[active],upper[active]),
                x_scale='jac',max_nfev=1500,ftol=1e-6,xtol=1e-7,gtol=1e-6)
                for initial in initial_points]
            result=min(fits,key=lambda item:item.fun@item.fun)
            solved[point_index]=(expanded(result.x),result)
        records=[]
        for point_index,fixed in enumerate(grid):
            full,result=solved[point_index]
            records.append(dict(fixed_value=float(fixed),chi_squared=float(result.fun@result.fun),
                optimizer_success=bool(result.success),optimizer_message=result.message,nfev=int(result.nfev),
                optimality=float(result.optimality),
                active_other_parameter_bounds=[PARAMETER_NAMES[index] for index in active
                    if min(full[index]-lower[index],upper[index]-full[index])<1e-4*(upper[index]-lower[index])],
                profiled_parameters=dict(zip(PARAMETER_NAMES,map(float,full)))))
        statistics=np.array([item['chi_squared'] for item in records])
        crossings={}
        for level in (1.,3.841458820694124):
            difference=statistics-minimum-level
            crossings[str(level)]=[float(grid[index]+(grid[index+1]-grid[index])
                *(-difference[index])/(difference[index+1]-difference[index]))
                for index in range(len(grid)-1) if difference[index]*difference[index+1]<0]
        profiles[COORDINATE_NAMES[coordinate]]=dict(points=records,
            delta_chi_squared_crossings_linear_grid=crossings,
            all_points_converged=all(item['optimizer_success'] for item in records),
            scan_lower_endpoint_delta_chi_squared=float(statistics[0]-minimum),
            scan_upper_endpoint_delta_chi_squared=float(statistics[-1]-minimum),
            minimum_profile_chi_squared_minus_reported_minimum=float(min(statistics)-minimum),
            numeric_endpoints_are_confidence_limits=False,
            confidence_intervals_certified=False,
            scan_qualification='Bounded coarse profile; per-point convergence and other active bounds retained; interpolated crossings are diagnostics',
            numerical_coefficient_domain=[float(lower[coordinate]),float(upper[coordinate])])
        print('Unified profile',COORDINATE_NAMES[coordinate],min(statistics),
            'endpoint delta',statistics[0]-minimum,statistics[-1]-minimum,flush=True)
    return profiles


def run_analysis():
    parameters,data=load_inputs()
    matrices=prepare_matrices(parameters,data)
    fit,values=fit_all_coefficients(parameters,data,matrices)
    print('Unified all-ten-coordinate fit',fit['chi_squared'],fit['parameters'],flush=True)
    # A direct finite-difference comparison checks the analytic likelihood
    # derivative used in the numerical optimization and rank calculation.
    residual,jacobian=likelihood(values,matrices,data)
    differences=[]
    for index in range(len(values)):
        step=1e-6*max(abs(values[index]),.01)
        plus=values.copy();minus=values.copy();plus[index]+=step;minus[index]-=step
        numerical=(likelihood(plus,matrices,data)[0]-likelihood(minus,matrices,data)[0])/(2*step)
        differences.append(float(np.linalg.norm(numerical-jacobian[:,index])/max(1.,np.linalg.norm(jacobian[:,index]))))
    pole_certificate=zero_width_amplitude_certificate(parameters)
    fit['conservative_rank_ceiling_from_ground_zero_width']=pole_certificate['ground_poles']['exact_rational_amplitude_rank']
    fit['width_or_excited_model_lifting_requires_separate_interpretation']=(
        max(fit['integrated_rate_rank_by_relative_tolerance'].values())
        >pole_certificate['ground_poles']['exact_rational_amplitude_rank'])
    profiles=profile_g_coordinates(values,matrices,data,fit['chi_squared'])
    fit['profile_likelihoods_completed']=all(item['all_points_converged'] for item in profiles.values())
    refined=prepare_matrices(parameters,data,energy_order=8,pair_order=72,angle_order=72)
    refinement={channel:float(np.max(np.abs(predictions_and_derivatives(values,refined,channel)[0]
        -predictions_and_derivatives(values,matrices,channel)[0]))) for channel in CHANNELS}
    return dict(bibtex_keys=['BaBar:2007ceh','Qin:2024ulb','Alda:2026nrd','Ruiz-Femenia:2003jdx'],
        scope='One inclusive likelihood with all ten c,d,g coefficient coordinates varied together',
        source_equations=['Qin:2024ulb eq:fr','Qin:2024ulb eq:kpkmpi0','Qin:2024ulb eq:kskpi',
            'Qin:2024ulb Eq:CRDR','Qin:2024ulb Eq:app:Gamma','Qin:2024ulb eq:bw'],
        source_tex_sha256=hashlib.sha256((ROOT/'memory/sources/Qin:2024ulb/KKpi.tex').read_bytes()).hexdigest(),
        conventions='Canonical epsilon^0123=+1; source WZW sign reversed and linear odd gradients preserved',
        parameter_order=list(PARAMETER_NAMES),coordinate_map_and_original_nulls=coordinate_nulls(),
        inputs=dict(preliminary_parameters=parameters,F_GeV=F,operator_mass_GeV=MV,cutoff_GeV=2.5,
            fixed_even_source_inputs={key:value for key,value in parameters.items() if not key.startswith(('c','d'))}),
        assumptions=['No OPE relation or odd-coefficient prior in this unified baseline',
            'FV,GV,zV and other even inputs fixed at current preliminary point; all ten odd coordinates varied simultaneously',
            'Ideal neutral-vector mixing; rho/omega source poles use F_l and phi source poles use F_s',
            'Two excited vector multiplets with real beta residues profiled separately per channel; independent PDG pole shapes fixed',
            'Open a2 tensor terms profiled through products gT*CT and gT*DT; closed K2* tails omitted',
            'Inclusive Tables2/3 only, separate 6.5%/5.1% Gaussian normalization nuisances; no Table7 likelihood rows',
            'Uniform energy-bin averages, physical phase-space masses, averaged vertex masses'],
        data={channel:{key:value.tolist() if isinstance(value,np.ndarray) else value for key,value in block.items()} for channel,block in data.items()},
        fit=fit,zero_width_pole_amplitude_certificate=pole_certificate,
        g_profile_likelihoods=profiles,quadrature_refinement_max_absolute_change_nb=refinement,
        coefficient_confidence_intervals_certified=False,
        analytic_jacobian_max_relative_difference=max(differences))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    arguments=parser.parse_args()
    arguments.output.write_text(json.dumps(run_analysis(),indent=2,allow_nan=False)+'\n')
    print(arguments.output)
