"""Two-pole pion TFF constraints, with the canonical signed scalar.

Invocation: .venv/bin/python scripts/derive_improved_pion_tff.py
Units: GeV; f_pi=0.093 GeV in the Alda:2026nrd convention; chiral pion, tree meromorphic approximation, zero widths.
Inputs: memory/rcht/odd-parity/radiative-tff-amplitude.md and
memory/extracts/Knecht:2001qf/pion-tff-short-distance.md.
The f_pi=93 MeV convention is from
memory/extracts/Alda:2026nrd/neutral-meson-generators.md.
The leading-twist convolution is Gerardin:2019vio, Eqs. (DIS),
(TFF_asympt), in memory/sources/Gerardin:2019vio/main.tex.
Eta and eta-prime asymptotic coefficients: Bali:2021qem Eq. (5.80),
memory/extracts/Bali:2021qem/eta-eta-prime-tff-hard-inputs.md.
Sourced dispersive, asymptotic and effective-pole components follow
Hoferichter:2018kwz Eqs. (4.1), (4.9), (4.10), as recorded in
memory/extracts/Hoferichter:2018kwz/pion-tff-dispersive-asymptotic-representation.md.
That source states the OPE domain |omega|<1/2; its singly-virtual
|omega|=1 evaluation is a formal endpoint, not a consequence within
that OPE domain. The separate Brodsky-Lepage constraint is stated explicitly.
This script derives a symbolic source-based hard completion and an
illustrative interpolation. Unspecified low-energy spectral inputs prevent
its use as a physical finite-virtuality predictor.
"""
import sympy as sp


def show(label, expression):
    """Print a deterministic intermediate result using SymPy LaTeX."""
    print(label + ':\n$$\n' + sp.latex(expression) + '\n$$\n')


def derive():
    """Derive the meromorphic constraints and low-energy coefficient maps."""
    first_virtuality, second_virtuality = sp.symbols('s_1 s_2', real=True)
    hard_scale, hard_ratio = sp.symbols('Q_squared r', positive=True)
    decay_constant, vector_decay_constant, vector_mass = sp.symbols('f_pi F_V M_V', positive=True)
    first_mass, second_mass = sp.symbols('M_1 M_2', positive=True)
    colors = sp.Symbol('N_C', positive=True)
    single_coefficient = sp.Symbol('C_pi_single', real=True)
    h1, h2, h5, h7 = sp.symbols('h_1 h_2 h_5 h_7', real=True)
    numerator = (first_virtuality*second_virtuality*(first_virtuality+second_virtuality)
                 + h1*(first_virtuality+second_virtuality)**2
                 + h2*first_virtuality*second_virtuality
                 + h5*(first_virtuality+second_virtuality) + h7)
    denominator = sp.prod((mass**2-virtuality) for mass in (first_mass, second_mass)
                          for virtuality in (first_virtuality, second_virtuality))
    form_factor = decay_constant*numerator/(3*denominator)
    anomaly = -colors/(12*sp.pi**2*decay_constant)
    anomaly_solution = sp.solve(form_factor.subs({first_virtuality:0, second_virtuality:0})-anomaly, h7)[0]
    single_constant = sp.limit(form_factor.subs({first_virtuality:-hard_scale, second_virtuality:0}), hard_scale, sp.oo)
    h1_solution = sp.solve(single_constant, h1)[0]
    single_limit = sp.limit(hard_scale*form_factor.subs({first_virtuality:-hard_scale, second_virtuality:0, h1:h1_solution}),hard_scale,sp.oo)
    h5_solution = sp.solve(single_limit-single_coefficient,h5)[0]
    constraints = {h1:h1_solution,h5:h5_solution,h7:anomaly_solution}
    constrained = form_factor.subs(constraints)
    symmetric_limit = sp.limit(hard_scale*constrained.subs({first_virtuality:-hard_scale,second_virtuality:-hard_scale}),hard_scale,sp.oo)
    ray_limit = sp.limit(hard_scale*constrained.subs({first_virtuality:-hard_ratio*hard_scale,second_virtuality:-hard_scale}),hard_scale,sp.oo)
    slope = sp.factor(sp.diff(constrained,first_virtuality).subs({first_virtuality:0,second_virtuality:0}))
    mixed = sp.factor(sp.diff(constrained,first_virtuality,second_virtuality).subs({first_virtuality:0,second_virtuality:0}))
    c1,c2,c5,c6,d3 = sp.symbols('c_1 c_2 c_5 c_6 d_3',real=True)
    flavor_trace = sp.sqrt(2)/3
    original = -sp.sqrt(2)*colors*flavor_trace/(8*sp.pi**2*decay_constant)
    original += 4*vector_decay_constant*flavor_trace/(decay_constant*vector_mass)*(
        ((-c1+c2+c5-2*c6)*first_virtuality+(c1-c2+c5)*second_virtuality)/(vector_mass**2-first_virtuality)
        + ((-c1+c2+c5-2*c6)*second_virtuality+(c1-c2+c5)*first_virtuality)/(vector_mass**2-second_virtuality))
    original -= 4*sp.sqrt(2)*vector_decay_constant**2*flavor_trace*d3*(first_virtuality+second_virtuality)/(decay_constant*(vector_mass**2-first_virtuality)*(vector_mass**2-second_virtuality))
    original_slope = sp.factor(sp.diff(original,first_virtuality).subs({first_virtuality:0,second_virtuality:0}))
    original_mixed = sp.factor(sp.diff(original,first_virtuality,second_virtuality).subs({first_virtuality:0,second_virtuality:0}))
    slope_relation = sp.Eq(original_slope,slope)
    h2_solution = sp.solve(mixed-original_mixed,h2)[0]
    # Gerardin:2019vio uses the positive scalar; reverse its overall sign
    # to retain this repository's negative anomaly and common tensor phase.
    momentum_fraction = sp.Symbol('x', real=True)
    distribution_amplitude = 6*momentum_fraction*(1-momentum_fraction)
    asymmetry = (1-hard_ratio)/(1+hard_ratio)
    convolution = sp.Integral(distribution_amplitude/(1+asymmetry*(1-2*momentum_fraction)),
                              (momentum_fraction,0,1))
    # The equivalent denominator x*r+1-x keeps positivity explicit for r>0.
    leading_twist_integral = sp.integrate(distribution_amplitude/
        (momentum_fraction*hard_ratio+1-momentum_fraction),(momentum_fraction,0,1))
    leading_twist_ray = sp.factor(-2*decay_constant*leading_twist_integral/3)
    hard_ray_difference = sp.factor(ray_limit-leading_twist_ray)
    leading_twist_symmetric = sp.limit(leading_twist_ray,hard_ratio,1)
    leading_twist_small_ratio = sp.limit(leading_twist_ray,hard_ratio,0,dir='+')
    leading_twist_large_ratio = sp.limit(leading_twist_ray,hard_ratio,sp.oo)
    leading_twist_rescaled_large_ratio = sp.limit(hard_ratio*leading_twist_ray,hard_ratio,sp.oo)
    two_pole_small_ratio = sp.limit(ray_limit,hard_ratio,0,dir='+')
    two_pole_large_ratio = sp.limit(ray_limit,hard_ratio,sp.oo)
    leading_twist_h5 = sp.factor(h5_solution.subs(single_coefficient,leading_twist_small_ratio))
    normalized_hard_ray = sp.factor(leading_twist_ray/leading_twist_small_ratio)
    eta_coefficient, eta_prime_coefficient = sp.symbols('C_eta_single C_etap_single',real=True)
    octet_decay_constant, singlet_decay_constant = sp.symbols('F_P_8 F_P_0_infinity_Nf4',real=True)
    bali_coefficient = -2*(octet_decay_constant+2*sp.sqrt(2)*singlet_decay_constant)/sp.sqrt(3)
    channel_hard_rays = sp.Tuple(single_coefficient*normalized_hard_ray,
        eta_coefficient*normalized_hard_ray, eta_prime_coefficient*normalized_hard_ray)
    eta_hard_ray = channel_hard_rays[1].subs(eta_coefficient,-sp.Rational(1605,10000))
    eta_prime_hard_ray = channel_hard_rays[2].subs(eta_prime_coefficient,-sp.Rational(2305,10000))
    eta_symmetric_coefficient = sp.limit(eta_hard_ray,hard_ratio,1)
    eta_prime_symmetric_coefficient = sp.limit(eta_prime_hard_ray,hard_ratio,1)
    eta_symmetric_uncertainty = sp.Rational(100,10000)*sp.limit(normalized_hard_ray,hard_ratio,1)
    eta_prime_symmetric_uncertainty = sp.Rational(101,10000)*sp.limit(normalized_hard_ray,hard_ratio,1)
    # Illustrative finite-virtuality interpolation, not a sourced hadronic
    # spectral density. M1 is an unfixed smoothing scale in this expression.
    shifted_convolution = -2*decay_constant*sp.Integral(distribution_amplitude/
        (first_mass**2-momentum_fraction*first_virtuality-(1-momentum_fraction)*second_virtuality),
        (momentum_fraction,0,1))/3
    shifted_normalization = -2*decay_constant*sp.integrate(distribution_amplitude/first_mass**2,(momentum_fraction,0,1))/3
    shifted_slope = -2*decay_constant*sp.integrate(distribution_amplitude*momentum_fraction/first_mass**4,(momentum_fraction,0,1))/3
    shifted_mixed = -2*decay_constant*sp.integrate(distribution_amplitude*2*momentum_fraction*(1-momentum_fraction)/first_mass**6,(momentum_fraction,0,1))/3
    correction_response = sp.prod((first_mass**2/(first_mass**2-virtuality))**3
        for virtuality in (first_virtuality,second_virtuality))
    response_slope = sp.diff(correction_response,first_virtuality).subs({first_virtuality:0,second_virtuality:0})
    response_mixed = sp.diff(correction_response,first_virtuality,second_virtuality).subs({first_virtuality:0,second_virtuality:0})
    normalization_difference = anomaly-shifted_normalization
    slope_difference = original_slope-shifted_slope
    mixed_difference = original_mixed-shifted_mixed
    constant_correction = normalization_difference
    linear_correction = sp.factor(slope_difference-response_slope*constant_correction)
    mixed_correction = sp.factor(mixed_difference-2*response_slope*linear_correction-response_mixed*constant_correction)
    finite_correction = correction_response*(constant_correction+linear_correction*(first_virtuality+second_virtuality)+mixed_correction*first_virtuality*second_virtuality)
    interpolated_form_factor = shifted_convolution+finite_correction
    normalization_residual = sp.factor(shifted_normalization+finite_correction.subs({first_virtuality:0,second_virtuality:0})-anomaly)
    slope_residual = sp.simplify(shifted_slope+sp.diff(finite_correction,first_virtuality).subs({first_virtuality:0,second_virtuality:0})-original_slope)
    mixed_residual = sp.simplify(shifted_mixed+sp.diff(finite_correction,first_virtuality,second_virtuality).subs({first_virtuality:0,second_virtuality:0})-original_mixed)
    correction_single_hard = sp.limit(hard_scale*finite_correction.subs({first_virtuality:-hard_scale,second_virtuality:0}),hard_scale,sp.oo)
    correction_double_hard = sp.limit(hard_scale*finite_correction.subs({first_virtuality:-hard_ratio*hard_scale,second_virtuality:-hard_scale}),hard_scale,sp.oo)
    odd_derivative_map = sp.Matrix([original_slope,original_mixed]).jacobian([c1,c2,c5,c6,d3])
    odd_derivative_rank = odd_derivative_map.rank()
    shifted_hard_integrand = sp.limit(hard_scale*distribution_amplitude/
        (first_mass**2+momentum_fraction*hard_ratio*hard_scale+(1-momentum_fraction)*hard_scale),hard_scale,sp.oo)
    shifted_hard_coefficient = sp.factor(-2*decay_constant*sp.integrate(shifted_hard_integrand,(momentum_fraction,0,1))/3)
    hard_coefficient_residual = sp.factor(shifted_hard_coefficient-leading_twist_ray)
    # Subtraction variations produce an invertible triangular map to the two
    # low-energy derivatives, independently of the odd-coefficient values.
    subtraction_derivative_map = sp.Matrix([[1,0],[2*response_slope,1]])
    compensated_odd_variations = -subtraction_derivative_map.inv()*odd_derivative_map
    return locals()


def main():
    """Display calculation inputs, solved constraints and their interpretation."""
    result = derive()
    print('Chiral pion, canonical epsilon^(0123)=+1; F_pi(0,0)<0. GeV units.')
    print('Alda:2026nrd convention: f_pi=93 MeV=0.093 GeV.')
    pion_inputs = {result['decay_constant']: sp.Rational(93,1000), result['colors']: 3}
    show('Anomaly at N_C=3, f_pi=0.093 GeV [GeV^-1]',
         sp.N(result['anomaly'].subs(pion_inputs),12))
    show('Symmetric hard coefficient at f_pi=0.093 GeV [GeV]',
         sp.N(result['symmetric_limit'].subs(pion_inputs),12))
    print('Two vector poles per photon; polynomial truncation is an explicit approximation.')
    for label,name in [('Bose-symmetric two-pole ansatz','form_factor'),('Anomaly','anomaly'),
        ('Single-virtual constant','single_constant'),('Solved constraints','constraints'),
        ('Single-virtual coefficient before fixing h5','single_limit'),('Symmetric double-virtual coefficient','symmetric_limit'),
        ('Unequal hard ray: Q_squared F(-r Q_squared,-Q_squared)','ray_limit'),
        ('Completed derivative dF/ds1 at zero','slope'),('Completed derivative d2F/ds1 ds2 at zero','mixed'),
        ('Original eleven-operator chiral pion expression','original'),
        ('Original slope','original_slope'),('Original mixed derivative','original_mixed'),
        ('Conditional slope-preserving equation','slope_relation'),('Conditional mixed-derivative-preserving h2','h2_solution')]:
        value = result[name]
        show(label,sp.Dict(value) if isinstance(value,dict) else value)
    print('C_pi_single is retained symbolically above; the leading-twist pion value is derived below. h2 remains free.')
    print('Preserving both derivatives is an optional low-energy matching assumption, not an OPE constraint.')
    print('The unequal-ray coefficient has endpoint singularities; symmetric matching alone is insufficient.')
    print('Pion matching does not infer eta, eta-prime or ALP QCD coefficients; the separate eta/eta-prime outputs below use Bali:2021qem, while the ALP coefficient remains uninferred.')
    print('Gerardin:2019vio leading twist, leading alpha_s order, asymptotic pion DA.')
    print('Qbar_squared=(r+1) Q_squared; omega=(1-r)/(1+r); canonical negative scalar.')
    for label,name in [('Asymptotic pion distribution amplitude','distribution_amplitude'),
        ('Convolution in the source asymmetry convention','convolution'),
        ('Integrated dimensionless convolution in x*r+1-x convention','leading_twist_integral'),
        ('Q_squared F_pi(-r Q_squared,-Q_squared), leading twist','leading_twist_ray'),
        ('Two-pole minus leading-twist hard-ray coefficient','hard_ray_difference'),
        ('Leading twist r to 1','leading_twist_symmetric'),
        ('Leading twist r to zero','leading_twist_small_ratio'),
        ('Leading twist r to infinity','leading_twist_large_ratio'),
        ('Leading twist r times coefficient, r to infinity','leading_twist_rescaled_large_ratio'),
        ('h5 after leading-twist single-virtual matching','leading_twist_h5'),
        ('Two poles r to zero','two_pole_small_ratio'),
        ('Two poles r to infinity','two_pole_large_ratio')]:
        show(label,result[name])
    print('Alpha_s corrections, nonasymptotic DA moments and higher twist are omitted explicitly.')
    print('Endpoint limits are taken after the fixed-positive-r hard limit.')
    print('General meson rays assume the asymptotic quark DA for each flavor component.')
    print('LO alpha_s; singlet evolution and threshold matching are absorbed into C_eta_single and C_etap_single.')
    print('Bali:2021qem inputs: C_eta_single=-0.1605 +/- 0.0100 GeV; C_etap_single=-0.2305 +/- 0.0101 GeV.')
    for label,name in [('Universal normalized leading-twist ray Phi(r)','normalized_hard_ray'),
        ('Canonical Bali coefficient in octet-singlet decay constants','bali_coefficient'),
        ('Pion, eta, eta-prime Q_squared F hard coefficients','channel_hard_rays'),
        ('Eta hard ray with central asymptotic coefficient [GeV]','eta_hard_ray'),
        ('Eta-prime hard ray with central asymptotic coefficient [GeV]','eta_prime_hard_ray'),
        ('Eta symmetric coefficient [GeV]','eta_symmetric_coefficient'),
        ('Eta-prime symmetric coefficient [GeV]','eta_prime_symmetric_coefficient'),
        ('Eta symmetric uncertainty [GeV]','eta_symmetric_uncertainty'),
        ('Eta-prime symmetric uncertainty [GeV]','eta_prime_symmetric_uncertainty')]:
        show(label,result[name])
    for label,central,error in [('Eta','eta_symmetric_coefficient','eta_symmetric_uncertainty'),
                               ('Eta-prime','eta_prime_symmetric_coefficient','eta_prime_symmetric_uncertainty')]:
        print(f"{label} symmetric coefficient [GeV]: {sp.N(result[central],12)} +/- {sp.N(result[error],12)}")
    print('These are conditional leading-twist limits, not an all-ratio vector representation or finite-virtuality fit.')
    print('Diagnostic finite-virtuality interpolation: M1 below is an unfixed smoothing scale, not an identified physical pole.')
    for label,name in [('Shifted leading-twist convolution','shifted_convolution'),
        ('Shifted normalization','shifted_normalization'),('Shifted slope','shifted_slope'),
        ('Shifted mixed derivative','shifted_mixed'),('Rapidly falling correction response','correction_response'),
        ('Constant subtraction coefficient','constant_correction'),('Linear subtraction coefficient','linear_correction'),
        ('Mixed subtraction coefficient','mixed_correction'),('Finite correction','finite_correction'),
        ('Normalization residual','normalization_residual'),('Slope residual','slope_residual'),
        ('Mixed derivative residual','mixed_residual'),('Q_squared times correction, single virtual','correction_single_hard'),
        ('Q_squared times correction, fixed positive ratio','correction_double_hard'),
        ('Odd-coupling derivative map','odd_derivative_map'),('Odd-coupling derivative rank','odd_derivative_rank')]:
        show(label,result[name])
    show('Shifted convolution hard coefficient minus QCD coefficient',result['hard_coefficient_residual'])
    show('Linear/mixed subtraction derivative map',result['subtraction_derivative_map'])
    show('Subtraction variations compensating every odd-coefficient variation',result['compensated_odd_variations'])
    print('The shifted convolution retains the leading-twist hard-ray function: its mass shift is suppressed at fixed r>0.')
    print('The rational correction preserves single and double hard coefficients while fixing anomaly and two derivatives.')
    print('Its repeated poles and shifted continuum are interpolation assumptions; no physical spectral discontinuity is inferred.')
    print('The two low-energy derivatives measure at most two odd-coupling combinations when the interpolation is fixed.')
    print('Independent linear and mixed subtraction terms can absorb both combinations; marginalizing them removes this coefficient identification.')
    print('A physical dispersive construction still requires channel spectral densities, subtraction information and matching without duplicated resonance strength.')
    derive_sourced_dispersive(result)


def derive_sourced_dispersive(result):
    """Derive Hoferichter:2018kwz continuum and symbolic moment matching.

    All low-energy moments are independent unspecified inputs. The canonical
    scalar reverses the source sign; f_pi retains the Alda:2026nrd value.
    The profile below is algebraic, without a likelihood or invented data.
    """
    first,second = result['first_virtuality'],result['second_virtuality']
    pion_decay_constant = result['decay_constant']
    threshold = sp.Symbol('s_m',positive=True)
    spectral_variable,scaled_variable = sp.symbols('x_spectral t',positive=True)
    hard_ratio,hard_scale = result['hard_ratio'],result['hard_scale']
    asymptotic_integrand = -2*pion_decay_constant*first*second/((spectral_variable-first)**2*(spectral_variable-second)**2)
    asymptotic_component = sp.Integral(asymptotic_integrand,(spectral_variable,threshold,sp.oo))
    origin = {first:0,second:0}
    continuum_normalization = sp.integrate(asymptotic_integrand.subs(origin),(spectral_variable,threshold,sp.oo))
    continuum_slope = sp.integrate(sp.diff(asymptotic_integrand,first).subs(origin),(spectral_variable,threshold,sp.oo))
    continuum_mixed = sp.integrate(sp.diff(asymptotic_integrand,first,second).subs(origin),(spectral_variable,threshold,sp.oo))
    # Scale x=Q_squared*t before taking the fixed-positive-r limit. Integrate
    # SymPy's partial fractions to keep the logarithmic result explicit.
    scaled_integrand = sp.factor(hard_scale**2*asymptotic_integrand.subs(
        {first:-hard_ratio*hard_scale,second:-hard_scale,spectral_variable:hard_scale*scaled_variable}))
    primitive = sp.integrate(sp.apart(scaled_integrand,scaled_variable),scaled_variable)
    continuum_hard = sp.simplify(sp.limit(primitive,scaled_variable,sp.oo)-sp.limit(primitive,scaled_variable,0,dir='+'))
    low_normalization,low_slope,low_mixed,low_single_hard = sp.symbols('F_disp_0 F_disp_s F_disp_12 C_disp_single',real=True)
    effective_normalization = sp.Symbol('F_eff_0',real=True)
    effective_mass_squared = sp.Symbol('M_rhop_squared',positive=True)
    effective_pole = effective_normalization*effective_mass_squared**2/((effective_mass_squared-first)*(effective_mass_squared-second))
    effective_slope = sp.diff(effective_pole,first).subs(origin)
    effective_mixed = sp.diff(effective_pole,first,second).subs(origin)
    effective_single = sp.limit(hard_scale*effective_pole.subs({first:-hard_scale,second:0}),hard_scale,sp.oo)
    physical_anomaly = result['anomaly'].subs(result['colors'],3)
    single_hard = -2*pion_decay_constant
    normalization_solution = sp.solve(low_normalization+effective_normalization-physical_anomaly,effective_normalization)[0]
    mass_solution = sp.solve(low_single_hard+effective_single.subs(effective_normalization,normalization_solution)-single_hard,effective_mass_squared)[0]
    matched_effective = {effective_normalization:normalization_solution,effective_mass_squared:mass_solution}
    low_slope_solution = sp.factor(sp.solve(low_slope+effective_slope.subs(matched_effective,simultaneous=True)-result['original_slope'],low_slope)[0])
    low_mixed_solution = sp.factor(sp.solve(low_mixed+effective_mixed.subs(matched_effective,simultaneous=True)+continuum_mixed-result['original_mixed'],low_mixed)[0])
    profiled_slope = low_slope+effective_slope.subs(matched_effective,simultaneous=True)
    profiled_mixed = low_mixed+effective_mixed.subs(matched_effective,simultaneous=True)+continuum_mixed
    nuisance_profile_map = sp.Matrix([profiled_slope,profiled_mixed]).jacobian(
        [low_slope,low_mixed,low_normalization,low_single_hard,threshold]).applyfunc(sp.factor)
    spectral_x,spectral_y = sp.symbols('x y',positive=True)
    source_density = sp.Function('rho_disp')(spectral_x,spectral_y)
    spectral_integrand = -source_density/sp.pi**2*(1/((spectral_x-first)*(spectral_y-second))+
        1/((spectral_x-second)*(spectral_y-first)))
    low_moment_integrands = sp.Tuple(spectral_integrand.subs(origin),
        sp.diff(spectral_integrand,first).subs(origin),
        sp.diff(spectral_integrand,first,second).subs(origin),
        sp.limit(hard_scale*spectral_integrand.subs({first:-hard_scale,second:0}),hard_scale,sp.oo))
    print('Hoferichter:2018kwz Eqs. (4.1), (4.9), (4.10), canonical negative scalar; f_pi=0.093 GeV.')
    for label,expression in [('Asymptotic continuum','asymptotic_component'),
        ('Continuum anomaly contribution','continuum_normalization'),('Continuum first slope','continuum_slope'),
        ('Continuum mixed derivative','continuum_mixed'),('Continuum fixed-r hard coefficient','continuum_hard'),
        ('Continuum hard coefficient minus sourced leading twist',sp.simplify(continuum_hard-result['leading_twist_ray'])),
        ('Continuum with second photon real',asymptotic_integrand.subs(second,0)),
        ('Distinct effective pole','effective_pole'),('Effective-pole first slope','effective_slope'),
        ('Effective-pole mixed derivative','effective_mixed'),('Effective-pole single-hard coefficient','effective_single'),
        ('Anomaly fixes effective normalization conditionally','normalization_solution'),
        ('Single-hard coefficient fixes effective mass squared conditionally','mass_solution'),
        ('Unsubstituted low-energy slope relation',sp.Eq(low_slope,result['original_slope']-effective_slope)),
        ('Unsubstituted low-energy mixed relation',sp.Eq(low_mixed,result['original_mixed']-effective_mixed-continuum_mixed)),
        ('Remaining low-energy slope moment','low_slope_solution'),('Remaining low-energy mixed moment','low_mixed_solution'),
        ('Low-energy spectral integrands for normalization, slope, mixed, single-hard moments',low_moment_integrands),
        ('Symbolic nuisance sensitivity of slope and mixed derivative','nuisance_profile_map')]:
        show(label,locals().get(expression,expression) if isinstance(expression,str) else expression)
    print('No spectral function, F_pi^V or f1 is supplied here. Low-energy normalization, single-hard moment and threshold remain free.')
    print('The moment integrands use the source positive rho_disp with a canonical overall minus sign, integrated over its finite low-energy domain.')
    print('The effective mass solution requires nonzero effective normalization and positive mass squared; normalization zero is a separate degenerate case.')
    print('Independent low-energy slope and mixed moments give an identity submatrix in the nuisance derivative map.')
    print('Thus profiling these unspecified moments leaves both original slope combinations unconstrained; no statistical interval follows.')
    print('The continuum contributes zero singly virtual at every finite virtuality, despite its nonzero fixed-r hard endpoint limit: the limits do not commute.')


if __name__ == '__main__':
    main()
