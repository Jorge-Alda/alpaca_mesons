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
This is a diagnostic completion, not a full arbitrary-ratio QCD match.
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


if __name__ == '__main__':
    main()
