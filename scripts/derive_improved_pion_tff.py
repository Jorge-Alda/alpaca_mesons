"""Two-pole pion TFF constraints, with the canonical signed scalar.

Invocation: .venv/bin/python scripts/derive_improved_pion_tff.py
Units: GeV; f_pi=0.093 GeV in the Alda:2026nrd convention; chiral pion, tree meromorphic approximation, zero widths.
Inputs: memory/rcht/odd-parity/radiative-tff-amplitude.md and
memory/extracts/Knecht:2001qf/pion-tff-short-distance.md.
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
    print('C_pi_single requires an independently sourced QCD value. h2 is not fixed by the leading limits.')
    print('Preserving both derivatives is an optional low-energy matching assumption, not an OPE constraint.')
    print('The unequal-ray coefficient has endpoint singularities; symmetric matching alone is insufficient.')
    print('No eta, eta-prime or ALP QCD coefficients are inferred from pion matching.')


if __name__ == '__main__':
    main()
