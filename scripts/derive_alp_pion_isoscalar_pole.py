#!/usr/bin/env python3
"""Derive pole and subtraction identities for the pion production wave.

Invocation: .venv/bin/python scripts/derive_alp_pion_isoscalar_pole.py
Natural units. The pole position is complex and away from the subtraction point.
The numerator and self-energy derivatives are analytic continuations to the
sheet of the pole. No numerical fit or assumption of a stable vector is made.
"""
import sympy as sp
import unittest


def source_residue():
    """Extract the coefficient of one neutral-vector response in N_j."""
    C, U, Y, T = (sp.MatrixSymbol(name, 3, 3) for name in ('C_j', 'U_j', 'Y_j', 'T_I'))
    r, y, mass_squared, p, e, d1, d2, d3 = sp.symbols('r_j y m_j^2 p_j e d_1 d_2 d_3')
    response, photon_residue = sp.symbols('W_I R_I')
    W = response * photon_residue * T
    axial_trace = sp.Trace((C * W + W * C) * U) / 2
    mass_trace = sp.Trace((C * W + W * C) * Y) / 2
    numerator = (p * 4 * sp.sqrt(2) * e * (((d1-d3)*mass_squared+d3*(r+y))*axial_trace-4*d2*mass_trace)).doit().expand()
    return numerator, sp.diff(numerator, response).expand()


def vppp_residue():
    """Coefficient of the isoscalar photon response in the indexed VPPP wave."""
    T, U, Y, plus, minus, chi = (sp.MatrixSymbol(name, 3, 3) for name in
                                ('T_I', 'U_a', 'Y_a', 'lambda_plus', 'lambda_minus', 'chi'))
    s, y, ma2, mpi2, MV, fpi, R = sp.symbols('s y m_a^2 m_pi^2 M_V f_pi R_I')
    g1, g2, g3, g4, g5 = sp.symbols('g_1 g_2 g_3 g_4 g_5')
    comm = lambda a, b: a*b-b*a
    anti = lambda a, b: a*b+b*a
    W = R*T
    first = ((g1-g3)*(ma2+4*mpi2-s-y)+2*g2*(y-ma2-s))*sp.Trace(W*U*comm(plus, minus))
    second = ((g1-g3)*(y-ma2-s)+2*g2*(2*mpi2-s))*sp.Trace(W*(plus*U*minus-minus*U*plus))
    fourth = g4*sp.Trace(anti(W, comm(plus, minus))*Y + anti(W, comm(plus, U))*anti(minus, chi) - anti(W, comm(minus, U))*anti(plus, chi))
    fifth = g5*sp.Trace((plus*W*minus-minus*W*plus)*Y + (plus*W*U-U*W*plus)*anti(minus, chi) - (minus*W*U-U*W*minus)*anti(plus, chi))
    return (2*sp.sqrt(2)/(MV*fpi**2)*(first+second+fourth+fifth)).doit().expand()


class DerivationTests(unittest.TestCase):
    def test_source_is_linear(self):
        numerator, coefficient = source_residue()
        response = sp.Symbol('W_I')
        self.assertEqual(sp.expand((numerator - response * coefficient).doit()), 0)

    def test_two_subtractions(self):
        y, pole, y0 = sp.symbols('y s_I y_0')
        value = 1/(y-pole)
        remainder = sp.factor(value - value.subs(y, y0) - (y-y0)*sp.diff(value, y).subs(y, y0))
        self.assertEqual(sp.simplify(remainder - (y-y0)**2/((y-pole)*(y0-pole)**2)), 0)


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    y, pole, subtraction_point = sp.symbols('y s_I y_0')
    mass_squared, pole_numerator, pole_numerator_derivative = sp.symbols(
        'M_I^2 N_I N_I_prime')
    self_energy, self_energy_derivative, self_energy_second_derivative = sp.symbols(
        'Delta_I Delta_I_prime Delta_I_doubleprime')
    electric_charge = sp.symbols('e', nonzero=True)
    # Linear expansion is sufficient for the Laurent residue and regular term.
    numerator = pole_numerator + pole_numerator_derivative * (y - pole)
    denominator = (mass_squared - y - self_energy - self_energy_derivative * (y - pole)
                   - self_energy_second_derivative * (y - pole)**2 / 2)
    denominator = denominator.subs(self_energy, mass_squared - pole)
    wave = numerator / (electric_charge * denominator)
    residue = sp.simplify(sp.limit((y - pole) * wave, y, pole))
    regular = sp.simplify(sp.limit(wave - residue / (y - pole), y, pole))
    # A subtraction at y_0 rearranges the same pole into a subtracted pole
    # and a constant. This identity preserves the Laurent residue.
    pole_term = residue / (y - pole)
    subtraction_constant = sp.simplify(pole_term.subs(y, subtraction_point))
    subtracted_pole = sp.factor(pole_term - subtraction_constant)
    show('Pole-expanded numerator', numerator)
    show('Pole-expanded inverse propagator', denominator)
    show('Residue in the y minus s_I convention', residue)
    show('Regular term including the self-energy second derivative', regular)
    show('Pole contribution to the subtraction value', subtraction_constant)
    show('Once-subtracted pole term', subtracted_pole)
    indexed_numerator, indexed_coefficient = source_residue()
    show('Indexed VVP numerator with a single isoscalar response', indexed_numerator)
    show('Extracted isoscalar numerator coefficient', indexed_coefficient)
    show('Extracted VPPP isoscalar coefficient in the charge-stripped wave', vppp_residue())
    # Substitution of the indexed table yields each actual channel numerator;
    # all traces and sources remain matrices, with no commutativity assumption.
    C, U, Y = (sp.MatrixSymbol(name, 3, 3) for name in ('C_j', 'U_j', 'Y_j'))
    plus, minus, ca, Ua, Ya, Uplus, Uminus, Yplus, Yminus = (
        sp.MatrixSymbol(name, 3, 3) for name in
        ('lambda_plus', 'lambda_minus', 'c_q^A', 'U_a', 'Y_a', 'U_pi_plus', 'U_pi_minus', 'Y_pi_plus', 'Y_pi_minus'))
    s, t, u, ma2, mpi2, GV, fpi, fa, Mrho2 = sp.symbols('s t u m_a^2 m_pi^2 G_V f_pi f_a M_rho^2')
    rj, mj2, pj = sp.symbols('r_j m_j^2 p_j')
    channel_data = [(s, ma2, -GV/(sp.sqrt(2)*fpi**2), plus*minus-minus*plus, Ua, Ya),
                    (t, mpi2, sp.sqrt(2)*GV/(fpi*fa), plus*ca-ca*plus, Uminus, Yminus),
                    (u, mpi2, -sp.sqrt(2)*GV/(fpi*fa), minus*ca-ca*minus, Uplus, Yplus)]
    channel_terms = [indexed_coefficient.xreplace({rj:r, mj2:m, pj:p, C:c, U:ax, Y:ps})/(r-Mrho2)
                     for r,m,p,c,ax,ps in channel_data]
    z = sp.Symbol('z')
    t_z, u_z = sp.Function('t')(z), sp.Function('u')(z)
    projection = sp.Rational(3,4)*sp.Integral((1-z**2)*sum(channel_terms).subs({t:t_z,u:u_z})/electric_charge, (z,-1,1)) + vppp_residue()
    show('Isoscalar projected coefficient including VVP and VPPP', projection)
    show('Tree-matched pole residue with continued inverse propagator', -projection.subs(y,pole)/(1+self_energy_derivative))
    # A twice-subtracted dispersion relation leaves the value and derivative
    # at y_0 as independent subtraction functions of the dipion invariant.
    twice_subtracted_pole = sp.factor(pole_term - pole_term.subs(y, subtraction_point)
                                    - (y-subtraction_point)*sp.diff(pole_term, y).subs(y, subtraction_point))
    show('Twice-subtracted pole contribution', twice_subtracted_pole)
    assert sp.simplify(subtraction_constant + subtracted_pole - pole_term) == 0
    assert sp.simplify(sp.limit((y-pole)*subtracted_pole, y, pole)-residue) == 0
    # Isoscalar photon pole projected from the indexed s,t,u numerators.
    # N_{j,I} denotes the coefficient of the isoscalar photon pole in the
    # existing N_j; it includes R_I T_I and is not a new model input.
    z, s, t, u, rho_mass_squared = sp.symbols('z s t u M_rho^2')
    N_s_I, N_t_I, N_u_I, local_I = sp.symbols('N_s_I N_t_I N_u_I F_VPPP_I')
    projected_pole_numerator = sp.Rational(3, 4) * sp.Integral(
        (1-z**2) * (N_s_I/(s-rho_mass_squared)
                       + N_t_I/(t-rho_mass_squared)
                       + N_u_I/(u-rho_mass_squared) + electric_charge*local_I),
        (z, -1, 1))
    show('Isoscalar pole numerator from the indexed production projection',
         projected_pole_numerator)
    # Elastic unitarity with a left-hand contribution: solving the boundary
    # relation explicitly establishes when Watson's phase applies.
    phase = sp.symbols('delta_1^1', real=True)
    left_hand = sp.symbols('f_hat')
    upper = sp.symbols('f_plus')
    lower = sp.symbols('f_minus')
    solution = sp.solve(sp.Eq(upper - lower,
        2*sp.I*sp.sin(phase)*sp.exp(-sp.I*phase)*(upper + left_hand)), upper)[0]
    show('Upper-rim wave from elastic unitarity with crossed contribution', solution)
    show('Homogeneous upper-to-lower ratio', sp.simplify(solution.subs(left_hand, 0)/lower))
    # Divide the boundary equation by Omega_+ with Omega_+=exp(2i delta) Omega_-.
    omnes_abs = sp.symbols('Abs_Omega', positive=True)
    omnes_plus = omnes_abs * sp.exp(sp.I*phase)
    omnes_minus = omnes_abs * sp.exp(-sp.I*phase)
    quotient_jump = sp.simplify(solution/omnes_plus-lower/omnes_minus)
    show('Discontinuity of the single-variable wave divided by Omnes', quotient_jump)
    x, s = sp.symbols('x s')
    n = sp.symbols('n', integer=True, positive=True)
    omnes = sp.Function('Omega')
    crossed = sp.Function('f_hat')
    delta = sp.Function('delta_1^1')
    polynomial = sp.Function('P_n_minus_1')(s, y)
    threshold = sp.Symbol('s_threshold', positive=True)
    subtracted_wave = omnes(s)*(polynomial + s**n/sp.pi*sp.Integral(
        sp.sin(delta(x))*crossed(x,y)/(sp.Abs(omnes(x))*x**n*(x-s)), (x,threshold,sp.oo)))
    show('n-subtracted Omnes representation (boundary i0 understood)', subtracted_wave)
    production, decay, pole_position = sp.symbols('A_gamma_I A_I_pipi s_I')
    show('Pole-factorized current amplitude on the continued sheet', production*decay/(y-pole_position))
    print('Assumptions: simple pole, analytic numerator/self-energy near pole;')
    print('elastic unitarity equation shown with stated phase convention.')


if __name__ == '__main__':
    import sys
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        main()
