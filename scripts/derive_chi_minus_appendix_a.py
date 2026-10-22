"""Derive the mass building blocks from Alda:2026nrd Appendix A.

Run: .venv/bin/python scripts/derive_chi_minus_appendix_a.py
Tests: .venv/bin/python scripts/derive_chi_minus_appendix_a.py --test

Inputs: Appendix A equations (A.1), (A.5)--(A.7), as indexed in
memory/extracts/Alda:2026nrd/appendix-a-expansion.md. The defining equation
is chi_+-=u^dagger chi u^dagger +- u chi^dagger u, u=exp(i phi),
phi=Phi/(2 f_pi). Hermitian chi_S, chi_P satisfy chi=chi_S+i chi_P.
chi_0=2 B_0 m_q; retain one ALP and quadratic total field order.
Natural units: phi and couplings dimensionless; f_pi, f_a, a and Phi
have mass dimension one; chi has mass dimension two.

Alda uses Phi=sum_P lambda_P P, tr(lambda_i lambda_j)=2 delta_ij,
u=exp(i Phi/(2 f_pi)), f_pi=0.093 GeV. An outgoing state
substitutes Phi=P lambda^dagger. Its complex generator is not assumed
Hermitian. No field equation or meson mass approximation is imposed.
"""

import argparse
import unittest

import sympy as sp


def show(label, expression):
    """Render the calculated expression using SymPy's LaTeX printer."""
    print(f'{label}:\n{sp.latex(expression)}\n')


def derive():
    """Multiply the defining noncommutative exponential series explicitly."""
    field_order, alp, alp_scale, decay_constant = sp.symbols(
        'epsilon a f_a f_pi', real=True, nonzero=True)
    phi, mass, scalar_coupling, pseudoscalar_coupling = sp.symbols(
        'phi chi_0 c_q^S c_q^P', commutative=False)
    scalar_source = mass + field_order*sp.I*alp/alp_scale*(
        scalar_coupling*mass-mass*scalar_coupling)
    pseudoscalar_source = field_order*alp/alp_scale*(
        pseudoscalar_coupling*mass+mass*pseudoscalar_coupling)
    source = scalar_source+sp.I*pseudoscalar_source
    conjugate_source = scalar_source-sp.I*pseudoscalar_source
    exponential = sum((sp.I*field_order*phi)**power/sp.factorial(power)
                      for power in range(3))
    conjugate_exponential = exponential.subs(field_order,-field_order)
    first_product = sp.expand(conjugate_exponential*source*conjugate_exponential)
    second_product = sp.expand(exponential*conjugate_source*exponential)
    plus_orders = tuple(sp.expand(first_product+second_product).coeff(field_order,power)
                        for power in range(3))
    minus_orders = tuple(sp.expand(first_product-second_product).coeff(field_order,power)
                         for power in range(3))
    outgoing_generator = sp.Symbol('lambda_dagger',commutative=False)
    outgoing_chi_minus = sp.expand(minus_orders[1].subs(alp,0).subs(
        phi,outgoing_generator/(2*decay_constant)))
    outgoing_mass_source = sp.expand(sp.I*outgoing_chi_minus)
    alp_mass_source = sp.expand(sp.I*minus_orders[1].subs(phi,0)/alp)
    return locals()


class MassBuildingBlockTests(unittest.TestCase):
    """Check Appendix trigonometric expressions and independently signed limits."""

    def test_appendix_a_trigonometric_form(self):
        result = derive()
        phi, mass = result['phi'], result['mass']
        anticommutator = lambda first,second: first*second+second*first
        expected_plus = (2*result['scalar_source']
            -result['field_order']**2*anticommutator(phi,anticommutator(phi,mass))
            +2*result['field_order']*anticommutator(phi,result['pseudoscalar_source']))
        expected_minus = (2*sp.I*result['pseudoscalar_source']
            -2*sp.I*result['field_order']*anticommutator(phi,result['scalar_source']))
        for name, expected in [('plus_orders',expected_plus),('minus_orders',expected_minus)]:
            for power, coefficient in enumerate(result[name]):
                self.assertEqual(sp.expand(coefficient-sp.expand(expected).coeff(
                    result['field_order'],power)),0)

    def test_commuting_exact_limit(self):
        # Exact scalar exponentials provide an independent sign and normalization.
        phase, mass, pseudoscalar = sp.symbols('x m p', real=True)
        exact = sp.exp(-2*sp.I*phase)*(mass+sp.I*pseudoscalar)-sp.exp(
            2*sp.I*phase)*(mass-sp.I*pseudoscalar)
        self.assertEqual(sp.simplify(sp.diff(exact,phase).subs(phase,0)),-4*sp.I*mass)
        result=derive()
        calculated=result['outgoing_mass_source'].subs({result['mass']:mass,
            result['outgoing_generator']:1})
        self.assertEqual(sp.simplify(calculated-2*mass/result['decay_constant']),0)

    def test_alda_canonical_pion_and_rcht_conversion(self):
        # The neutral pion has lambda_3=diag(1,-1,0). Its kinetic
        # coefficient is Tr(lambda_3^2)/4=1/2 in Alda conventions.
        generator = sp.diag(1,-1,0)
        self.assertEqual(sp.trace(generator**2)/4,sp.Rational(1,2))
        result = derive()
        rcht_generator = result['outgoing_generator']/sp.sqrt(2)
        rcht_phi = rcht_generator/(sp.sqrt(2)*result['decay_constant'])
        inherited = sp.expand(result['minus_orders'][1].subs(
            result['alp'],0).subs(result['phi'],rcht_phi))
        self.assertEqual(sp.expand(inherited-result['outgoing_chi_minus']),0)

    def test_no_mixing_basis_preserves_meson_sign(self):
        result=derive()
        self.assertFalse(result['outgoing_mass_source'].has(result['scalar_coupling']))
        self.assertFalse(result['outgoing_mass_source'].has(result['pseudoscalar_coupling']))
        self.assertEqual(sp.expand(result['minus_orders'][2].subs(
            result['scalar_coupling'],0)),0)


def main():
    """Print each total-field coefficient before extracting source vertices."""
    result=derive()
    print(__doc__)
    show('u through quadratic field order',result['exponential'])
    show('chi_S',result['scalar_source'])
    show('chi_P',result['pseudoscalar_source'])
    for name in ('plus_orders','minus_orders'):
        for power,expression in enumerate(result[name]):
            show(f'{name}, total field degree {power}',expression)
    show('chi_- coefficient of outgoing P',result['outgoing_chi_minus'])
    show('i chi_- coefficient of outgoing P',result['outgoing_mass_source'])
    show('i chi_- coefficient of a',result['alp_mass_source'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test',action='store_true')
    arguments=parser.parse_args()
    if arguments.test:
        unittest.main(argv=['derive_chi_minus_appendix_a.py'])
    else:
        main()
