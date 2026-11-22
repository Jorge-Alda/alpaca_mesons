"""Symbolic cancellations in the dressed-frame loop counting.

Run: .venv/bin/python scripts/derive_alp_loop_cancellations.py
All masses and couplings are symbolic. Use consistent natural units.
The class contributions are formal signed coefficients, including loop weights.
"""

import unittest

import sympy as sp


def derive():
    B_0, m_0_squared, m_a_squared, N_f = sp.symbols(
        'B_0 m_0^2 m_a^2 N_f', nonzero=True, real=True
    )
    epsilon = sp.symbols('epsilon', real=True)
    m_q = sp.MatrixSymbol('m_q', 3, 3)
    scalar_first = sp.MatrixSymbol(r'\widehat{S}_1', 3, 3)
    inverse_on_C = sp.MatrixSymbol(
        r'\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})', 3, 3
    )
    theta_first, theta_second = sp.symbols(
        r'\widehat{\Theta}_1 \widehat{\Theta}_2', real=True
    )
    identity = sp.Identity(3)

    # In the commuting example Mhat(x)=m_q exp(ix), so P1=m_q and S2=-m_q.
    # M_(ma²)(I/2) is computed before solving the source equation for Theta1.
    inverse_at_special_source = identity / 2
    operator_on_special_source = (
        B_0 * (m_q * inverse_at_special_source + inverse_at_special_source * m_q)
        + m_0_squared / N_f * sp.Trace(inverse_at_special_source) * identity
        - m_a_squared * inverse_at_special_source
    ).doit()
    theta_for_special_source = sp.solve(
        sp.Eq(-m_0_squared / N_f * theta_first,
              3 * m_0_squared / (2 * N_f) - m_a_squared / 2),
        theta_first,
    )[0]
    mixing_special_source = B_0 * m_q - m_0_squared / N_f * theta_for_special_source * identity
    inverse_equation_residual = sp.simplify((operator_on_special_source - mixing_special_source).doit())
    scalar_contraction = (
        -m_q + 4 * inverse_at_special_source * m_q
        - 4 * (identity / 4) * m_q
    ).doit()

    # General second-order source choice, retaining noncommuting matrices.
    pseudoscalar_second = (
        2 * (inverse_on_C * scalar_first + scalar_first * inverse_on_C)
        + m_0_squared / (N_f * B_0) * theta_second * identity
    )
    cubic_contraction = (
        B_0 * (pseudoscalar_second - 2 * (
            inverse_on_C * scalar_first + scalar_first * inverse_on_C
        )) - m_0_squared / N_f * theta_second * identity
    ).doit()
    pseudoscalar_first = m_0_squared / (N_f * B_0) * theta_first * identity
    mixing_source = (B_0 * pseudoscalar_first - m_0_squared / N_f * theta_first * identity).doit()

    # Formal class coefficients include their respective loop weights.
    class_2, class_3, class_5, class_6 = sp.symbols('class_2 class_3 class_5 class_6')
    normalization = 1 + 2 * epsilon**2 * sp.Trace(inverse_on_C**2)
    leading_sum = class_2 + class_3 + class_5 + class_6
    condition = sp.solve(leading_sum, class_5)[0]
    finite_sum = sp.factor(((class_2 + class_3 + class_6) / normalization + class_5).subs(class_5, condition))
    return {
        'Theta1 for the inverse equation': theta_for_special_source,
        'Inverse equation residual': inverse_equation_residual,
        'Class 2 scalar contraction': scalar_contraction,
        'Class 6 cubic contraction': cubic_contraction,
        'Class 5 mixing source': mixing_source,
        'Condition on the class 5 coefficient': condition,
        'Finite-normalization sum under the condition': finite_sum,
    }


def main():
    for label, expression in derive().items():
        print(f'{label}:\n{sp.latex(expression)}\n')


class SymbolicCancellationsTest(unittest.TestCase):
    def test_invariant_contractions(self):
        results = derive()
        for label in ('Inverse equation residual', 'Class 2 scalar contraction',
                      'Class 6 cubic contraction', 'Class 5 mixing source'):
            self.assertEqual(results[label].as_explicit().applyfunc(sp.simplify), sp.zeros(3))

    def test_finite_normalization(self):
        results = derive()
        finite_sum = results['Finite-normalization sum under the condition']
        self.assertNotEqual(finite_sum, 0)
        epsilon = next(symbol for symbol in finite_sum.free_symbols if symbol.name == 'epsilon')
        self.assertEqual(finite_sum.subs(epsilon, 0), 0)


if __name__ == '__main__':
    main()
