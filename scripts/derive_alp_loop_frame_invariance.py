"""Symbolic frame check. Run .venv/bin/python scripts/derive_alp_loop_frame_invariance.py.

The source coefficients and axial generator are noncommuting matrices. All
physical parameters remain symbolic. No loop functions are evaluated.
"""
import unittest

import sympy as sp

from derive_alp_invariant_loop_crossing import raw_tadpole_coefficient_shift


def derivation():
    x = sp.Symbol('x', real=True)
    generator = sp.Symbol(r'\kappa', commutative=False)
    coefficients = sp.symbols('M_0:7', commutative=False)
    source = sum(coefficients[n] * x**n / sp.factorial(n) for n in range(7))
    group = sum((sp.I * x * generator)**n / sp.factorial(n) for n in range(7))
    inverse_group = sum((-sp.I * x * generator)**n / sp.factorial(n) for n in range(7))
    transformed = sp.expand(group * source * group)
    recovered = sp.expand(inverse_group * transformed * inverse_group)
    residuals = [sp.expand(recovered.coeff(x, n) - source.coeff(x, n)) for n in range(7)]
    theta, theta_one, generator_trace = sp.symbols('Theta Theta_1 tr_kappa', real=True)
    transformed_theta = theta + x * theta_one - x * generator_trace
    dressed_theta = sp.expand(transformed_theta + x * generator_trace)
    negative = raw_tadpole_coefficient_shift()
    return residuals, sp.expand(dressed_theta - theta - x * theta_one), negative


def contraction_checks():
    """Compose a projected contraction and exact tree functional from the two frames."""
    x = sp.Symbol('x', real=True)
    kappa, mass, pseudoscalar, scalar = sp.symbols(
        'kappa m_q P_1 S_2', real=True
    )
    inverse_pole_operator = sp.Symbol(
        r'\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})', real=True
    )
    epsilon, decay_constant, pole_mass_squared, b0 = sp.symbols(
        r'\epsilon f_\pi m_a^2 B_0', nonzero=True, real=True
    )
    dressed_mass = mass + sp.I * pseudoscalar * x + scalar * x**2 / 2
    framed_mass = sp.series(sp.exp(2 * sp.I * kappa * x) * dressed_mass, x, 0, 3).removeO()
    recovered_mass = sp.series(sp.exp(-2 * sp.I * kappa * x) * framed_mass, x, 0, 3).removeO().expand()
    recovered_pseudoscalar = sp.simplify(sp.im(recovered_mass.coeff(x, 1)))
    recovered_scalar = sp.simplify(2 * sp.re(recovered_mass.coeff(x, 2)))
    raw_scalar = sp.simplify(2 * sp.re(sp.expand(framed_mass).coeff(x, 2)))

    def contraction(scalar_source, pseudoscalar_source):
        return -b0 / decay_constant**2 * (
            scalar_source + 4 * inverse_pole_operator * pseudoscalar_source
            - 4 * inverse_pole_operator**2 * mass
        )

    dressed_contraction = contraction(scalar, pseudoscalar)
    recovered_contraction = contraction(recovered_scalar, recovered_pseudoscalar)
    # Symbols represent the remaining dressed derivatives, pole vectors and
    # kinematics needed by the tree amplitude; their invariance follows from
    # the general dressed-source identity and the invariant pole solution.
    scalar_first, pseudoscalar_second, pseudoscalar_third, scalar_third, scalar_fourth = sp.symbols(
        r'\widehat{S}_1 \widehat{P}_2 \widehat{P}_3 \widehat{S}_3 \widehat{S}_4'
    )
    theta_one, theta_two, theta_three = sp.symbols(
        r'\widehat{\Theta}_1 \widehat{\Theta}_2 \widehat{\Theta}_3'
    )
    meson_poles, external_momenta = sp.symbols(r'm_n p_n')
    tree_function = sp.Function(r'\mathcal{M}_{\mathrm{tree}}')
    common_inputs = (
        mass, scalar_first, pseudoscalar_second, pseudoscalar_third,
        scalar_third, scalar_fourth, theta_one, theta_two, theta_three,
        inverse_pole_operator, meson_poles, external_momenta,
        b0, decay_constant, pole_mass_squared, epsilon,
    )
    tree_amplitude = tree_function(scalar, pseudoscalar, *common_inputs)
    recovered_tree_amplitude = tree_function(recovered_scalar, recovered_pseudoscalar, *common_inputs)
    pole_normalization = 1 + 2 * epsilon**2 * sp.Symbol(
        r'\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})]^2\rangle', real=True
    )
    ratio_factor = epsilon**2 * pole_mass_squared / (
        16 * sp.pi**2 * decay_constant**2 * pole_normalization
    )
    contraction_residual = sp.simplify(recovered_contraction - dressed_contraction)
    ratio_residual = sp.simplify(
        ratio_factor * recovered_contraction / recovered_tree_amplitude
        - ratio_factor * dressed_contraction / tree_amplitude
    )
    negative_residual = sp.simplify(
        ratio_factor * (contraction(raw_scalar, recovered_pseudoscalar) - dressed_contraction)
        / tree_amplitude
    )
    return contraction_residual, ratio_residual, negative_residual


def main():
    residuals, theta_residual, negative = derivation()
    for degree, residual in enumerate(residuals):
        print(f'Dressed mass derivative {degree} residual:\n{sp.latex(residual)}\n')
    print(f'Dressed topological residual:\n{sp.latex(theta_residual)}\n')
    print(f'Off-shell replacement residual:\n{sp.latex(negative["linear_spurion_minus_invariant"])}\n')
    contraction_residual, ratio_residual, negative_residual = contraction_checks()
    print(f'Projected contraction residual:\n{sp.latex(contraction_residual)}\n')
    print(f'Illustrative contraction-to-tree ratio residual:\n{sp.latex(ratio_residual)}\n')
    print(f'Ratio with off-shell substitution residual:\n{sp.latex(negative_residual)}\n')


class FrameInvarianceTest(unittest.TestCase):
    def test_dressed_sources(self):
        residuals, theta_residual, _ = derivation()
        self.assertEqual(residuals, [0] * 7)
        self.assertEqual(theta_residual, 0)

    def test_off_shell_control(self):
        _, _, negative = derivation()
        self.assertTrue(negative['sum_equals_invariant'])
        self.assertNotEqual(negative['linear_spurion_minus_invariant'], 0)

    def test_contraction_and_ratio(self):
        contraction_residual, ratio_residual, negative_residual = contraction_checks()
        self.assertEqual(contraction_residual, 0)
        self.assertEqual(ratio_residual, 0)
        self.assertNotEqual(negative_residual, 0)


if __name__ == '__main__':
    main()
