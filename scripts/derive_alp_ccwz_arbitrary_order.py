"""Derive CCWZ expansions at arbitrary fixed ALP and meson orders.

Run .venv/bin/python scripts/derive_alp_ccwz_arbitrary_order.py --order 3 --mesons 3
Run .venv/bin/python scripts/derive_alp_ccwz_arbitrary_order.py --test

Natural units. x=a/f_a, u=exp(i phi), phi=Phi/(2 f_pi). A commuting
variable t counts meson fields, including a differentiated meson field.
The symbol dx denotes partial_mu x and counts as one inverse power of f_a.
All flavor symbols are noncommuting. Taylor coefficients C_eta,m and
X_eta-zeta,m are derivatives at x=0, so their series include 1/m!.
X_LR,m is the Hermitian conjugate of X_RL,m in the physical source; it
is kept as a distinct noncommuting symbol during the calculation.

The script derives finite exponential series and differentiates each matrix
power by the product rule. SymPy expands and extracts all coefficients.
Unit tests compare these direct definitions with nested commutator and
anticommutator formulas; no matrices are assumed to commute.
"""

import argparse
import unittest

import sympy as sp


ALP_VARIABLE, MESON_VARIABLE = sp.symbols("x t", real=True)
ALP_DERIVATIVE = sp.Symbol("partial_mu_x", real=True)
PHI, PHI_DERIVATIVE = sp.symbols("phi partial_mu_phi", commutative=False)
B_ZERO = sp.Symbol("B_0", real=True)


def exponential_series(sign, meson_order):
    """Expand exp(sign*i*t*phi) with SymPy's ordinary series coefficients."""
    return sp.Add(*(sign**degree * sp.I**degree * MESON_VARIABLE**degree
                    * PHI**degree / sp.factorial(degree)
                    for degree in range(meson_order + 1)))


def differentiated_exponential(sign, meson_order):
    """Apply the product rule to each noncommuting power of phi."""
    return sp.Add(*(sign**degree * sp.I**degree * MESON_VARIABLE**degree
                    / sp.factorial(degree)
                    * sp.Add(*(PHI**position * PHI_DERIVATIVE
                               * PHI**(degree - position - 1)
                               for position in range(degree)))
                    for degree in range(1, meson_order + 1)))


def source_coefficients(label, order):
    """Create the independent noncommuting flavor coefficients."""
    return tuple(sp.Symbol(f"{label}_{degree}", commutative=False)
                 for degree in range(order + 1))


def taylor_series(coefficients):
    """Build the Taylor series with derivatives at zero as coefficients."""
    return sp.Add(*(coefficient * ALP_VARIABLE**degree / sp.factorial(degree)
                    for degree, coefficient in enumerate(coefficients)))


def derive_blocks(alp_order, meson_order):
    """Derive the building blocks from their defining matrix products.

    Current series include C_eta,alp_order for explicit bookkeeping, even
    though its dx prefactor puts it at inverse-f_a order alp_order+1.
    """
    right = source_coefficients("C_R", alp_order)
    left = source_coefficients("C_L", alp_order)
    mass_right_left = source_coefficients("X_RL", alp_order)
    mass_left_right = source_coefficients("X_LR", alp_order)
    representative = exponential_series(1, meson_order)
    representative_adjoint = exponential_series(-1, meson_order)
    derivative = differentiated_exponential(1, meson_order)
    adjoint_derivative = differentiated_exponential(-1, meson_order)
    right_connection = ALP_DERIVATIVE * taylor_series(right)
    left_connection = ALP_DERIVATIVE * taylor_series(left)
    right_product = (representative_adjoint * derivative
                     - sp.I * representative_adjoint * right_connection
                     * representative)
    left_product = (representative * adjoint_derivative
                    - sp.I * representative * left_connection
                    * representative_adjoint)
    right_mass_product = (2 * B_ZERO * representative_adjoint
                          * taylor_series(mass_right_left)
                          * representative_adjoint)
    left_mass_product = (2 * B_ZERO * representative
                         * taylor_series(mass_left_right) * representative)
    return {
        "u_mu": sp.expand(sp.I * (right_product - left_product)),
        "Gamma_mu": sp.expand((right_product + left_product) / 2),
        "chi_plus": sp.expand(right_mass_product + left_mass_product),
        "chi_minus": sp.expand(right_mass_product - left_mass_product),
    }


def fixed_coefficient(expression, alp_degree, meson_degree):
    """Extract the coefficient of x^m t^k without Taylor factorials."""
    return sp.expand(expression).coeff(ALP_VARIABLE, alp_degree).coeff(
        MESON_VARIABLE, meson_degree)


def nested_product(expression, order, anticommutator=False):
    """Expand repeated [phi,...] or {phi,...} with SymPy products."""
    for _ in range(order):
        expression = PHI * expression + (1 if anticommutator else -1) * expression * PHI
    return sp.expand(expression)


def show(label, expression):
    """Render each derivation result with SymPy's standard LaTeX printer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


class CCWZExpansionTests(unittest.TestCase):
    """Compare direct matrix definitions with independent operator series."""

    @classmethod
    def setUpClass(cls):
        cls.blocks = derive_blocks(3, 5)

    def test_current_sources_all_retained_coefficients(self):
        for alp_degree in range(4):
            current_right = sp.Symbol(f"C_R_{alp_degree}", commutative=False)
            current_left = sp.Symbol(f"C_L_{alp_degree}", commutative=False)
            for meson_degree in range(6):
                right_rotation = ((-sp.I)**meson_degree
                                  * nested_product(current_right, meson_degree))
                left_rotation = (sp.I**meson_degree
                                 * nested_product(current_left, meson_degree))
                denominator = sp.factorial(alp_degree) * sp.factorial(meson_degree)
                expected_u = ALP_DERIVATIVE * (right_rotation - left_rotation) / denominator
                expected_gamma = -sp.I * ALP_DERIVATIVE * (right_rotation + left_rotation) / (2 * denominator)
                for name, expected in (("u_mu", expected_u), ("Gamma_mu", expected_gamma)):
                    actual = fixed_coefficient(self.blocks[name], alp_degree, meson_degree)
                    actual = sp.expand(actual).coeff(ALP_DERIVATIVE) * ALP_DERIVATIVE
                    self.assertEqual(sp.expand(actual - expected), 0)

    def test_mass_sources_all_retained_coefficients(self):
        for alp_degree in range(4):
            mass_rl = sp.Symbol(f"X_RL_{alp_degree}", commutative=False)
            mass_lr = sp.Symbol(f"X_LR_{alp_degree}", commutative=False)
            for meson_degree in range(6):
                right_rotation = (-sp.I)**meson_degree * nested_product(mass_rl, meson_degree, True)
                left_rotation = sp.I**meson_degree * nested_product(mass_lr, meson_degree, True)
                denominator = sp.factorial(alp_degree) * sp.factorial(meson_degree)
                for name, sign in (("chi_plus", 1), ("chi_minus", -1)):
                    expected = 2 * B_ZERO * (right_rotation + sign * left_rotation) / denominator
                    actual = fixed_coefficient(self.blocks[name], alp_degree, meson_degree)
                    self.assertEqual(sp.expand(actual - expected), 0)

    def test_source_free_vielbein_and_connection(self):
        for meson_degree in range(1, 6):
            nested_derivative = nested_product(PHI_DERIVATIVE, meson_degree - 1)
            expected_u = (2 * (-1)**((meson_degree + 1) // 2) * nested_derivative
                          / sp.factorial(meson_degree) if meson_degree % 2 else 0)
            expected_gamma = (-sp.I**meson_degree * nested_derivative
                              / sp.factorial(meson_degree) if meson_degree % 2 == 0 else 0)
            for name, expected in (("u_mu", expected_u), ("Gamma_mu", expected_gamma)):
                actual = fixed_coefficient(self.blocks[name], 0, meson_degree).subs(ALP_DERIVATIVE, 0)
                self.assertEqual(sp.expand(actual - expected), 0)

    def test_alda_vector_axial_and_mass_decomposition(self):
        """Recover Appendix A's sine/cosine series without meson commutativity."""
        vector, axial, scalar, pseudoscalar = sp.symbols(
            "c_q_V c_q_A chi_S chi_P", commutative=False)
        replacements = {
            sp.Symbol("C_R_0", commutative=False): vector + axial,
            sp.Symbol("C_L_0", commutative=False): vector - axial,
            sp.Symbol("X_RL_1", commutative=False): (scalar + sp.I * pseudoscalar) / (2 * B_ZERO),
            sp.Symbol("X_LR_1", commutative=False): (scalar - sp.I * pseudoscalar) / (2 * B_ZERO),
        }
        for degree in range(6):
            if degree % 2 == 0:
                factor = 2 * (-1)**(degree // 2) / sp.factorial(degree)
                expected_u = factor * nested_product(axial, degree)
                expected_plus = factor * nested_product(scalar, degree, True)
                expected_minus = sp.I * factor * nested_product(pseudoscalar, degree, True)
            else:
                factor = 2 * (-1)**((degree - 1) // 2) / sp.factorial(degree)
                expected_u = -sp.I * factor * nested_product(vector, degree)
                expected_plus = factor * nested_product(pseudoscalar, degree, True)
                expected_minus = -sp.I * factor * nested_product(scalar, degree, True)
            actual_u = fixed_coefficient(self.blocks["u_mu"], 0, degree).coeff(ALP_DERIVATIVE).subs(replacements)
            self.assertEqual(sp.expand(actual_u - expected_u), 0)
            for name, expected in (("chi_plus", expected_plus), ("chi_minus", expected_minus)):
                actual = fixed_coefficient(self.blocks[name], 1, degree).subs(replacements)
                self.assertEqual(sp.expand(actual - expected), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=2, help="Highest inverse-f_a order retained")
    parser.add_argument("--mesons", type=int, default=3, help="Highest number of pseudoscalar mesons")
    parser.add_argument("--test", action="store_true", help="Run the embedded unit tests")
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=[__file__])
    if arguments.order < 0 or arguments.mesons < 0:
        parser.error("orders must be nonnegative")
    blocks = derive_blocks(arguments.order, arguments.mesons)
    print("Natural units; x=a/f_a; dx=partial_mu a/f_a; phi=Phi/(2 f_pi).")
    print("All flavor products are noncommuting. Coefficients are derivatives at x=0.")
    print(f"Retain inverse-f_a orders 0..{arguments.order} and meson degrees 0..{arguments.mesons}.")
    show("u=exp(i t phi), expanded", exponential_series(1, arguments.mesons))
    show("partial_mu u from the product rule", differentiated_exponential(1, arguments.mesons))
    show("C_R(x), input through requested Taylor order", taylor_series(source_coefficients("C_R", arguments.order)))
    show("C_L(x), input through requested Taylor order", taylor_series(source_coefficients("C_L", arguments.order)))
    show("X_RL(x), input through requested Taylor order", taylor_series(source_coefficients("X_RL", arguments.order)))
    show("X_LR(x), input through requested Taylor order", taylor_series(source_coefficients("X_LR", arguments.order)))
    for name, expression in blocks.items():
        for inverse_fa_order in range(arguments.order + 1):
            for meson_degree in range(arguments.mesons + 1):
                if name in ("u_mu", "Gamma_mu"):
                    if inverse_fa_order == 0:
                        coefficient = fixed_coefficient(expression, 0, meson_degree).subs(ALP_DERIVATIVE, 0)
                        prefactor = sp.Integer(1)
                    else:
                        coefficient = fixed_coefficient(expression, inverse_fa_order - 1, meson_degree)
                        coefficient = sp.expand(coefficient).coeff(ALP_DERIVATIVE)
                        prefactor = ALP_VARIABLE**(inverse_fa_order - 1) * ALP_DERIVATIVE
                else:
                    coefficient = fixed_coefficient(expression, inverse_fa_order, meson_degree)
                    prefactor = ALP_VARIABLE**inverse_fa_order
                if coefficient != 0:
                    show(f"{name}, inverse-f_a order {inverse_fa_order}, {meson_degree} mesons",
                         prefactor * coefficient)
    print(f"C_eta,{arguments.order} is present in the input but first contributes at inverse-f_a order {arguments.order + 1}.")


if __name__ == "__main__":
    main()
