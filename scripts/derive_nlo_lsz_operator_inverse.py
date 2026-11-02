"""Derive the local NLO extension of Alda:2026nrd's meson operator inverse.

Run .venv/bin/python scripts/derive_nlo_lsz_operator_inverse.py [--test].
Natural units; general Hermitian U(3) matrices; local LECs retained linearly.
An exact inverse of that truncated operator is distinguished from its strict
first-order expansion. Loops and O(1/f_a**2) effects are omitted.
"""
import sys
import unittest

import sympy as sp

from derive_nlo_lsz_source_projection import derive as derive_sources
from derive_nlo_lsz_source_projection import expand_standard, trace_gradient


def derive():
    sources = derive_sources()
    chi = sources["chi"]
    flavor_matrix = sources["phi"]
    f_pi = sources["f_pi"]
    momentum_squared = sp.Symbol("w", real=True)
    couplings = {name: sp.Symbol(name, real=True) for name in sources["meson_kinetic"] if name != "LO"}
    couplings.update({name: sp.Symbol(name, real=True) for name in sources["meson_potential"] if name != "LO"})
    kinetic_operators = {
        name: expand_standard(2 * trace_gradient(expression, sources["derivative_phi"]).subs(sources["derivative_phi"], flavor_matrix))
        for name, expression in sources["meson_kinetic"].items()
    }
    mass_operators = {
        name: expand_standard(-2 * trace_gradient(expression, flavor_matrix))
        for name, expression in sources["meson_potential"].items()
    }
    total_kinetic_operator = expand_standard(
        kinetic_operators["LO"] + sum(
            (couplings[name] * expression for name, expression in kinetic_operators.items() if name != "LO"),
            sp.ZeroMatrix(3, 3),
        )
    )
    total_mass_operator = expand_standard(
        mass_operators["LO"] + sum(
            (couplings[name] * expression for name, expression in mass_operators.items() if name != "LO"),
            sp.ZeroMatrix(3, 3),
        )
    )
    inverse_operator_contributions = {
        name: expand_standard(
            mass_operators.get(name, sp.ZeroMatrix(3, 3))
            - momentum_squared * kinetic_operators.get(name, sp.ZeroMatrix(3, 3))
        )
        for name in ("LO", *couplings)
    }
    leading_operator = expand_standard(mass_operators["LO"] - momentum_squared * kinetic_operators["LO"])
    correction_operator = sum((couplings[name] * expression for name, expression in mass_operators.items() if name != "LO"), sp.ZeroMatrix(3, 3))
    correction_operator -= momentum_squared * sum((couplings[name] * expression for name, expression in kinetic_operators.items() if name != "LO"), sp.ZeroMatrix(3, 3))
    total_operator = expand_standard(leading_operator + correction_operator)
    # Only three traces involving Phi occur. Scalar markers isolate their
    # coefficient matrices using standard scalar differentiation in SymPy.
    trace_moments = sp.symbols("T_0:3", real=True)
    trace_substitution = {}
    for power, marker in enumerate(trace_moments):
        trace_substitution[sp.Trace(flavor_matrix * chi**power)] = marker
        trace_substitution[sp.Trace(chi**power * flavor_matrix)] = marker
    trace_substitution[sp.Trace(flavor_matrix)] = trace_moments[0]
    marked_operator = total_operator.xreplace(trace_substitution)
    diagonal_operator = marked_operator.subs(dict.fromkeys(trace_moments, 0))
    trace_operator = expand_standard(total_operator - diagonal_operator)
    trace_columns = [sp.diff(marked_operator, marker) for marker in trace_moments]
    polynomial_variable = sp.Symbol("z", real=True)
    trace_coefficients = sp.Matrix(3, 3, lambda row, column: sp.expand(trace_columns[column].subs(chi, polynomial_variable * sp.Identity(3)).doit().as_explicit()[0, 0]).coeff(polynomial_variable, row))
    eigenvalue_left, eigenvalue_right = sp.symbols("chi_i chi_j", real=True)
    trace_chi = sp.Symbol("tr_chi_0", real=True)
    # D_a is a polynomial in left/right multiplication by chi. Its component
    # denominator follows by applying the derived operator to a matrix unit.
    matrix_unit = sp.zeros(3)
    matrix_unit[0, 1] = 1
    denominator_expression = diagonal_operator.subs({chi: sp.diag(eigenvalue_left, eigenvalue_right, trace_chi - eigenvalue_left - eigenvalue_right), flavor_matrix: matrix_unit}).doit().as_explicit()[0, 1]
    denominator_expression = sp.factor(denominator_expression)
    matrix_unit_action = diagonal_operator.subs({
        chi: sp.diag(eigenvalue_left, eigenvalue_right, trace_chi - eigenvalue_left - eigenvalue_right),
        flavor_matrix: matrix_unit,
    }).doit().as_explicit().applyfunc(sp.factor)
    # Derive the LO trace equation and its solution, rather than insert it.
    trace_inverse_source, trace_inverse_identity, trace_solution = sp.symbols("tr_Dinv_C tr_Dinv_1 tr_S")
    anomaly_coefficient = sources["m_zero_squared"] / 3
    trace_equation = sp.Eq(trace_solution, trace_inverse_source - anomaly_coefficient * trace_solution * trace_inverse_identity)
    leading_trace_solution = sp.solve(trace_equation, trace_solution)[0]
    inverse_leading_source = sp.MatrixSymbol("D_{a,0}^{-1}[Y]", 3, 3)
    inverse_leading_identity = sp.MatrixSymbol("D_{a,0}^{-1}[I]", 3, 3)
    leading_inverse = inverse_leading_source - anomaly_coefficient * leading_trace_solution.subs({
        trace_inverse_source: sp.Trace(inverse_leading_source),
        trace_inverse_identity: sp.Trace(inverse_leading_identity),
    }) * inverse_leading_identity
    moment_matrix = sp.MatrixSymbol("G", 3, 3)
    coefficient_matrix = sp.MatrixSymbol("H", 3, 3)
    source_moments = sp.MatrixSymbol("b", 3, 1)
    solution_moments = sp.MatrixSymbol("T", 3, 1)
    trace_system_equation = sp.Eq(
        solution_moments,
        source_moments - moment_matrix * coefficient_matrix * solution_moments,
    )
    collected_trace_equation = sp.Eq(
        (sp.Identity(3) + moment_matrix * coefficient_matrix) * solution_moments,
        source_moments,
    )
    solved_moments = (sp.Identity(3) + moment_matrix * coefficient_matrix).I * source_moments
    trace_weights = coefficient_matrix * solved_moments
    source_matrix = sp.MatrixSymbol("Y", 3, 3)
    identity_weight, linear_weight, quadratic_weight = sp.symbols("q_0 q_1 q_2")
    inverse_source_argument = source_matrix - identity_weight * sp.Identity(3) - linear_weight * chi - quadratic_weight * chi**2
    # The undefined function is only a display of the linear operator applied
    # to its derived matrix argument. The executable inverse is below.
    exact_truncated_inverse = sp.Function(r"\mathbb{D}_a^{-1}")(inverse_source_argument)
    # The first-order solution is obtained by equating coefficients in the
    # operator equation, retaining the order of all matrix products.
    expansion_parameter = sp.Symbol("epsilon")
    leading_matrix = sp.MatrixSymbol("M_a", 9, 9)
    correction_matrix = sp.MatrixSymbol("Delta_M_a", 9, 9)
    leading_source = sp.MatrixSymbol("C", 9, 1)
    correction_source = sp.MatrixSymbol("Delta_C", 9, 1)
    leading_solution = leading_matrix.I * leading_source
    source_remainder = correction_source - correction_matrix * leading_solution
    correction_solution = leading_matrix.I * source_remainder
    expanded_solution = leading_solution + expansion_parameter * correction_solution
    expansion_residual = sp.expand((leading_matrix + expansion_parameter * correction_matrix) * expanded_solution - leading_source - expansion_parameter * correction_source)
    return locals()


def inverse_from_trace_system(source, chi_values, momentum_squared, anomaly_coefficient, coefficient_values=None):
    """Solve the derived diagonal operator plus its three trace terms."""
    result = derive()
    substitutions = dict.fromkeys(result["couplings"].values(), 0)
    substitutions.update(coefficient_values or {})
    substitutions.update({result["f_pi"]: 1, result["sources"]["m_zero_squared"]: 3 * anomaly_coefficient, result["momentum_squared"]: momentum_squared, result["trace_chi"]: sum(chi_values)})
    denominator = sp.Matrix(3, 3, lambda row, column: result["denominator_expression"].subs(substitutions).subs({result["eigenvalue_left"]: chi_values[row], result["eigenvalue_right"]: chi_values[column]}))
    trace_coefficients = result["trace_coefficients"].subs(substitutions)
    inverse_source = sp.Matrix(3, 3, lambda row, column: source[row, column] / denominator[row, column])
    trace_matrix = sp.Matrix(3, 3, lambda row, column: sum(chi_values[index]**(row + column) / denominator[index, index] for index in range(3)))
    source_traces = sp.Matrix([sum(chi_values[index]**power * inverse_source[index, index] for index in range(3)) for power in range(3)])
    solved_traces = (sp.eye(3) + trace_matrix * trace_coefficients).inv() * source_traces
    correction_weights = trace_coefficients * solved_traces
    solution = inverse_source - sp.diag(*[sum(chi_values[index]**power * correction_weights[power] for power in range(3)) / denominator[index, index] for index in range(3)])
    return solution.applyfunc(sp.simplify), result, substitutions


class OperatorInverseTests(unittest.TestCase):
    def test_operator_normalization_recovers_quadratic_forms(self):
        result = derive()
        for name, expression in result["sources"]["meson_kinetic"].items():
            variable = result["sources"]["derivative_phi"]
            reconstructed = sp.Trace(variable * result["kinetic_operators"][name].subs(result["flavor_matrix"], variable)) / 4
            self.assertEqual(expand_standard(reconstructed - expression), 0)
        for name, expression in result["sources"]["meson_potential"].items():
            reconstructed = -sp.Trace(result["flavor_matrix"] * result["mass_operators"][name]) / 4
            self.assertEqual(expand_standard(reconstructed - expression), 0)

    def test_leading_diagonal_reduction(self):
        source = sp.diag(2, 4, 7)
        chi_values = [sp.Integer(3), sp.Integer(5), sp.Integer(8)]
        solution, _, _ = inverse_from_trace_system(source, chi_values, 1, sp.Rational(2, 3))
        trace_value = sum(source[index, index] / (chi_values[index] - 1) for index in range(3)) / (1 + sp.Rational(2, 3) * sum(1 / (value - 1) for value in chi_values))
        expected = sp.diag(*[(source[index, index] - sp.Rational(2, 3) * trace_value) / (chi_values[index] - 1) for index in range(3)])
        self.assertEqual((solution - expected).applyfunc(sp.simplify), sp.zeros(3))

    def test_general_hermitian_source_solves_full_operator(self):
        source = sp.Matrix([[1, 2 + sp.I, 3], [2 - sp.I, 4, 1 - sp.I], [3, 1 + sp.I, 2]])
        initial = derive()
        values = {symbol: sp.Rational(index + 1, 1000) for index, symbol in enumerate(initial["couplings"].values())}
        solution, result, substitutions = inverse_from_trace_system(source, [2, 3, 5], sp.Rational(1, 2), sp.Rational(2, 3), values)
        applied = result["total_operator"].subs(substitutions).subs({result["chi"]: sp.diag(2, 3, 5), result["flavor_matrix"]: solution}).doit().as_explicit()
        self.assertEqual((applied - source).applyfunc(sp.simplify), sp.zeros(3))
        self.assertEqual((solution - solution.H).applyfunc(sp.simplify), sp.zeros(3))

    def test_first_order_operator_equation(self):
        result = derive()
        for power in (0, 1):
            coefficient = sp.expand(result["expansion_residual"]).coeff(result["expansion_parameter"], power)
            self.assertTrue(coefficient == 0 or coefficient == sp.ZeroMatrix(9, 1))


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    result = derive()
    print("Natural units; general U(3) matrices; exact inverse of retained local terms and strict first-order expansion; no loops.")
    print("1. Derive the kinetic and mass operators from the quadratic Lagrangian.")
    print("The U(3) normalization is L_kin = tr(Z K[Z])/4 and L_pot = -tr(Phi M^2[Phi])/4.")
    print("Thus K[Z] is twice the kinetic trace gradient, and M^2[Phi] is minus twice the potential trace gradient.")
    for name, expression in result["sources"]["meson_kinetic"].items():
        show(f"Input quadratic kinetic term, coefficient of {name}", expression)
        show(f"Twice its trace gradient, K_{name}[Phi]", result["kinetic_operators"][name])
    show("Full kinetic operator K[Phi]", result["total_kinetic_operator"])
    for name, expression in result["sources"]["meson_potential"].items():
        show(f"Input quadratic potential term, coefficient of {name}", expression)
        show(f"Minus twice its trace gradient, M^2_{name}[Phi]", result["mass_operators"][name])
    show("Full mass operator M^2[Phi]", result["total_mass_operator"])
    print("2. At p^2=w the operator to invert is M^2-wK.")
    for name, expression in result["inverse_operator_contributions"].items():
        show(f"Coefficient of {name} in (M^2-wK)[Phi]", expression)
    show("Full local NLO operator M_a[Phi], in matrices and traces", result["total_operator"])
    print("3. Separate terms with Phi outside all traces from terms depending on its three trace moments.")
    show("Terms with Phi outside all traces, D_a[Phi]", result["diagonal_operator"])
    show("Remaining trace operator (M_a-D_a)[Phi]", result["trace_operator"])
    show("Trace coefficient matrix H", result["trace_coefficients"])
    print("4. Component-level derivation: apply D_a to the matrix unit E_12 in the quark-mass basis.")
    show("Component-level representative matrix unit E_12", result["matrix_unit"])
    show("Component-level D_a[E_12]", result["matrix_unit_action"])
    show("Component-level denominator d_ij of D_a inverse", result["denominator_expression"])
    print("Consequently D_a is inverted by dividing each flavor component by d_ij; the trace operator is the only remaining coupling between components.")
    print("5. Solve the trace moments after applying D_a^{-1} to the matrix equation.")
    show("LO trace equation", result["trace_equation"])
    show("LO solved trace", result["leading_trace_solution"])
    show("LO rank-one inverse", result["leading_inverse"])
    show("Three-trace equation before collecting unknown moments", result["trace_system_equation"])
    show("Collected three-trace equation", result["collected_trace_equation"])
    show("Three trace moments solved", result["solved_moments"])
    show("Trace-system coefficient vector q", result["trace_weights"])
    show("Exact inverse of the truncated local operator", result["exact_truncated_inverse"])
    print("6. Expand the operator equation once for a strict first-order result.")
    show("Leading solution", result["leading_solution"])
    show("First correction source", result["source_remainder"])
    show("Strict first-order solution", result["expanded_solution"])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
