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
    # Derive the LO trace equation and its solution, rather than insert it.
    trace_inverse_source, trace_inverse_identity, trace_solution = sp.symbols("tr_Dinv_C tr_Dinv_1 tr_S")
    anomaly_coefficient = sources["m_zero_squared"] / 3
    trace_equation = sp.Eq(trace_solution, trace_inverse_source - anomaly_coefficient * trace_solution * trace_inverse_identity)
    leading_trace_solution = sp.solve(trace_equation, trace_solution)[0]
    inverse_diagonal_source = sp.MatrixSymbol("D_a^{-1}(C)", 3, 3)
    inverse_diagonal_powers = [sp.MatrixSymbol(f"D_a^{{-1}}(chi_0^{power})", 3, 3) for power in range(3)]
    leading_inverse = inverse_diagonal_source - anomaly_coefficient * leading_trace_solution * inverse_diagonal_powers[0]
    moment_matrix = sp.MatrixSymbol("G", 3, 3)
    coefficient_matrix = sp.MatrixSymbol("H", 3, 3)
    source_moments = sp.MatrixSymbol("b", 3, 1)
    solved_moments = (sp.Identity(3) + moment_matrix * coefficient_matrix).I * source_moments
    moment_solution_symbol = sp.MatrixSymbol("T", 3, 1)
    trace_weights = coefficient_matrix * moment_solution_symbol
    exact_truncated_inverse = inverse_diagonal_source - sum((inverse_diagonal_powers[power] * trace_weights[power, 0] for power in range(3)), sp.ZeroMatrix(3, 3))
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
    show("LO operator", result["leading_operator"])
    for name, expression in result["mass_operators"].items():
        if name != "LO":
            show(f"Mass correction operator, {name}", expression)
    for name, expression in result["kinetic_operators"].items():
        if name != "LO":
            show(f"Kinetic correction operator, {name}", expression)
    show("Diagonal anticommutator part D_a[Phi]", result["diagonal_operator"])
    show("Trace coefficient matrix H", result["trace_coefficients"])
    show("Denominator d_ij of D_a inverse", result["denominator_expression"])
    show("LO trace equation", result["trace_equation"])
    show("LO solved trace", result["leading_trace_solution"])
    show("LO rank-one inverse", result["leading_inverse"])
    show("Three trace moments solved", result["solved_moments"])
    show("Exact inverse of the truncated local operator", result["exact_truncated_inverse"])
    show("Leading solution", result["leading_solution"])
    show("First correction source", result["source_remainder"])
    show("Strict first-order solution", result["expanded_solution"])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
