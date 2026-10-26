"""Expand the gauged WZW linear source and a candidate odd V3pi operator.

Run with .venv/bin/python scripts/derive_wzw_gamma_three_pion_and_vector_contact.py.
Natural units; Phi uses the Alda:2026nrd normalization U=exp(i Phi/f_pi).
The overall WZW orientation is denoted by ``orientation``; its magnitude
is fixed, but a signed matching to the archived functional is not attempted.
Only the source-linear cubic-pion part of the full functional is needed:
the ungauged five-form starts at five pions; other gauged terms contain
at least two external sources. The source-linear piece is represented as
orientation*N_c/(48*pi**2) epsilon tr[v(R R R + L L L)],
R=dU U^dagger, L=U^dagger dU, v=-e Q A.

Inputs: memory/extracts/Alda:2026nrd/appendix-a-expansion.md;
memory/rcht/foundations/rcht-basics.md;
memory/rcht/odd-parity/levi-civita-convention.md.
The VPPP operator definitions are GomezDumm:2010ra, Eq. (13).
Explicit charged-channel matrices are necessary to extract the charged
flavor coefficient. The covariant reduction uses integration by parts,
zero external sources, nabla_[alpha u_beta]=0, and the SU(N_f) leading
Goldstone EOM nabla.u=i(chi_minus-tr(chi_minus)/N_f)/2. Its double-trace
remainder is retained; singlet/topological U(3) EOM terms are not assumed.
"""

import itertools
import unittest

import sympy as sp
from sympy.combinatorics import Permutation


def show(label, expression):
    """Print deterministic symbolic derivation output in LaTeX."""
    print(f"{label}:\n{sp.latex(expression)}\n")


number_colors, pion_decay_constant, electric_charge = sp.symbols(
    "N_c f_pi e", positive=True, real=True
)
orientation = sp.Symbol("orientation", real=True)
non_strange_entry, strange_entry = sp.symbols("t_n t_s", real=True)
vector_mass = sp.Symbol("M_i", positive=True)
pion_mass_squared, total_invariant = sp.symbols("m_pi_squared s", real=True)
pair_invariant_12, pair_invariant_13, pair_invariant_23 = sp.symbols(
    "s_12 s_13 s_23", real=True
)

# Phi = lambda_a pi_a; the charged fields have the conventional sqrt(2).
charged_plus = sp.sqrt(2) * sp.Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
charged_minus = charged_plus.T
neutral_pion = sp.diag(1, -1, 0)
pion_generators = (charged_plus, charged_minus, neutral_pion)
charge_generator = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
isoscalar_generator = sp.diag(non_strange_entry, non_strange_entry, strange_entry)


def antisymmetric_trace(source):
    """Extract the coefficient of d pi+ wedge d pi- wedge d pi0.

    SymPy does not expand matrix-valued differential forms directly;
    the six permutations implement this finite exterior-product coefficient.
    """
    result = sp.S.Zero
    for assignment in itertools.permutations(range(3)):
        result += Permutation(assignment).signature() * sp.trace(
            source * pion_generators[assignment[0]]
            * pion_generators[assignment[1]] * pion_generators[assignment[2]]
        )
    return sp.simplify(result)


flavor_wedge = antisymmetric_trace(charge_generator)
# Each Maurer-Cartan current starts at i*dPhi/f_pi. Their two placements
# contribute equally at cubic pion order; derive their coefficient directly.
field_order = sp.Symbol("field_order", real=True)
pion_matrix, derivative_pion_matrix = sp.symbols("Phi dPhi", commutative=False)
unitary_field_linear = 1 + sp.I * field_order * pion_matrix / pion_decay_constant
unitary_field_adjoint_linear = 1 - sp.I * field_order * pion_matrix / pion_decay_constant
derivative_unitary_field_linear = sp.I * field_order * derivative_pion_matrix / pion_decay_constant
right_current_linear = sp.expand(
    derivative_unitary_field_linear * unitary_field_adjoint_linear
).coeff(field_order, 1)
left_current_linear = sp.expand(
    unitary_field_adjoint_linear * derivative_unitary_field_linear
).coeff(field_order, 1)
current_product_coefficient = (
    right_current_linear.coeff(derivative_pion_matrix) ** 3
    + left_current_linear.coeff(derivative_pion_matrix) ** 3
)
wzw_charged_coefficient = sp.simplify(
    orientation * number_colors / (48 * sp.pi**2)
    * (-electric_charge) * current_product_coefficient * flavor_wedge
)

# Test i epsilon tr[nabla_alpha V^{mu nu}{u^alpha,[u^rho,u^sigma]}].
# C: V -> -V^T, u -> u^T. The three-u combination reverses sign on
# transposition, so the two signs compensate. Hermitian conjugation also
# changes its sign; the explicit i makes the operator Hermitian.
first_u, second_u, third_u, vector_field = sp.symbols("A B C V", commutative=False)
three_u_expression = first_u * (second_u * third_u - third_u * second_u) + (
    second_u * third_u - third_u * second_u
) * first_u


def reverse_words(expression):
    """Use SymPy adjoint to reverse products of Hermitian abstract blocks."""
    replacements = {sp.adjoint(item): item for item in (first_u, second_u, third_u, vector_field)}
    return sp.expand(sp.adjoint(expression).doit().xreplace(replacements))


charge_conjugation_remainder = sp.expand(reverse_words(three_u_expression) + three_u_expression)
hermiticity_remainder = sp.expand(sp.conjugate(sp.I) * reverse_words(three_u_expression) - sp.I * three_u_expression)

q_dot_pion = (
    (pair_invariant_12 + pair_invariant_13) / 2 - pion_mass_squared,
    (pair_invariant_12 + pair_invariant_23) / 2 - pion_mass_squared,
    (pair_invariant_13 + pair_invariant_23) / 2 - pion_mass_squared,
)

def contact_tensor_coefficient(source):
    """Coefficient of epsilon(p+,p-,p0,epsilon_i), before Fourier factors.

    Project V^{mu nu}=(q^mu epsilon^nu-q^nu epsilon^mu)/M_i. The
    epsilon contraction selects the pion not in the rho,sigma slots.
    Include the three leading u=-dPhi/f_pi factors and the explicit i.
    The connection starts at two pions and cannot enter this cubic vertex.
    """
    result = sp.S.Zero
    for assignment in itertools.permutations(range(3)):
        first, second, third = (pion_generators[index] for index in assignment)
        commutator = second * third - third * second
        flavor = sp.trace(source * (first * commutator + commutator * first))
        # epsilon(q,polarization,p_second,p_third) in the displayed basis.
        sign = Permutation([assignment[0], 3, assignment[1], assignment[2]]).signature()
        result += 2 * sign * q_dot_pion[assignment[0]] * flavor
    result *= -sp.I / (vector_mass * pion_decay_constant**3)
    result = sp.expand(result).subs(
        pair_invariant_23,
        total_invariant + 3 * pion_mass_squared - pair_invariant_12 - pair_invariant_13,
    )
    return sp.simplify(result)


contact_coefficient = contact_tensor_coefficient(isoscalar_generator)

# Covariant reduction at zero external sources. Preserve flavor ordering;
# cyclically place V first in every trace and antisymmetrize alpha,beta.
u_alpha, u_beta, u_gamma, h_alpha_gamma, h_beta_gamma = sp.symbols(
    "u_alpha u_beta u_gamma h_alpha_gamma h_beta_gamma", commutative=False
)
divergence_u, chi_minus = sp.symbols("divergence_u chi_minus", commutative=False)
trace_chi_minus = sp.Symbol("trace_chi_minus", commutative=True)
flavor_number = sp.Symbol("N_f", positive=True, integer=True)
index_exchange = {
    u_alpha: u_beta, u_beta: u_alpha,
    h_alpha_gamma: h_beta_gamma, h_beta_gamma: h_alpha_gamma,
}


def epsilon_antisymmetrize(expression):
    """Project the alpha,beta slots with the antisymmetric epsilon tensor."""
    return sp.expand((expression - expression.xreplace(index_exchange)) / 2)


commutator_alpha_beta = u_alpha * u_beta - u_beta * u_alpha
# Integration by parts differentiates all three u factors. At sources=0,
# MC gives nabla_gamma u_alpha = h_alpha_gamma/2 (and likewise beta).
derivative_commutator = (
    h_alpha_gamma * u_beta - u_beta * h_alpha_gamma
    + u_alpha * h_beta_gamma - h_beta_gamma * u_alpha
) / 2
contact_after_ibp = -sp.I * vector_field * (
    divergence_u * commutator_alpha_beta + commutator_alpha_beta * divergence_u
    + u_gamma * derivative_commutator + derivative_commutator * u_gamma
)
contact_after_su3_eom = epsilon_antisymmetrize(contact_after_ibp.subs(
    divergence_u, sp.I * (chi_minus - trace_chi_minus / flavor_number) / 2
))

literature_operators = (
    sp.I * vector_field * (h_alpha_gamma * u_gamma * u_beta - u_beta * u_gamma * h_alpha_gamma),
    sp.I * vector_field * (h_alpha_gamma * u_beta * u_gamma - u_gamma * u_beta * h_alpha_gamma),
    sp.I * vector_field * (u_gamma * h_alpha_gamma * u_beta - u_beta * h_alpha_gamma * u_gamma),
    vector_field * (u_alpha * u_beta * chi_minus + chi_minus * u_alpha * u_beta),
    vector_field * u_beta * chi_minus * u_alpha,
)
double_trace_operator = trace_chi_minus * vector_field * u_alpha * u_beta
reduction_operators = tuple(epsilon_antisymmetrize(item) for item in literature_operators) + (
    epsilon_antisymmetrize(double_trace_operator),
)
reduction_coefficients = sp.symbols("c_1:7", real=True)
matching_expression = sp.expand(contact_after_su3_eom - sum(
    coefficient * operator for coefficient, operator in zip(reduction_coefficients, reduction_operators)
))
# SymPy does not offer Poly for noncommuting generators. args_cnc exposes
# its native noncommuting monomials; scalar matching remains ordinary solve.
noncommuting_words = sorted({
    sp.Mul(*term.args_cnc()[1]) for term in sp.Add.make_args(matching_expression)
}, key=str)
matching_equations = []
for word in noncommuting_words:
    matching_equations.append(sp.Add(*(
        sp.Mul(*term.args_cnc()[0])
        for term in sp.Add.make_args(matching_expression)
        if sp.Mul(*term.args_cnc()[1]) == word
    )))
reduction_solution = sp.solve(matching_equations, reduction_coefficients, dict=True)[0]
reduction_remainder = sp.expand(matching_expression.subs(reduction_solution))


class DerivationChecks(unittest.TestCase):
    """Check anomaly magnitude, symmetries and physical flavor limits."""

    def test_anomaly_magnitude(self):
        self.assertEqual(sp.simplify(wzw_charged_coefficient / (sp.I * orientation * electric_charge)),
                         number_colors / (12 * sp.pi**2 * pion_decay_constant**3))

    def test_operator_charge_conjugation_and_hermiticity(self):
        self.assertEqual(charge_conjugation_remainder, 0)
        self.assertEqual(hermiticity_remainder, 0)

    def test_nonzero_isoscalar_and_ideal_phi(self):
        self.assertNotEqual(contact_coefficient, 0)
        self.assertEqual(contact_coefficient.subs(non_strange_entry, 0), 0)
        self.assertEqual(sp.diff(contact_coefficient, strange_entry), 0)

    def test_isovector_contact_absent_in_exact_isospin(self):
        self.assertEqual(contact_tensor_coefficient(sp.diag(1, -1, 0) / sp.sqrt(2)), 0)

    def test_covariant_reduction_with_trace_projected_eom(self):
        self.assertEqual(reduction_remainder, 0)
        self.assertEqual(reduction_solution[reduction_coefficients[5]], -2 / flavor_number)

    def test_isoscalar_pion_double_trace_remainder_absent(self):
        for first, second in itertools.product(pion_generators, repeat=2):
            self.assertEqual(sp.trace(isoscalar_generator * (first * second - second * first)), 0)


if __name__ == "__main__":
    show("Linear unitary field", unitary_field_linear)
    show("Linear right Maurer-Cartan current", right_current_linear)
    show("Linear left Maurer-Cartan current", left_current_linear)
    show("Charged flavor exterior-product trace", flavor_wedge)
    show("Charged WZW coefficient of epsilon A d pi+ d pi- d pi0", wzw_charged_coefficient)
    show("Three-u flavor ordering", three_u_expression)
    show("Charge-conjugation remainder", charge_conjugation_remainder)
    show("Hermiticity remainder", hermiticity_remainder)
    show("Candidate contact tensor coefficient before Fourier factors", contact_coefficient)
    show("Candidate after IBP and zero-source Maurer-Cartan relation (trace understood)", contact_after_ibp)
    show("Candidate after trace-projected Goldstone EOM (epsilon and trace understood)", contact_after_su3_eom)
    show("Derived coefficients for O1,O2,O3,O4,O5 and epsilon trace(chi_minus) trace(V u_alpha u_beta)",
         sp.Tuple(*(reduction_solution[item] for item in reduction_coefficients)))
    show("Covariant reduction remainder", reduction_remainder)
    unittest.main(argv=[__file__])
