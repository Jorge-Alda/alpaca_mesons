"""First-order dressed quark mass source in the Alda:2026nrd parametrization.

Run (from the repository root):
    .venv/bin/python scripts/derive_alp_dressed_mass_first_order.py

Conventions and units
---------------------
Natural units; x = a/f_a is dimensionless. All flavor matrices are
noncommuting. The quark mass function M(x) = S(x) + i P(x) multiplies
bar(q_R) M q_L and has mass dimension one; C_R = C_V + C_A and
C_L = C_V - C_A are the dimensionless Hermitian coefficients of the
right/left derivative currents (d_mu a / f_a) bar(q) gamma^mu C q.
The transport matrices solve dW_eta/dx = i C_eta W_eta with W_eta(0) = 1,
and the dressed mass source is hat(M) = W_R^dagger M W_L. Taylor
coefficients are derivatives at x = 0.

Inputs
------
* memory/alp/shared/alp-general-chiral-redefinition-invariants.md:
  definitions of W_R, W_L, hat(M), M = S + i P, C_R, C_L, and the
  first-order transformation laws of S_1, P_1, C_V0, C_A0 under the constant
  redefinition g_R = exp(i x (delta + kappa)), g_L = exp(i x (delta - kappa)).
* memory/extracts/Alda:2026nrd/appendix-a-expansion.md: source
  parametrization S(x) = m + i x [c^S, m], P(x) = x {c^P, m},
  C_V(0) = c^V, C_A(0) = c^A, with m, c^V, c^A, c^S, c^P Hermitian, the
  invariants c^VS = c^V - c^S, c^AP = c^A - c^P, and the representative
  chosen by delta = -c^S, kappa = -c^A - M_a^{-1}(C).

Assumption: only the O(x) Taylor coefficient is studied; the transport
matrices are truncated at O(x), which is exact for this coefficient.

Each check prints PASS or FAIL; the process exits with status 1 on failure.
"""

import random
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_alp_invariant_coefficient_functions import transport_series  # noqa: E402

ALP_FIELD_RATIO = sp.Symbol("x", real=True)


def hermitian_matrix_symbol(name):
    return sp.Symbol(name, commutative=False, hermitian=True)


QUARK_MASS = hermitian_matrix_symbol("m_q")
VECTOR_COUPLING = hermitian_matrix_symbol("c^V")
AXIAL_COUPLING = hermitian_matrix_symbol("c^A")
SCALAR_COUPLING = hermitian_matrix_symbol("c^S")
PSEUDOSCALAR_COUPLING = hermitian_matrix_symbol("c^P")
VECTOR_SCALAR_INVARIANT = hermitian_matrix_symbol("c^{VS}")
AXIAL_PSEUDOSCALAR_INVARIANT = hermitian_matrix_symbol("c^{AP}")
VECTOR_GENERATOR = hermitian_matrix_symbol(r"\delta")
AXIAL_GENERATOR = hermitian_matrix_symbol(r"\kappa")


def commutator(first, second):
    return first * second - second * first


def anticommutator(first, second):
    return first * second + second * first


def dressed_mass_first_order(right_current_zero, left_current_zero, mass_zero, mass_one):
    """Return (hat M_0, hat M_1) from W_R^dagger M W_L truncated at O(x)."""
    right_transport = transport_series([right_current_zero], ALP_FIELD_RATIO, 1)
    left_transport = transport_series([left_current_zero], ALP_FIELD_RATIO, 1)
    mass_function = mass_zero + ALP_FIELD_RATIO * mass_one
    dressed = sp.expand(sp.adjoint(right_transport) * mass_function * left_transport)
    return (sp.expand(dressed.subs(ALP_FIELD_RATIO, 0)),
            sp.expand(sp.diff(dressed, ALP_FIELD_RATIO).subs(ALP_FIELD_RATIO, 0)))


def alda_sources(vector, axial, scalar, pseudoscalar):
    """Return (C_R0, C_L0, M_0, M_1) in the Alda:2026nrd parametrization."""
    scalar_one = sp.I * commutator(scalar, QUARK_MASS)
    pseudoscalar_one = anticommutator(pseudoscalar, QUARK_MASS)
    return (vector + axial, vector - axial, QUARK_MASS,
            scalar_one + sp.I * pseudoscalar_one)


def closed_form_dressed_mass_one(vector_scalar, axial_pseudoscalar):
    return sp.expand(-sp.I * commutator(vector_scalar, QUARK_MASS)
                     - sp.I * anticommutator(axial_pseudoscalar, QUARK_MASS))


def alda_dressed_mass_one():
    return dressed_mass_first_order(*alda_sources(
        VECTOR_COUPLING, AXIAL_COUPLING, SCALAR_COUPLING, PSEUDOSCALAR_COUPLING))


def check_general_and_alda_form():
    """Check 1: general first-order dressing and its Alda specialization."""
    right_zero = hermitian_matrix_symbol("C_R0")
    left_zero = hermitian_matrix_symbol("C_L0")
    mass_zero, mass_one = sp.symbols("M_0 M_1", commutative=False)
    dressed_zero, dressed_one = dressed_mass_first_order(right_zero, left_zero,
                                                         mass_zero, mass_one)
    general_residual = sp.expand(
        dressed_one - (mass_one - sp.I * (right_zero * mass_zero - mass_zero * left_zero)))
    alda_zero, alda_one = alda_dressed_mass_one()
    alda_residual = sp.expand(alda_one - closed_form_dressed_mass_one(
        VECTOR_COUPLING - SCALAR_COUPLING, AXIAL_COUPLING - PSEUDOSCALAR_COUPLING))
    return {"general_zero_residual": sp.expand(dressed_zero - mass_zero),
            "general_one_residual": general_residual,
            "alda_zero_residual": sp.expand(alda_zero - QUARK_MASS),
            "alda_one_residual": alda_residual,
            "alda_one": alda_one}


def hermitian_parts(matrix_expression):
    adjoint_expression = sp.expand(sp.adjoint(matrix_expression))
    return (sp.expand((matrix_expression + adjoint_expression) / 2),
            sp.expand((matrix_expression - adjoint_expression) / (2 * sp.I)))


def check_hermitian_split():
    """Check 2: hat S_1 = -i[c^VS, m] and hat P_1 = -{c^AP, m}."""
    dressed_one = closed_form_dressed_mass_one(VECTOR_SCALAR_INVARIANT,
                                               AXIAL_PSEUDOSCALAR_INVARIANT)
    scalar_one, pseudoscalar_one = hermitian_parts(dressed_one)
    return {"scalar_one": scalar_one, "pseudoscalar_one": pseudoscalar_one,
            "scalar_residual": sp.expand(
                scalar_one + sp.I * commutator(VECTOR_SCALAR_INVARIANT, QUARK_MASS)),
            "pseudoscalar_residual": sp.expand(
                pseudoscalar_one + anticommutator(AXIAL_PSEUDOSCALAR_INVARIANT, QUARK_MASS))}


def check_constant_redefinition_invariance():
    """Check 3: hat M_1 is unchanged by the first-order constant redefinition."""
    right_zero, left_zero, mass_zero, mass_one = alda_sources(
        VECTOR_COUPLING, AXIAL_COUPLING, SCALAR_COUPLING, PSEUDOSCALAR_COUPLING)
    _, dressed_one = dressed_mass_first_order(right_zero, left_zero, mass_zero, mass_one)
    scalar_one = sp.I * commutator(SCALAR_COUPLING, QUARK_MASS)
    pseudoscalar_one = anticommutator(PSEUDOSCALAR_COUPLING, QUARK_MASS)
    # Laws from the invariants note with S_0 = m, P_0 = 0.
    scalar_one_new = scalar_one + sp.I * commutator(VECTOR_GENERATOR, QUARK_MASS)
    pseudoscalar_one_new = pseudoscalar_one + anticommutator(AXIAL_GENERATOR, QUARK_MASS)
    vector_zero_new = VECTOR_COUPLING + VECTOR_GENERATOR
    axial_zero_new = AXIAL_COUPLING + AXIAL_GENERATOR
    _, dressed_one_from_laws = dressed_mass_first_order(
        vector_zero_new + axial_zero_new, vector_zero_new - axial_zero_new,
        QUARK_MASS, scalar_one_new + sp.I * pseudoscalar_one_new)
    # Same transformation written as source shifts c^{V,S} += delta, c^{A,P} += kappa.
    _, dressed_one_from_shifts = dressed_mass_first_order(*alda_sources(
        VECTOR_COUPLING + VECTOR_GENERATOR, AXIAL_COUPLING + AXIAL_GENERATOR,
        SCALAR_COUPLING + VECTOR_GENERATOR, PSEUDOSCALAR_COUPLING + AXIAL_GENERATOR))
    # Independent: M' = g_R M g_L^dagger with g = 1 + i x G, at O(x).
    right_group = 1 + sp.I * ALP_FIELD_RATIO * (VECTOR_GENERATOR + AXIAL_GENERATOR)
    left_group = 1 + sp.I * ALP_FIELD_RATIO * (VECTOR_GENERATOR - AXIAL_GENERATOR)
    transformed_mass = sp.expand(right_group * (mass_zero + ALP_FIELD_RATIO * mass_one)
                                 * sp.adjoint(left_group))
    transformed_mass_one = transformed_mass.coeff(ALP_FIELD_RATIO, 1)
    law_residual = sp.expand(transformed_mass_one
                             - (scalar_one_new + sp.I * pseudoscalar_one_new))
    return {"law_residual": law_residual,
            "from_laws_residual": sp.expand(dressed_one_from_laws - dressed_one),
            "from_shifts_residual": sp.expand(dressed_one_from_shifts - dressed_one)}


def random_rational_hermitian(generator, size):
    entries = sp.zeros(size, size)
    for row in range(size):
        entries[row, row] = sp.Rational(generator.randint(-9, 9), generator.randint(1, 5))
        for column in range(row + 1, size):
            value = (sp.Rational(generator.randint(-9, 9), generator.randint(1, 5))
                     + sp.I * sp.Rational(generator.randint(-9, 9), generator.randint(1, 5)))
            entries[row, column] = value
            entries[column, row] = sp.conjugate(value)
    return entries


def numeric_couplings(seed=20261111, size=3):
    generator = random.Random(seed)
    auxiliary = random_rational_hermitian(generator, size)
    quark_mass = (auxiliary.H * auxiliary + sp.eye(size)).applyfunc(sp.expand)  # Hermitian positive definite
    couplings = [random_rational_hermitian(generator, size) for _ in range(4)]
    return quark_mass, couplings


def check_numeric_rational(seed=20261111, size=3):
    """Check 4: exact rational matrices, direct dressing versus closed form."""
    quark_mass, (vector, axial, scalar, pseudoscalar) = numeric_couplings(seed, size)
    identity = sp.eye(size)
    right_transport = identity + sp.I * ALP_FIELD_RATIO * (vector + axial)
    left_transport = identity + sp.I * ALP_FIELD_RATIO * (vector - axial)
    mass_function = quark_mass + ALP_FIELD_RATIO * (
        sp.I * (scalar * quark_mass - quark_mass * scalar)
        + sp.I * (pseudoscalar * quark_mass + quark_mass * pseudoscalar))
    dressed = (right_transport.H * mass_function * left_transport).applyfunc(sp.expand)
    dressed_one = dressed.applyfunc(lambda entry: entry.coeff(ALP_FIELD_RATIO, 1))
    vector_scalar = vector - scalar
    axial_pseudoscalar = axial - pseudoscalar
    closed_form = (-sp.I * (vector_scalar * quark_mass - quark_mass * vector_scalar)
                   - sp.I * (axial_pseudoscalar * quark_mass + quark_mass * axial_pseudoscalar))
    residual = (dressed_one - closed_form).applyfunc(sp.expand)
    return {"residual_is_zero": residual == sp.zeros(size, size),
            "quark_mass": quark_mass}


def linear_map_matrix(quark_mass, use_anticommutator):
    """Matrix of c -> [c, m] or {c, m} acting on row-major vec(c), complex N^2."""
    size = quark_mass.shape[0]
    columns = []
    for index in range(size * size):
        basis = sp.zeros(size, size)
        basis[index // size, index % size] = 1
        image = (basis * quark_mass + quark_mass * basis if use_anticommutator
                 else basis * quark_mass - quark_mass * basis)
        columns.append(sp.Matrix(size * size, 1, list(image)))
    return sp.Matrix.hstack(*columns)


def kernel_dimension(quark_mass, use_anticommutator):
    """Complex kernel dimension. Both maps commute with the adjoint up to a sign
    (after the factor -i for the commutator), so it equals the real kernel
    dimension on Hermitian matrices."""
    return len(linear_map_matrix(quark_mass, use_anticommutator).nullspace())


def check_kernel_structure():
    """Check 5: {., m} invertible for positive m; [., m] has a kernel."""
    up_mass, down_mass, strange_mass = sp.symbols("m_u m_d m_s", positive=True)
    symbolic_diagonal = sp.diag(up_mass, down_mass, strange_mass)
    anticommutator_determinant = sp.factor(
        linear_map_matrix(symbolic_diagonal, True).det())
    diagonal_coupling = sp.diag(*sp.symbols("c_1 c_2 c_3"))
    diagonal_scalar_one = (-sp.I * (diagonal_coupling * symbolic_diagonal
                                    - symbolic_diagonal * diagonal_coupling))
    nondegenerate = sp.diag(1, 2, 7)
    isospin_degenerate = sp.diag(1, 1, 7)
    generic_mass, _ = numeric_couplings()
    return {
        "anticommutator_determinant": anticommutator_determinant,
        "diagonal_scalar_one_is_zero": diagonal_scalar_one == sp.zeros(3, 3),
        "commutator_kernel_nondegenerate": kernel_dimension(nondegenerate, False),
        "commutator_kernel_isospin_degenerate": kernel_dimension(isospin_degenerate, False),
        "commutator_kernel_generic_hermitian": kernel_dimension(generic_mass, False),
        "anticommutator_kernel_nondegenerate": kernel_dimension(nondegenerate, True),
        "anticommutator_kernel_isospin_degenerate": kernel_dimension(isospin_degenerate, True),
        "anticommutator_kernel_generic_hermitian": kernel_dimension(generic_mass, True),
    }


def check_representative_consistency():
    """Check 6: Alda's final representative gives the same hat M_1."""
    axial_hat = hermitian_matrix_symbol(r"\hat c^A")
    _, dressed_from_representative = dressed_mass_first_order(*alda_sources(
        VECTOR_SCALAR_INVARIANT, axial_hat,
        sp.Integer(0), axial_hat - AXIAL_PSEUDOSCALAR_INVARIANT))
    representative_residual = sp.expand(
        dressed_from_representative
        - closed_form_dressed_mass_one(VECTOR_SCALAR_INVARIANT, AXIAL_PSEUDOSCALAR_INVARIANT))
    # Explicit redefinition delta = -c^S, kappa = -c^A - M_a^{-1}(C).
    inverse_mixing_term = hermitian_matrix_symbol(r"\mathbb{M}_a^{-1}(C)")
    vector_shift = -SCALAR_COUPLING
    axial_shift = -AXIAL_COUPLING - inverse_mixing_term
    _, dressed_after_shift = dressed_mass_first_order(*alda_sources(
        VECTOR_COUPLING + vector_shift, AXIAL_COUPLING + axial_shift,
        SCALAR_COUPLING + vector_shift, PSEUDOSCALAR_COUPLING + axial_shift))
    _, dressed_original = alda_dressed_mass_one()
    return {"representative_residual": representative_residual,
            "explicit_shift_residual": sp.expand(dressed_after_shift - dressed_original),
            "independent_of_axial_hat": not dressed_from_representative.has(axial_hat)}


def report(label, passed):
    print(f"[{'PASS' if passed else 'FAIL'}] {label}")
    return passed


def main():
    results = []
    first = check_general_and_alda_form()
    print("hat M_1 (Alda parametrization):\n" + sp.latex(first["alda_one"]) + "\n")
    results.append(report("1a hat M_0 = M_0 and hat M_1 = M_1 - i(C_R0 M_0 - M_0 C_L0)",
                          first["general_zero_residual"] == 0
                          and first["general_one_residual"] == 0))
    results.append(report("1b Alda: hat M_0 = m, hat M_1 = -i[c^VS,m] - i{c^AP,m}",
                          first["alda_zero_residual"] == 0
                          and first["alda_one_residual"] == 0))

    split = check_hermitian_split()
    print("\nhat S_1:\n" + sp.latex(split["scalar_one"]))
    print("hat P_1:\n" + sp.latex(split["pseudoscalar_one"]) + "\n")
    results.append(report("2 hat S_1 = -i[c^VS,m], hat P_1 = -{c^AP,m}",
                          split["scalar_residual"] == 0
                          and split["pseudoscalar_residual"] == 0))

    invariance = check_constant_redefinition_invariance()
    results.append(report("3a first-order laws reproduce g_R M g_L^dagger at O(x)",
                          invariance["law_residual"] == 0))
    results.append(report("3b hat M_1 invariant under first-order laws",
                          invariance["from_laws_residual"] == 0))
    results.append(report("3c hat M_1 invariant under c^{V,S}+=delta, c^{A,P}+=kappa",
                          invariance["from_shifts_residual"] == 0))

    numeric = check_numeric_rational()
    print("\nnumeric m (Hermitian positive definite, seed 20261111):\n"
          + sp.latex(numeric["quark_mass"]))
    results.append(report("4 exact rational 3x3 dressing equals closed form",
                          numeric["residual_is_zero"]))

    kernel = check_kernel_structure()
    print("\ndet of {., diag(m_u,m_d,m_s)} on 3x3 matrices:\n"
          + sp.latex(kernel["anticommutator_determinant"]))
    for key in ("commutator_kernel_nondegenerate", "commutator_kernel_isospin_degenerate",
                "commutator_kernel_generic_hermitian", "anticommutator_kernel_nondegenerate",
                "anticommutator_kernel_isospin_degenerate",
                "anticommutator_kernel_generic_hermitian"):
        print(f"{key}: {kernel[key]}")
    results.append(report("5a {., m} invertible for positive m (kernel 0)",
                          kernel["anticommutator_kernel_nondegenerate"] == 0
                          and kernel["anticommutator_kernel_isospin_degenerate"] == 0
                          and kernel["anticommutator_kernel_generic_hermitian"] == 0))
    results.append(report("5b diagonal c^VS with diagonal m gives hat S_1 = 0",
                          kernel["diagonal_scalar_one_is_zero"]))
    results.append(report("5c [., m] kernel dimension: 3 (nondegenerate), 5 (m,m,m_s), 3 (generic)",
                          kernel["commutator_kernel_nondegenerate"] == 3
                          and kernel["commutator_kernel_isospin_degenerate"] == 5
                          and kernel["commutator_kernel_generic_hermitian"] == 3))

    representative = check_representative_consistency()
    results.append(report("6a Alda representative (hat c^S=0, hat c^V=c^VS) gives same hat M_1",
                          representative["representative_residual"] == 0
                          and representative["independent_of_axial_hat"]))
    results.append(report("6b explicit delta=-c^S, kappa=-c^A-M_a^{-1}(C) leaves hat M_1 unchanged",
                          representative["explicit_shift_residual"] == 0))

    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
