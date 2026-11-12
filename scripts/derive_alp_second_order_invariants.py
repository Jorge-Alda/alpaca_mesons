"""Second-order coefficients of the dressed ALP sources hat(M), hat(Theta), hat(Theta_gamma).

Run (from the repository root):
    .venv/bin/python scripts/derive_alp_second_order_invariants.py

Conventions and units
---------------------
Natural units. x = a/f_a is dimensionless and real. A subscript n is the n-th
DERIVATIVE at x = 0, F_n = d^n F/dx^n |_0, so that F(x) = sum_n F_n x^n / n!.
This is the convention of derive_alp_invariant_coefficient_functions.py,
derive_alp_dressed_mass_first_order.py, derive_alp_dressed_topological_invariant.py
and derive_alp_dressed_photon_invariant.py.

The quark mass function M(x) = S(x) + i P(x) multiplies bar(q_R) M q_L and has
mass dimension one (S, P Hermitian). The current coefficients C_R(x), C_L(x) of
(d_mu a/f_a) bar(q_eta) gamma^mu C_eta q_eta are dimensionless Hermitian flavor
matrices, with C_R = C_V + C_A and C_L = C_V - C_A. Theta(x) multiplies G tilde(G)
in the Alda:2026nrd normalization (Theta_1 = c_G); Theta_gamma(x) multiplies
(alpha_em/4pi) F tilde(F) (Theta_gamma,1 = c_gamma). Q = diag(2/3,-1/3,-1/3),
N_c = 3. All flavor matrices are noncommuting.

Definitions (memory/alp/shared/alp-general-chiral-redefinition-invariants.md):
    dW_eta/dx = i C_eta W_eta, W_eta(0) = 1,
    hat(M) = W_R^dagger M W_L, hat(S) = (hat M + hat M^dagger)/2,
    hat(P) = (hat M - hat M^dagger)/(2i),
    hat(Theta) = Theta + int_0^x tr C_A        (alp-dressed-topological-invariant.md),
    hat(Theta_gamma) = Theta_gamma + 2 N_c int_0^x tr(Q^2 C_A)
                                               (memory/alp/photon/alp-dressed-photon-invariant.md).
Transformation laws under q_eta' = g_eta(x) q_eta, g_eta(0) = 1 (same records):
    M' = g_R M g_L^dagger,
    C_eta' = g_eta C_eta g_eta^dagger + i g_eta d_x g_eta^dagger,
    Theta' = Theta - (1/(2i)) log(det g_R / det g_L),
    Theta_gamma' = Theta_gamma - 2 N_c sum_b Q_b^2 (1/(2i)) log(det_b g_R / det_b g_L)
                 = Theta_gamma - N_c tr[Q^2 (G_R - G_L)]  for g = exp(i G), [g, Q] = 0.
The first-order Alda:2026nrd parametrization (alp-dressed-mass-first-order-invariants.md):
    S_0 = m_q, P_0 = 0, S_1 = i[c^S, m_q], P_1 = {c^P, m_q}, C_{V,0} = c^V, C_{A,0} = c^A.

Derivations (noncommuting SymPy symbols)
----------------------------------------
D1  W_{eta,1} = i C_{eta,0}, W_{eta,2} = i C_{eta,1} - C_{eta,0}^2 from the transport ODE.
D2  hat M_2 = M_2 - 2i(C_{R,0}M_1 - M_1C_{L,0}) - i(C_{R,1}M_0 - M_0C_{L,1})
              - (C_{R,0}^2 M_0 - 2C_{R,0}M_0C_{L,0} + M_0C_{L,0}^2);
    hat S_2, hat P_2 written with S_n, P_n and C_{R,n}, C_{L,n} only. Several groupings
    with commutators and anticommutators are verified; the one with fewest terms is kept.
    The same is done with C_{V,n}, C_{A,n} for comparison.
D3  hat Theta_2 = Theta_2 + tr C_{A,1}; hat Theta_{gamma,2} = Theta_{gamma,2} + 2N_c tr(Q^2 C_{A,1}).
D4  Alda:2026nrd specialization of hat S_2, hat P_2.
D5  Real quark masses: S_0 = m_q, P_0 = 0, all other sources generic Hermitian
    (noncommuting). Candidates for hat S_2, hat P_2 in the C_V, C_A and C_R, C_L bases are
    compared with the transport, symbolically and on exact 3x3 matrices (m_q real positive
    diagonal). Rewritings of the hat P_2 cross term
    i([C_V0, {C_A0, m_q}] + {C_A0, [C_V0, m_q]}) via the Leibniz rule
    [X, {Y, Z}] = {[X, Y], Z} + {Y, [X, Z]} are verified (valid for any m_q); a monomial
    count shows that no single nested bracket equals it.

Checks (each prints PASS/FAIL; exit status 1 on any failure)
------------------------------------------------------------
Agreement with the existing transport recursion (not independent):
    hat M_2 equals invariant_mass_coefficients(x, 2)[2] of
    derive_alp_invariant_coefficient_functions.py.
Independent checks:
C1  invariance from the transformation laws only (Cayley-transform g_eta(x), noncommuting,
    tr H != 0; photon with [g, Q] = 0), exact rationals;
C2  constant connections, hat M = exp(-i x C_R) M exp(i x C_L);
C3  numerical integration of the transport ODE with mpmath;
C4  every candidate grouping of hat S_2, hat P_2 on exact random rational 3x3 matrices;
C5  first-order regression (hat S_1, hat P_1, hat Theta_1, hat Theta_{gamma,1});
C6  vacuum angle Theta + (1/2) arg det M through O(x^2);
C7  Hermiticity of hat S_2, hat P_2;
C8  chiral-current forms of hat Theta_2, hat Theta_{gamma,2} and the N_c = 3 flavor
    components of hat Theta_{gamma,2} (generic symbolic 3x3 matrices).

Assumptions and approximations
------------------------------
* Same restricted partonic truncation as the general record (at most one ALP
  derivative, no dipoles, no electroweak fields besides the photon with fixed Q).
* Photon: sources and redefinitions commute with Q (assumptions 1 and 2 of the photon record).
* All series checks are exact through O(x^2) only (Cayley transforms truncated at O(x^3));
  they are not all-order proofs. C3 is a floating-point check with a stated tolerance.
* "Shortest" means fewest top-level terms among the candidate groupings coded below;
  it is not a proof of minimality over all possible groupings.
"""

import random
import sys
from itertools import permutations, product
from math import comb
from pathlib import Path

import mpmath
import sympy as sp
from sympy.physics.quantum import AntiCommutator, Commutator

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_alp_invariant_coefficient_functions import invariant_mass_coefficients  # noqa: E402
from derive_alp_dressed_mass_first_order import random_rational_hermitian  # noqa: E402
from derive_alp_dressed_topological_invariant import (  # noqa: E402
    determinant_series,
    logarithm_series,
    random_rational_complex,
    report,
    taylor_polynomial,
)
from derive_alp_dressed_photon_invariant import (  # noqa: E402
    block_anomaly_shift,
    project_to_charge_commutant,
)

ALP_FIELD_RATIO = sp.Symbol("x", real=True)
COLOR_NUMBER = 3
FLAVOR_NUMBER = 3
CHARGE_MATRIX = sp.diag(sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(-1, 3))
DEFAULT_SEED = 20261112
# Highest derivative kept for the redefinitions g_eta(x) in exact checks.
GROUP_SERIES_ORDER = 3
# Sampling step and relative tolerance of the numerical transport check C3.
NUMERICAL_STEP = sp.Rational(1, 100)
NUMERICAL_TOLERANCE = mpmath.mpf("1e-12")


def hermitian_matrix_symbol(name):
    return sp.Symbol(name, commutative=False, hermitian=True)


# Derivatives at x = 0 of the partonic sources.
RIGHT_CURRENT = [hermitian_matrix_symbol(f"C_{{R,{n}}}") for n in range(2)]
LEFT_CURRENT = [hermitian_matrix_symbol(f"C_{{L,{n}}}") for n in range(2)]
VECTOR_CURRENT = [hermitian_matrix_symbol(f"C_{{V,{n}}}") for n in range(2)]
AXIAL_CURRENT = [hermitian_matrix_symbol(f"C_{{A,{n}}}") for n in range(2)]
SCALAR = [hermitian_matrix_symbol(f"S_{{{n}}}") for n in range(3)]
PSEUDOSCALAR = [hermitian_matrix_symbol(f"P_{{{n}}}") for n in range(3)]
MASS = [SCALAR[n] + sp.I * PSEUDOSCALAR[n] for n in range(3)]
THETA = [sp.Symbol(rf"\Theta_{{{n}}}", real=True) for n in range(3)]
THETA_PHOTON = [sp.Symbol(rf"\Theta_{{\gamma,{n}}}", real=True) for n in range(3)]
AXIAL_TRACE = [sp.Symbol(rf"\operatorname{{tr}}C_{{A,{n}}}", real=True) for n in range(2)]
AXIAL_CHARGE_TRACE = [sp.Symbol(rf"\operatorname{{tr}}(Q^2C_{{A,{n}}})", real=True)
                      for n in range(2)]

# Alda:2026nrd first-order couplings.
QUARK_MASS = hermitian_matrix_symbol("m_q")
VECTOR_COUPLING = hermitian_matrix_symbol("c^{V}")
AXIAL_COUPLING = hermitian_matrix_symbol("c^{A}")
SCALAR_COUPLING = hermitian_matrix_symbol("c^{S}")
PSEUDOSCALAR_COUPLING = hermitian_matrix_symbol("c^{P}")
VECTOR_SCALAR_INVARIANT = hermitian_matrix_symbol("c^{VS}")
AXIAL_PSEUDOSCALAR_INVARIANT = hermitian_matrix_symbol("c^{AP}")


# ---------------------------------------------------------------------------
# Generic algebra (works for noncommuting symbols and for explicit matrices)
# ---------------------------------------------------------------------------

def commutator(first, second):
    return first * second - second * first


def anticommutator(first, second):
    return first * second + second * first


def adjoint(value):
    return value.H if isinstance(value, sp.MatrixBase) else sp.adjoint(value)


def expand_any(value):
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(sp.expand)
    return sp.expand(value)


def is_zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(sp.expand(entry) == 0 for entry in value)
    return sp.expand(value) == 0


def hermitian_and_antihermitian_parts(value):
    """(hat S, hat P) = ((X + X^dagger)/2, (X - X^dagger)/(2i))."""
    value_adjoint = adjoint(value)
    return (expand_any((value + value_adjoint) / 2),
            expand_any((value - value_adjoint) / (2 * sp.I)))


def leibniz_product(first, second, order):
    """Derivatives at 0 of a product: (f g)_n = sum_r binom(n, r) f_r g_{n-r}.

    Entries beyond the given lists are zero.
    """
    result = []
    for degree in range(order + 1):
        terms = [comb(degree, r) * first[r] * second[degree - r] for r in range(degree + 1)
                 if r < len(first) and degree - r < len(second)]
        result.append(expand_any(sum(terms[1:], start=terms[0])))
    return result


def derivatives_of_scalar_polynomial(expression, order):
    expression = sp.expand(expression)
    return [sp.factorial(degree) * expression.coeff(ALP_FIELD_RATIO, degree)
            for degree in range(order + 1)]


# ---------------------------------------------------------------------------
# D1: transport matrices
# ---------------------------------------------------------------------------

def transport_derivatives(current_derivatives, identity, order):
    """W_n for dW/dx = i C W, W(0) = 1.

    Differentiating the ODE n times at 0: W_{n+1} = i sum_{r=0}^{n} binom(n, r) C_r W_{n-r}.
    Derivatives C_r beyond the list are zero.
    """
    transport = [identity]
    for degree in range(order):
        terms = [comb(degree, r) * current_derivatives[r] * transport[degree - r]
                 for r in range(min(degree, len(current_derivatives) - 1) + 1)]
        transport.append(expand_any(sp.I * sum(terms[1:], start=terms[0])))
    return transport


def dressed_mass_derivatives(right_current, left_current, mass, order):
    """hat M_n = derivatives of W_R^dagger M W_L (general Leibniz rule)."""
    identity = sp.eye(mass[0].shape[0]) if isinstance(mass[0], sp.MatrixBase) else sp.Integer(1)
    right_adjoint = [adjoint(w) for w in transport_derivatives(right_current, identity, order)]
    left_transport = transport_derivatives(left_current, identity, order)
    return leibniz_product(leibniz_product(right_adjoint, mass, order), left_transport, order)


def derive_transport():
    """D1: W_eta,1, W_eta,2, ODE residuals at x^0, x^1 and unitarity through x^2."""
    results = {}
    for chirality, current in (("R", RIGHT_CURRENT), ("L", LEFT_CURRENT)):
        transport = transport_derivatives(current, sp.Integer(1), 2)
        transport_polynomial = taylor_polynomial(transport)
        current_polynomial = current[0] + ALP_FIELD_RATIO * current[1]
        ode_residual = sp.expand(sp.diff(transport_polynomial, ALP_FIELD_RATIO)
                                 - sp.I * current_polynomial * transport_polynomial)
        unitarity = sp.expand(sp.adjoint(transport_polynomial) * transport_polynomial)
        results[chirality] = {
            "transport": transport,
            "ode_residuals": [ode_residual.coeff(ALP_FIELD_RATIO, degree) for degree in range(2)],
            "unitarity_residuals": [sp.expand(unitarity.coeff(ALP_FIELD_RATIO, degree)
                                              - (1 if degree == 0 else 0))
                                    for degree in range(3)],
            "closed_form_residuals": [
                sp.expand(transport[1] - sp.I * current[0]),
                sp.expand(transport[2] - (sp.I * current[1] - current[0]**2))],
        }
    return results


# ---------------------------------------------------------------------------
# D2: closed form of hat M_2 and groupings of hat S_2, hat P_2
# ---------------------------------------------------------------------------

def right_minus_left_multiplication(right_current, left_current, matrix):
    """C_R X - X C_L."""
    return right_current * matrix - matrix * left_current


def closed_form_dressed_mass_one(mass, right_current, left_current):
    """hat M_1 = M_1 - i (C_{R,0} M_0 - M_0 C_{L,0})."""
    return expand_any(mass[1] - sp.I * right_minus_left_multiplication(
        right_current[0], left_current[0], mass[0]))


def closed_form_dressed_mass_two(mass, right_current, left_current):
    """hat M_2 = M_2 - 2i(C_{R,0}M_1 - M_1C_{L,0}) - i(C_{R,1}M_0 - M_0C_{L,1})
                 - (C_{R,0}^2 M_0 - 2 C_{R,0} M_0 C_{L,0} + M_0 C_{L,0}^2)."""
    first_order_dressing = right_minus_left_multiplication(
        right_current[0], left_current[0], mass[0])
    return expand_any(
        mass[2]
        - 2 * sp.I * right_minus_left_multiplication(right_current[0], left_current[0], mass[1])
        - sp.I * right_minus_left_multiplication(right_current[1], left_current[1], mass[0])
        - right_minus_left_multiplication(right_current[0], left_current[0],
                                          first_order_dressing))


# Each grouping below takes the derivative lists of S, P and of (C_R, C_L) or (C_V, C_A),
# plus the bracket functions, and returns (hat S_2, hat P_2).

def grouping_single_current_brackets(s, p, right, left, bracket, anti_bracket):
    """Chiral basis; each bracket contains C_R or C_L alone, or C_R^2 +- C_L^2."""
    r0, r1, l0, l1 = right[0], right[1], left[0], left[1]
    half, i_half = sp.Rational(1, 2), sp.I / 2
    scalar = (s[2]
              - sp.I * bracket(r0, s[1]) - sp.I * bracket(l0, s[1])
              + anti_bracket(r0, p[1]) - anti_bracket(l0, p[1])
              - i_half * bracket(r1, s[0]) - i_half * bracket(l1, s[0])
              + half * anti_bracket(r1, p[0]) - half * anti_bracket(l1, p[0])
              - half * bracket(r0, bracket(l0, s[0])) + half * anti_bracket(r0, anti_bracket(l0, s[0]))
              - half * anti_bracket(r0**2 + l0**2, s[0])
              + i_half * bracket(r0, anti_bracket(l0, p[0]))
              - i_half * anti_bracket(r0, bracket(l0, p[0]))
              - i_half * bracket(r0**2 - l0**2, p[0]))
    pseudoscalar = (p[2]
                    - anti_bracket(r0, s[1]) + anti_bracket(l0, s[1])
                    - sp.I * bracket(r0, p[1]) - sp.I * bracket(l0, p[1])
                    - half * anti_bracket(r1, s[0]) + half * anti_bracket(l1, s[0])
                    - i_half * bracket(r1, p[0]) - i_half * bracket(l1, p[0])
                    - i_half * bracket(r0, anti_bracket(l0, s[0]))
                    + i_half * anti_bracket(r0, bracket(l0, s[0]))
                    + i_half * bracket(r0**2 - l0**2, s[0])
                    - half * bracket(r0, bracket(l0, p[0]))
                    + half * anti_bracket(r0, anti_bracket(l0, p[0]))
                    - half * anti_bracket(r0**2 + l0**2, p[0]))
    return scalar, pseudoscalar


def grouping_chiral_by_order(s, p, right, left, bracket, anti_bracket):
    """Chiral basis; one bracket per S_n or P_n structure, C_R +- C_L inside the brackets."""
    sum0, difference0 = right[0] + left[0], right[0] - left[0]
    sum1, difference1 = right[1] + left[1], right[1] - left[1]
    quarter, i_quarter = sp.Rational(1, 4), sp.I / 4
    scalar = (s[2] - sp.I * bracket(sum0, s[1]) + anti_bracket(difference0, p[1])
              - sp.I / 2 * bracket(sum1, s[0]) + anti_bracket(difference1, p[0]) / 2
              - quarter * bracket(sum0, bracket(sum0, s[0]))
              - quarter * anti_bracket(difference0, anti_bracket(difference0, s[0]))
              - i_quarter * bracket(sum0, anti_bracket(difference0, p[0]))
              - i_quarter * anti_bracket(difference0, bracket(sum0, p[0])))
    pseudoscalar = (p[2] - sp.I * bracket(sum0, p[1]) - anti_bracket(difference0, s[1])
                    - sp.I / 2 * bracket(sum1, p[0]) - anti_bracket(difference1, s[0]) / 2
                    - quarter * bracket(sum0, bracket(sum0, p[0]))
                    - quarter * anti_bracket(difference0, anti_bracket(difference0, p[0]))
                    + i_quarter * bracket(sum0, anti_bracket(difference0, s[0]))
                    + i_quarter * anti_bracket(difference0, bracket(sum0, s[0])))
    return scalar, pseudoscalar


def grouping_chiral_by_order_mixed(s, p, right, left, bracket, anti_bracket):
    """As grouping_chiral_by_order, but the terms with C_{.,0}^2 and P_0 (in hat S_2) or
    S_0 (in hat P_2) use {C_{R,0}, [C_{R,0} + C_{L,0}, X]} - [C_{R,0} + C_{L,0}, {C_{L,0}, X}]."""
    r0, l0 = right[0], left[0]
    sum0, difference0 = r0 + l0, r0 - l0
    sum1, difference1 = right[1] + left[1], right[1] - left[1]
    quarter = sp.Rational(1, 4)
    scalar = (s[2] - sp.I * bracket(sum0, s[1]) + anti_bracket(difference0, p[1])
              - sp.I / 2 * bracket(sum1, s[0]) + anti_bracket(difference1, p[0]) / 2
              - quarter * bracket(sum0, bracket(sum0, s[0]))
              - quarter * anti_bracket(difference0, anti_bracket(difference0, s[0]))
              - sp.I / 2 * anti_bracket(r0, bracket(sum0, p[0]))
              + sp.I / 2 * bracket(sum0, anti_bracket(l0, p[0])))
    pseudoscalar = (p[2] - sp.I * bracket(sum0, p[1]) - anti_bracket(difference0, s[1])
                    - sp.I / 2 * bracket(sum1, p[0]) - anti_bracket(difference1, s[0]) / 2
                    - quarter * bracket(sum0, bracket(sum0, p[0]))
                    - quarter * anti_bracket(difference0, anti_bracket(difference0, p[0]))
                    + sp.I / 2 * anti_bracket(r0, bracket(sum0, s[0]))
                    - sp.I / 2 * bracket(sum0, anti_bracket(l0, s[0])))
    return scalar, pseudoscalar


def grouping_chiral_nested(s, p, right, left, bracket, anti_bracket):
    """Chiral basis, nested: the C_{.,0} terms acting on S_1, P_1 and on the first-order
    dressing of S_0, P_0 are collected inside one commutator and one anticommutator.

    hat S_2 = S_2 - (i/4)[C_R0 + C_L0, 4S_1 - i[C_R0 + C_L0, S_0] + {C_R0 - C_L0, P_0}]
                  + (1/4){C_R0 - C_L0, 4P_1 - i[C_R0 + C_L0, P_0] - {C_R0 - C_L0, S_0}}
                  - (i/2)[C_R1 + C_L1, S_0] + (1/2){C_R1 - C_L1, P_0},
    hat P_2 = P_2 - (i/4)[C_R0 + C_L0, 4P_1 - i[C_R0 + C_L0, P_0] - {C_R0 - C_L0, S_0}]
                  - (1/4){C_R0 - C_L0, 4S_1 - i[C_R0 + C_L0, S_0] + {C_R0 - C_L0, P_0}}
                  - (i/2)[C_R1 + C_L1, P_0] - (1/2){C_R1 - C_L1, S_0}.
    """
    sum0, difference0 = right[0] + left[0], right[0] - left[0]
    sum1, difference1 = right[1] + left[1], right[1] - left[1]
    scalar_inner = 4 * s[1] - sp.I * bracket(sum0, s[0]) + anti_bracket(difference0, p[0])
    pseudoscalar_inner = 4 * p[1] - sp.I * bracket(sum0, p[0]) - anti_bracket(difference0, s[0])
    scalar = (s[2] - sp.I / 4 * bracket(sum0, scalar_inner)
              + anti_bracket(difference0, pseudoscalar_inner) / 4
              - sp.I / 2 * bracket(sum1, s[0]) + anti_bracket(difference1, p[0]) / 2)
    pseudoscalar = (p[2] - sp.I / 4 * bracket(sum0, pseudoscalar_inner)
                    - anti_bracket(difference0, scalar_inner) / 4
                    - sp.I / 2 * bracket(sum1, p[0]) - anti_bracket(difference1, s[0]) / 2)
    return scalar, pseudoscalar


def grouping_vector_axial_by_order(s, p, vector, axial, bracket, anti_bracket):
    """Vector/axial basis, one bracket per structure."""
    v0, v1, a0, a1 = vector[0], vector[1], axial[0], axial[1]
    scalar = (s[2] - 2 * sp.I * bracket(v0, s[1]) + 2 * anti_bracket(a0, p[1])
              - sp.I * bracket(v1, s[0]) + anti_bracket(a1, p[0])
              - bracket(v0, bracket(v0, s[0])) - anti_bracket(a0, anti_bracket(a0, s[0]))
              - sp.I * bracket(v0, anti_bracket(a0, p[0]))
              - sp.I * anti_bracket(a0, bracket(v0, p[0])))
    pseudoscalar = (p[2] - 2 * sp.I * bracket(v0, p[1]) - 2 * anti_bracket(a0, s[1])
                    - sp.I * bracket(v1, p[0]) - anti_bracket(a1, s[0])
                    - bracket(v0, bracket(v0, p[0])) - anti_bracket(a0, anti_bracket(a0, p[0]))
                    + sp.I * bracket(v0, anti_bracket(a0, s[0]))
                    + sp.I * anti_bracket(a0, bracket(v0, s[0])))
    return scalar, pseudoscalar


def grouping_vector_axial_nested(s, p, vector, axial, bracket, anti_bracket):
    """Vector/axial basis, nested:

    hat S_2 = S_2 - i[C_V0, 2S_1 - i[C_V0, S_0] + {C_A0, P_0}]
                  + {C_A0, 2P_1 - i[C_V0, P_0] - {C_A0, S_0}} - i[C_V1, S_0] + {C_A1, P_0},
    hat P_2 = P_2 - i[C_V0, 2P_1 - i[C_V0, P_0] - {C_A0, S_0}]
                  - {C_A0, 2S_1 - i[C_V0, S_0] + {C_A0, P_0}} - i[C_V1, P_0] - {C_A1, S_0}.
    """
    v0, v1, a0, a1 = vector[0], vector[1], axial[0], axial[1]
    scalar_inner = 2 * s[1] - sp.I * bracket(v0, s[0]) + anti_bracket(a0, p[0])
    pseudoscalar_inner = 2 * p[1] - sp.I * bracket(v0, p[0]) - anti_bracket(a0, s[0])
    scalar = (s[2] - sp.I * bracket(v0, scalar_inner) + anti_bracket(a0, pseudoscalar_inner)
              - sp.I * bracket(v1, s[0]) + anti_bracket(a1, p[0]))
    pseudoscalar = (p[2] - sp.I * bracket(v0, pseudoscalar_inner)
                    - anti_bracket(a0, scalar_inner)
                    - sp.I * bracket(v1, p[0]) - anti_bracket(a1, s[0]))
    return scalar, pseudoscalar


CHIRAL_GROUPINGS = {
    "chiral, single-current brackets": grouping_single_current_brackets,
    "chiral, one bracket per structure": grouping_chiral_by_order,
    "chiral, one bracket per structure (mixed quadratic)": grouping_chiral_by_order_mixed,
    "chiral, nested": grouping_chiral_nested,
}
VECTOR_AXIAL_GROUPINGS = {
    "vector/axial, one bracket per structure": grouping_vector_axial_by_order,
    "vector/axial, nested": grouping_vector_axial_nested,
}


def vector_axial_as_chiral(right, left):
    """C_V = (C_R + C_L)/2, C_A = (C_R - C_L)/2."""
    return ([(r + l) / 2 for r, l in zip(right, left)],
            [(r - l) / 2 for r, l in zip(right, left)])


def evaluate_grouping(grouping, scalar, pseudoscalar, right, left, chiral,
                      bracket=commutator, anti_bracket=anticommutator):
    """Apply a grouping to chiral currents; vector/axial groupings get C_V, C_A from C_R, C_L."""
    currents = (right, left) if chiral else vector_axial_as_chiral(right, left)
    return grouping(scalar, pseudoscalar, *currents, bracket, anti_bracket)


def count_terms(expression):
    """(number of top-level terms, number of commutators and anticommutators)."""
    brackets = sum(1 for node in sp.preorder_traversal(expression)
                   if isinstance(node, (Commutator, AntiCommutator)))
    return len(sp.Add.make_args(expression)), brackets


def transport_scalar_pseudoscalar_two():
    """hat S_2, hat P_2 from the D1 transport (symbols)."""
    dressed = dressed_mass_derivatives(RIGHT_CURRENT, LEFT_CURRENT, MASS, 2)
    return dressed, hermitian_and_antihermitian_parts(dressed[2])


def check_groupings(chiral_groupings=None, vector_axial_groupings=None):
    """D2: verify each grouping symbolically, count terms and pick the shortest.

    Residuals are (grouping - transport) for hat S_2 and hat P_2 in the noncommuting
    symbols S_n, P_n, C_{R,n}, C_{L,n}; vector/axial groupings are compared after
    C_V = (C_R + C_L)/2, C_A = (C_R - C_L)/2.
    """
    chiral_groupings = CHIRAL_GROUPINGS if chiral_groupings is None else chiral_groupings
    vector_axial_groupings = (VECTOR_AXIAL_GROUPINGS if vector_axial_groupings is None
                              else vector_axial_groupings)
    _, (transport_scalar, transport_pseudoscalar) = transport_scalar_pseudoscalar_two()
    table = {}
    for chiral, groupings in ((True, chiral_groupings), (False, vector_axial_groupings)):
        for name, grouping in groupings.items():
            scalar, pseudoscalar = evaluate_grouping(grouping, SCALAR, PSEUDOSCALAR,
                                                     RIGHT_CURRENT, LEFT_CURRENT, chiral)
            currents = (RIGHT_CURRENT, LEFT_CURRENT) if chiral else (VECTOR_CURRENT,
                                                                     AXIAL_CURRENT)
            display = grouping(SCALAR, PSEUDOSCALAR, *currents, Commutator, AntiCommutator)
            table[name] = {
                "chiral": chiral,
                "scalar_residual": sp.expand(scalar - transport_scalar),
                "pseudoscalar_residual": sp.expand(pseudoscalar - transport_pseudoscalar),
                "display": display,
                "scalar_count": count_terms(display[0]),
                "pseudoscalar_count": count_terms(display[1]),
            }

    def shortest(chiral):
        verified = [name for name, row in table.items() if row["chiral"] == chiral
                    and row["scalar_residual"] == 0 and row["pseudoscalar_residual"] == 0]
        return min(verified, key=lambda name: (
            table[name]["scalar_count"][0] + table[name]["pseudoscalar_count"][0],
            table[name]["scalar_count"][1] + table[name]["pseudoscalar_count"][1]),
            default=None)

    return {"table": table, "shortest_chiral": shortest(True),
            "shortest_vector_axial": shortest(False)}


def check_existing_transport_agreement():
    """hat M_2 equals the order-2 result of derive_alp_invariant_coefficient_functions.py.

    Both use the derivative convention; this is agreement with the existing transport
    recursion, not an independent check.
    """
    _, _, derivatives = invariant_mass_coefficients(ALP_FIELD_RATIO, 2)
    by_name = {symbol.name: symbol for symbol in derivatives[2].free_symbols}
    substitution = {by_name["C_R0"]: RIGHT_CURRENT[0], by_name["C_R1"]: RIGHT_CURRENT[1],
                    by_name["C_L0"]: LEFT_CURRENT[0], by_name["C_L1"]: LEFT_CURRENT[1],
                    by_name["M_0"]: MASS[0], by_name["M_1"]: MASS[1], by_name["M_2"]: MASS[2]}
    existing = sp.expand(derivatives[2].xreplace(substitution))
    dressed, _ = transport_scalar_pseudoscalar_two()
    closed = closed_form_dressed_mass_two(MASS, RIGHT_CURRENT, LEFT_CURRENT)
    return {"existing_minus_transport": sp.expand(existing - dressed[2]),
            "existing_minus_closed_form": sp.expand(existing - closed),
            "transport_one_residual": sp.expand(dressed[1] - closed_form_dressed_mass_one(
                MASS, RIGHT_CURRENT, LEFT_CURRENT))}


# ---------------------------------------------------------------------------
# D3: hat Theta_2 and hat Theta_{gamma,2}
# ---------------------------------------------------------------------------

def closed_form_theta_two(theta, axial_trace):
    """hat Theta_2 = Theta_2 + tr C_{A,1}."""
    return theta[2] + axial_trace[1]


def closed_form_theta_photon_two(theta_photon, axial_charge_trace):
    """hat Theta_{gamma,2} = Theta_{gamma,2} + 2 N_c tr(Q^2 C_{A,1})."""
    return theta_photon[2] + 2 * COLOR_NUMBER * axial_charge_trace[1]


def derive_theta_two():
    """D3: second derivative of Theta + int_0^x tr C_A (and the photon analogue)."""
    integration_variable = sp.Symbol("t", real=True)
    axial_trace_series = AXIAL_TRACE[0] + AXIAL_TRACE[1] * integration_variable
    axial_charge_series = AXIAL_CHARGE_TRACE[0] + AXIAL_CHARGE_TRACE[1] * integration_variable
    dressed_theta = taylor_polynomial(THETA) + sp.integrate(
        axial_trace_series, (integration_variable, 0, ALP_FIELD_RATIO))
    dressed_theta_photon = taylor_polynomial(THETA_PHOTON) + 2 * COLOR_NUMBER * sp.integrate(
        axial_charge_series, (integration_variable, 0, ALP_FIELD_RATIO))
    theta = derivatives_of_scalar_polynomial(dressed_theta, 2)
    theta_photon = derivatives_of_scalar_polynomial(dressed_theta_photon, 2)
    return {"theta_one": theta[1], "theta_photon_one": theta_photon[1],
            "theta_two": theta[2], "theta_photon_two": theta_photon[2],
            "theta_residual": sp.expand(theta[2] - closed_form_theta_two(THETA, AXIAL_TRACE)),
            "theta_photon_residual": sp.expand(theta_photon[2] - closed_form_theta_photon_two(
                THETA_PHOTON, AXIAL_CHARGE_TRACE))}


# ---------------------------------------------------------------------------
# D4: Alda:2026nrd specialization
# ---------------------------------------------------------------------------

def alda_second_order_sources(vector, axial, scalar_coupling, pseudoscalar_coupling):
    """S_n, P_n, C_{V,n}, C_{A,n} with the Alda:2026nrd parametrization at orders 0 and 1.

    S_0 = m_q, P_0 = 0, S_1 = i[c^S, m_q], P_1 = {c^P, m_q}, C_{V,0} = c^V, C_{A,0} = c^A;
    C_{V,1}, C_{A,1}, S_2, P_2 remain generic.
    """
    scalar = [QUARK_MASS, sp.I * commutator(scalar_coupling, QUARK_MASS), SCALAR[2]]
    pseudoscalar = [sp.Integer(0), anticommutator(pseudoscalar_coupling, QUARK_MASS),
                    PSEUDOSCALAR[2]]
    return scalar, pseudoscalar, [vector, VECTOR_CURRENT[1]], [axial, AXIAL_CURRENT[1]]


def alda_compact_forms(bracket=commutator, anti_bracket=anticommutator):
    """hat S_2 = S_2 - i[C_V1, m_q] + [c^V, [2c^S - c^V, m_q]] + {c^A, {2c^P - c^A, m_q}},
    hat P_2 = P_2 - {C_A1, m_q} - i[c^V, {2c^P - c^A, m_q}] - i{c^A, [2c^S - c^V, m_q]}."""
    mass, vector, axial = QUARK_MASS, VECTOR_COUPLING, AXIAL_COUPLING
    vector_combination = 2 * SCALAR_COUPLING - vector
    axial_combination = 2 * PSEUDOSCALAR_COUPLING - axial
    scalar = (SCALAR[2] - sp.I * bracket(VECTOR_CURRENT[1], mass)
              + bracket(vector, bracket(vector_combination, mass))
              + anti_bracket(axial, anti_bracket(axial_combination, mass)))
    pseudoscalar = (PSEUDOSCALAR[2] - anti_bracket(AXIAL_CURRENT[1], mass)
                    - sp.I * bracket(vector, anti_bracket(axial_combination, mass))
                    - sp.I * anti_bracket(axial, bracket(vector_combination, mass)))
    return scalar, pseudoscalar


def alda_invariant_coupling_forms(bracket=commutator, anti_bracket=anticommutator):
    """Same as alda_compact_forms with 2c^S - c^V = c^S - c^VS, 2c^P - c^A = c^P - c^AP
    and c^V = c^VS + c^S, c^A = c^AP + c^P."""
    mass = QUARK_MASS
    vector = VECTOR_SCALAR_INVARIANT + SCALAR_COUPLING
    axial = AXIAL_PSEUDOSCALAR_INVARIANT + PSEUDOSCALAR_COUPLING
    vector_combination = SCALAR_COUPLING - VECTOR_SCALAR_INVARIANT
    axial_combination = PSEUDOSCALAR_COUPLING - AXIAL_PSEUDOSCALAR_INVARIANT
    scalar = (SCALAR[2] - sp.I * bracket(VECTOR_CURRENT[1], mass)
              + bracket(vector, bracket(vector_combination, mass))
              + anti_bracket(axial, anti_bracket(axial_combination, mass)))
    pseudoscalar = (PSEUDOSCALAR[2] - anti_bracket(AXIAL_CURRENT[1], mass)
                    - sp.I * bracket(vector, anti_bracket(axial_combination, mass))
                    - sp.I * anti_bracket(axial, bracket(vector_combination, mass)))
    return scalar, pseudoscalar


def derive_alda_specialization(compact_forms=alda_compact_forms):
    """D4: insert the Alda sources into the chosen chiral grouping and compare.

    The chiral currents are C_R = C_V + C_A, C_L = C_V - C_A with C_V, C_A from
    alda_second_order_sources. Also returns whether hat S_2, hat P_2 depend on c^S, c^P
    at fixed c^VS = c^V - c^S, c^AP = c^A - c^P.
    """
    scalar, pseudoscalar, vector, axial = alda_second_order_sources(
        VECTOR_COUPLING, AXIAL_COUPLING, SCALAR_COUPLING, PSEUDOSCALAR_COUPLING)
    right = [v + a for v, a in zip(vector, axial)]
    left = [v - a for v, a in zip(vector, axial)]
    chosen_scalar, chosen_pseudoscalar = grouping_chiral_nested(
        scalar, pseudoscalar, right, left, commutator, anticommutator)
    chosen_scalar, chosen_pseudoscalar = sp.expand(chosen_scalar), sp.expand(chosen_pseudoscalar)
    compact_scalar, compact_pseudoscalar = compact_forms()
    to_invariants = {VECTOR_COUPLING: VECTOR_SCALAR_INVARIANT + SCALAR_COUPLING,
                     AXIAL_COUPLING: AXIAL_PSEUDOSCALAR_INVARIANT + PSEUDOSCALAR_COUPLING}
    invariant_scalar, invariant_pseudoscalar = alda_invariant_coupling_forms()
    raw_couplings = (SCALAR_COUPLING, PSEUDOSCALAR_COUPLING)
    scalar_in_invariants = sp.expand(chosen_scalar.xreplace(to_invariants))
    pseudoscalar_in_invariants = sp.expand(chosen_pseudoscalar.xreplace(to_invariants))
    return {
        "display": compact_forms(Commutator, AntiCommutator),
        "display_invariant_couplings": alda_invariant_coupling_forms(Commutator,
                                                                     AntiCommutator),
        "scalar_residual": sp.expand(sp.expand(compact_scalar) - chosen_scalar),
        "pseudoscalar_residual": sp.expand(sp.expand(compact_pseudoscalar)
                                           - chosen_pseudoscalar),
        "scalar_invariant_residual": sp.expand(sp.expand(invariant_scalar)
                                               - scalar_in_invariants),
        "pseudoscalar_invariant_residual": sp.expand(sp.expand(invariant_pseudoscalar)
                                                     - pseudoscalar_in_invariants),
        "scalar_depends_on_raw": any(scalar_in_invariants.has(c) for c in raw_couplings),
        "pseudoscalar_depends_on_raw": any(pseudoscalar_in_invariants.has(c)
                                           for c in raw_couplings),
    }


# ---------------------------------------------------------------------------
# Exact rational test data
# ---------------------------------------------------------------------------

def scalar_pseudoscalar_lists(mass):
    parts = [hermitian_and_antihermitian_parts(coefficient) for coefficient in mass]
    return [part[0] for part in parts], [part[1] for part in parts]


def all_invariants_from_derivatives(right, left, mass, theta, theta_photon,
                                    grouping=None):
    """(hat S_2, hat P_2, hat Theta_2, hat Theta_{gamma,2}) from derivatives at x = 0.

    hat S_2, hat P_2 use the chosen chiral grouping unless another is given.
    """
    grouping = grouping_chiral_nested if grouping is None else grouping
    scalar, pseudoscalar = scalar_pseudoscalar_lists(mass)
    dressed_scalar, dressed_pseudoscalar = grouping(scalar, pseudoscalar, right, left,
                                                    commutator, anticommutator)
    _, axial = vector_axial_as_chiral(right, left)
    return {"scalar_two": expand_any(dressed_scalar),
            "pseudoscalar_two": expand_any(dressed_pseudoscalar),
            "theta_two": sp.expand(closed_form_theta_two(
                theta, [coefficient.trace() for coefficient in axial])),
            "theta_photon_two": sp.expand(closed_form_theta_photon_two(
                theta_photon, [(CHARGE_MATRIX**2 * coefficient).trace()
                               for coefficient in axial]))}


def random_rational(generator):
    return sp.Rational(generator.randint(-9, 9), generator.randint(1, 5))


def random_sources(seed, charge_conserving=False, size=FLAVOR_NUMBER):
    """Derivatives 0..2 of C_R, C_L, M, Theta and Theta_gamma (exact rationals)."""
    generator = random.Random(seed)
    project = project_to_charge_commutant if charge_conserving else (lambda matrix: matrix)
    right = [project(random_rational_hermitian(generator, size)) for _ in range(3)]
    left = [project(random_rational_hermitian(generator, size)) for _ in range(3)]
    mass = [project(random_rational_complex(generator, size)) for _ in range(3)]
    theta = [random_rational(generator) for _ in range(3)]
    theta_photon = [random_rational(generator) for _ in range(3)]
    return right, left, mass, theta, theta_photon


def cayley_group_derivatives(hermitian, order=GROUP_SERIES_ORDER):
    """g = (1 - iH)^{-1}(1 + iH) = 1 + 2 sum_{k>=1} (iH)^k, H(0) = 0, derivatives through order."""
    size = hermitian[1].shape[0]
    zero = sp.zeros(size, size)
    step = [sp.I * coefficient for coefficient in hermitian[:order + 1]]
    step += [zero] * (order + 1 - len(step))
    group = [sp.eye(size)] + [zero] * order
    power = [sp.eye(size)] + [zero] * order
    for _ in range(order):
        power = leibniz_product(power, step, order)
        group = [expand_any(group[n] + 2 * power[n]) for n in range(order + 1)]
    return group


def log_determinant_derivatives(group, order=GROUP_SERIES_ORDER):
    """Derivatives at 0 of log det g(x), continuous branch with log det g(0) = 0."""
    determinant = determinant_series(taylor_polynomial(group), order)
    return derivatives_of_scalar_polynomial(logarithm_series(determinant, order), order)


def random_redefinition_generators(seed, charge_conserving=False, size=FLAVOR_NUMBER):
    """Hermitian H_eta with derivatives (0, H_1, H_2), different for R and L, tr H_1 != 0."""
    generator = random.Random(seed)
    generators = {}
    for chirality in ("right", "left"):
        coefficients = [sp.zeros(size, size)]
        for _ in range(2):
            matrix = random_rational_hermitian(generator, size)
            if charge_conserving:
                matrix = project_to_charge_commutant(matrix)
            coefficients.append(matrix)
        generators[chirality] = coefficients
    return generators["right"], generators["left"]


def transform_sources(right, left, mass, theta, theta_photon, right_group, left_group,
                      order=2):
    """Apply the transformation laws only, as derivatives at 0 through `order`."""
    def transformed_current(group, current):
        group_adjoint = [coefficient.H for coefficient in group]
        derivative_adjoint = group_adjoint[1:]
        homogeneous = leibniz_product(leibniz_product(group, current, order),
                                      group_adjoint, order)
        inhomogeneous = leibniz_product(group, derivative_adjoint, order)
        return [expand_any(homogeneous[n] + sp.I * inhomogeneous[n]) for n in range(order + 1)]

    left_adjoint = [coefficient.H for coefficient in left_group]
    new_mass = leibniz_product(leibniz_product(right_group, mass, order), left_adjoint, order)
    log_right = log_determinant_derivatives(right_group)
    log_left = log_determinant_derivatives(left_group)
    new_theta = [sp.expand(theta[n] - (log_right[n] - log_left[n]) / (2 * sp.I))
                 for n in range(order + 1)]
    photon_shift = derivatives_of_scalar_polynomial(
        block_anomaly_shift(taylor_polynomial(right_group), taylor_polynomial(left_group)),
        order)
    new_theta_photon = [sp.expand(theta_photon[n] - photon_shift[n]) for n in range(order + 1)]
    return (transformed_current(right_group, right), transformed_current(left_group, left),
            new_mass, new_theta, new_theta_photon)


# ---------------------------------------------------------------------------
# C1: invariance from the transformation laws
# ---------------------------------------------------------------------------

def check_invariance(seed=DEFAULT_SEED, redefinition_seed=7, charge_conserving=False):
    """C1: chosen invariants unchanged under a Cayley-transform redefinition."""
    right, left, mass, theta, theta_photon = random_sources(seed, charge_conserving)
    right_generator, left_generator = random_redefinition_generators(redefinition_seed,
                                                                     charge_conserving)
    right_group = cayley_group_derivatives(right_generator)
    left_group = cayley_group_derivatives(left_generator)
    identity, zero = sp.eye(FLAVOR_NUMBER), sp.zeros(FLAVOR_NUMBER, FLAVOR_NUMBER)
    unitarity = all(is_zero(term - (identity if n == 0 else zero))
                    for group in (right_group, left_group)
                    for n, term in enumerate(leibniz_product(
                        group, [c.H for c in group], GROUP_SERIES_ORDER)))
    transformed = transform_sources(right, left, mass, theta, theta_photon,
                                    right_group, left_group)
    before = all_invariants_from_derivatives(right, left, mass, theta, theta_photon)
    after = all_invariants_from_derivatives(*transformed)
    trace_form_residual = None
    if charge_conserving:
        # Trace form of the photon anomaly: G = 2 arctan H = 2 (H - H^3/3) + O(x^5).
        def angle_derivatives(hermitian):
            cube = leibniz_product(leibniz_product(hermitian, hermitian, 3), hermitian, 3)
            padded = list(hermitian) + [zero] * (4 - len(hermitian))
            return [expand_any(2 * padded[n] - sp.Rational(2, 3) * cube[n]) for n in range(4)]

        right_angle = angle_derivatives(right_generator)
        left_angle = angle_derivatives(left_generator)
        trace_shift = [sp.expand(COLOR_NUMBER * (CHARGE_MATRIX**2
                                                 * (right_angle[n] - left_angle[n])).trace())
                       for n in range(3)]
        trace_form_residual = [sp.expand(theta_photon[n] - transformed[4][n] - trace_shift[n])
                               for n in range(3)]
    return {
        "group_is_unitary": unitarity,
        "generators_noncommuting": not is_zero(commutator(right_generator[1],
                                                          right_generator[2])),
        "generator_traces_nonzero": all(generator[1].trace() != 0
                                        for generator in (right_generator, left_generator)),
        "raw_scalar_changes": not is_zero(
            hermitian_and_antihermitian_parts(transformed[2][2])[0]
            - hermitian_and_antihermitian_parts(mass[2])[0]),
        "raw_theta_changes": sp.expand(transformed[3][2] - theta[2]) != 0,
        "raw_theta_photon_changes": sp.expand(transformed[4][2] - theta_photon[2]) != 0,
        "residuals": {key: expand_any(after[key] - before[key]) for key in before},
        "photon_trace_form_residual": trace_form_residual,
    }


# ---------------------------------------------------------------------------
# C2: constant connections
# ---------------------------------------------------------------------------

def check_constant_connections(seed=DEFAULT_SEED + 1):
    """C2: hat M = exp(-i x C_R) M exp(i x C_L) for constant C_eta; compare the 2nd derivative."""
    right, left, mass, _, _ = random_sources(seed)
    order = 2

    def exponential_derivatives(generator):
        return [generator**n for n in range(order + 1)]

    dressed = leibniz_product(leibniz_product(exponential_derivatives(-sp.I * right[0]),
                                              mass, order),
                              exponential_derivatives(sp.I * left[0]), order)
    zero = sp.zeros(FLAVOR_NUMBER, FLAVOR_NUMBER)
    constant_right, constant_left = [right[0], zero], [left[0], zero]
    scalar, pseudoscalar = scalar_pseudoscalar_lists(mass)
    chosen = grouping_chiral_nested(scalar, pseudoscalar, constant_right, constant_left,
                                    commutator, anticommutator)
    exponential_scalar, exponential_pseudoscalar = hermitian_and_antihermitian_parts(dressed[2])
    return {"closed_form_residual": expand_any(dressed[2] - closed_form_dressed_mass_two(
                mass, constant_right, constant_left)),
            "chosen_scalar_residual": expand_any(chosen[0] - exponential_scalar),
            "chosen_pseudoscalar_residual": expand_any(chosen[1] - exponential_pseudoscalar)}


# ---------------------------------------------------------------------------
# C3: numerical transport integration
# ---------------------------------------------------------------------------

def check_numerical_transport(seed=DEFAULT_SEED + 2, precision=30, step=NUMERICAL_STEP,
                              point_number=12, tolerance=NUMERICAL_TOLERANCE):
    """C3: integrate dW/dx = i (C_0 + x C_1) W with mpmath.odefun; fit the 2nd derivative.

    hat M(x) is sampled at x_k = k * step (k = 0..point_number); an interpolating
    polynomial of degree point_number is fitted and 2 times its x^2 coefficient is
    compared with the closed-form hat M_2. Source matrices are scaled by 1/5. The fit
    error is the interpolation remainder, of order step^(point_number - 1).
    """
    right, left, mass, _, _ = random_sources(seed)
    scale = sp.Rational(1, 5)
    right = [scale * coefficient for coefficient in right[:2]]
    left = [scale * coefficient for coefficient in left[:2]]
    size = FLAVOR_NUMBER
    with mpmath.workdps(precision):
        def to_mp(matrix):
            def exact_to_mp(value):
                rational = sp.Rational(value)
                return mpmath.mpf(rational.p) / rational.q

            return mpmath.matrix([[mpmath.mpc(exact_to_mp(sp.re(matrix[r, c])),
                                              exact_to_mp(sp.im(matrix[r, c])))
                                   for c in range(size)] for r in range(size)])

        def transport_solution(current):
            current_zero, current_one = to_mp(current[0]), to_mp(current[1])

            def derivative(position, flattened):
                transport = mpmath.matrix(size, size)
                for index, value in enumerate(flattened):
                    transport[index // size, index % size] = value
                product = 1j * (current_zero + position * current_one) * transport
                return [product[index // size, index % size] for index in range(size * size)]

            identity = [mpmath.mpc(1 if index // size == index % size else 0)
                        for index in range(size * size)]
            return mpmath.odefun(derivative, 0, identity)

        right_solution, left_solution = transport_solution(right), transport_solution(left)
        mass_mp = [to_mp(coefficient) for coefficient in mass]
        points = [mpmath.mpf(k) * mpmath.mpf(step.p) / step.q for k in range(point_number + 1)]
        samples = []
        for point in points:
            right_transport = mpmath.matrix(size, size)
            left_transport = mpmath.matrix(size, size)
            for index, (right_value, left_value) in enumerate(zip(right_solution(point),
                                                                  left_solution(point))):
                right_transport[index // size, index % size] = right_value
                left_transport[index // size, index % size] = left_value
            mass_value = mass_mp[0] + point * mass_mp[1] + point**2 / 2 * mass_mp[2]
            samples.append(right_transport.transpose_conj() * mass_value * left_transport)
        vandermonde = mpmath.matrix([[point**power for power in range(point_number + 1)]
                                     for point in points])
        fitted = mpmath.matrix(size, size)
        for row in range(size):
            for column in range(size):
                values = mpmath.matrix([sample[row, column] for sample in samples])
                fitted[row, column] = 2 * mpmath.lu_solve(vandermonde, values)[2]
        closed_mp = to_mp(closed_form_dressed_mass_two(mass, right, left))
        deviation = max(abs(fitted[r, c] - closed_mp[r, c])
                        for r in range(size) for c in range(size))
        scale_of_result = max(abs(closed_mp[r, c]) for r in range(size) for c in range(size))
        return {"max_deviation": deviation, "relative_deviation": deviation / scale_of_result,
                "tolerance": tolerance, "passed": deviation / scale_of_result < tolerance}


# ---------------------------------------------------------------------------
# C4: groupings on exact random matrices
# ---------------------------------------------------------------------------

def check_groupings_on_matrices(seed=DEFAULT_SEED + 5, groupings=None):
    """C4: every grouping against the matrix transport hat M_2 (exact rationals, 3x3)."""
    if groupings is None:
        groupings = {name: (grouping, True) for name, grouping in CHIRAL_GROUPINGS.items()}
        groupings.update({name: (grouping, False)
                          for name, grouping in VECTOR_AXIAL_GROUPINGS.items()})
    right, left, mass, _, _ = random_sources(seed)
    transport_scalar, transport_pseudoscalar = hermitian_and_antihermitian_parts(
        dressed_mass_derivatives(right[:2], left[:2], mass, 2)[2])
    scalar, pseudoscalar = scalar_pseudoscalar_lists(mass)
    passed = {}
    for name, (grouping, chiral) in groupings.items():
        grouped_scalar, grouped_pseudoscalar = evaluate_grouping(
            grouping, scalar, pseudoscalar, right[:2], left[:2], chiral)
        passed[name] = (is_zero(grouped_scalar - transport_scalar)
                        and is_zero(grouped_pseudoscalar - transport_pseudoscalar))
    return passed


# ---------------------------------------------------------------------------
# C5: first-order regression
# ---------------------------------------------------------------------------

def check_first_order_regression():
    """C5: hat S_1 = -i[c^VS, m], hat P_1 = -{c^AP, m}; hat Theta_1, hat Theta_{gamma,1}."""
    scalar, pseudoscalar, vector, axial = alda_second_order_sources(
        VECTOR_COUPLING, AXIAL_COUPLING, SCALAR_COUPLING, PSEUDOSCALAR_COUPLING)
    mass = [scalar[n] + sp.I * pseudoscalar[n] for n in range(2)]
    dressed_one = dressed_mass_derivatives([vector[0] + axial[0]], [vector[0] - axial[0]],
                                           mass, 1)[1]
    scalar_one, pseudoscalar_one = hermitian_and_antihermitian_parts(dressed_one)
    vs = VECTOR_COUPLING - SCALAR_COUPLING
    ap = AXIAL_COUPLING - PSEUDOSCALAR_COUPLING
    theta = derive_theta_two()
    gluon_coupling, photon_coupling = sp.symbols("c_G c_gamma", real=True)
    axial_trace_zero = sp.Symbol(r"\operatorname{tr}c^A", real=True)
    axial_charge_zero = sp.Symbol(r"\operatorname{tr}(Q^2c^A)", real=True)
    substitution = {THETA[1]: gluon_coupling, THETA_PHOTON[1]: photon_coupling,
                    AXIAL_TRACE[0]: axial_trace_zero, AXIAL_CHARGE_TRACE[0]: axial_charge_zero}
    return {
        "scalar_one_residual": sp.expand(scalar_one + sp.I * commutator(vs, QUARK_MASS)),
        "pseudoscalar_one_residual": sp.expand(pseudoscalar_one + anticommutator(ap, QUARK_MASS)),
        "theta_one_residual": sp.expand(theta["theta_one"].xreplace(substitution)
                                        - gluon_coupling - axial_trace_zero),
        "theta_photon_one_residual": sp.expand(
            theta["theta_photon_one"].xreplace(substitution)
            - photon_coupling - 2 * COLOR_NUMBER * axial_charge_zero),
    }


# ---------------------------------------------------------------------------
# C6: vacuum angle
# ---------------------------------------------------------------------------

def imaginary_log_determinant(mass_derivatives, order=2):
    """Derivatives at 0 of Im log det M(x) through `order`, with M_0 positive diagonal."""
    determinant = determinant_series(taylor_polynomial(mass_derivatives), order)
    normalized = sp.expand(determinant / determinant.coeff(ALP_FIELD_RATIO, 0))
    logarithm = derivatives_of_scalar_polynomial(logarithm_series(normalized, order), order)
    return [sp.im(coefficient) for coefficient in logarithm]


def check_vacuum_angle(seed=DEFAULT_SEED + 3):
    """C6: hat Theta + (1/2) arg det hat M = Theta + (1/2) arg det M at O(x) and O(x^2)."""
    generator = random.Random(seed)
    right, left, mass, theta, _ = random_sources(seed)
    mass[0] = sp.diag(*[sp.Rational(generator.randint(1, 9), generator.randint(1, 5))
                        for _ in range(FLAVOR_NUMBER)])
    _, axial = vector_axial_as_chiral(right, left)
    scalar, pseudoscalar = scalar_pseudoscalar_lists(mass)
    dressed_scalar_two, dressed_pseudoscalar_two = grouping_chiral_nested(
        scalar, pseudoscalar, right[:2], left[:2], commutator, anticommutator)
    dressed = [mass[0], closed_form_dressed_mass_one(mass, right, left),
               expand_any(dressed_scalar_two + sp.I * dressed_pseudoscalar_two)]
    dressed_theta = [theta[0], sp.expand(theta[1] + axial[0].trace()),
                     sp.expand(closed_form_theta_two(theta, [c.trace() for c in axial]))]
    phase = imaginary_log_determinant(mass)
    dressed_phase = imaginary_log_determinant(dressed)
    return {"residuals": [sp.simplify(dressed_theta[n] + dressed_phase[n] / 2
                                      - theta[n] - phase[n] / 2) for n in (1, 2)],
            "second_order_phase_changes": sp.simplify(dressed_phase[2] - phase[2]) != 0}


# ---------------------------------------------------------------------------
# C7: Hermiticity
# ---------------------------------------------------------------------------

def check_hermiticity(seed=DEFAULT_SEED + 4):
    """C7: chosen hat S_2 and hat P_2 Hermitian, symbolically and on exact matrices."""
    scalar_two, pseudoscalar_two = (sp.expand(value) for value in grouping_chiral_nested(
        SCALAR, PSEUDOSCALAR, RIGHT_CURRENT, LEFT_CURRENT, commutator, anticommutator))
    right, left, mass, theta, theta_photon = random_sources(seed)
    numeric = all_invariants_from_derivatives(right, left, mass, theta, theta_photon)
    return {"symbolic_scalar": sp.expand(sp.adjoint(scalar_two) - scalar_two),
            "symbolic_pseudoscalar": sp.expand(sp.adjoint(pseudoscalar_two) - pseudoscalar_two),
            "numeric_scalar": is_zero(numeric["scalar_two"].H - numeric["scalar_two"]),
            "numeric_pseudoscalar": is_zero(numeric["pseudoscalar_two"].H
                                            - numeric["pseudoscalar_two"])}


# ---------------------------------------------------------------------------
# C8: chiral-current form of hat Theta_2 and hat Theta_{gamma,2}
# ---------------------------------------------------------------------------

def symbolic_hermitian_matrix(label, size=FLAVOR_NUMBER):
    """Generic Hermitian matrix with real diagonal and complex off-diagonal symbols."""
    matrix = sp.zeros(size, size)
    for row in range(size):
        matrix[row, row] = sp.Symbol(f"{label}_{row}{row}", real=True)
        for column in range(row + 1, size):
            entry = sp.Symbol(f"{label}_{row}{column}")
            matrix[row, column] = entry
            matrix[column, row] = sp.conjugate(entry)
    return matrix


def check_chiral_trace_forms(theta_coefficient=sp.Rational(1, 2),
                             photon_coefficient=sp.Integer(COLOR_NUMBER),
                             component_coefficients=(sp.Rational(8, 3), sp.Rational(2, 3))):
    """C8: chiral-current and component forms of hat Theta_2, hat Theta_{gamma,2} against D3.

    Forms tested (derivative convention):
        hat Theta_2 = Theta_2 + (1/2) tr(C_{R,1} - C_{L,1}),
        hat Theta_{gamma,2} = Theta_{gamma,2} + N_c tr[Q^2 (C_{R,1} - C_{L,1})],
        hat Theta_{gamma,2} = Theta_{gamma,2} + (8/3)(C_{A,1})_uu + (2/3)[(C_{A,1})_dd + (C_{A,1})_ss]
    for N_c = 3, Q = diag(2/3, -1/3, -1/3), generic symbolic Hermitian 3x3 matrices.
    The coefficients are arguments so that the unit tests can check that wrong values fail.
    """
    right_one = symbolic_hermitian_matrix("r")
    left_one = symbolic_hermitian_matrix("l")
    theta = derive_theta_two()

    def d3_with_axial(axial):
        return {AXIAL_TRACE[1]: axial.trace(),
                AXIAL_CHARGE_TRACE[1]: (CHARGE_MATRIX**2 * axial).trace()}

    chiral_substitution = d3_with_axial((right_one - left_one) / 2)
    difference = right_one - left_one
    theta_residual = sp.expand(THETA[2] + theta_coefficient * difference.trace()
                               - theta["theta_two"].xreplace(chiral_substitution))
    photon_residual = sp.expand(
        THETA_PHOTON[2] + photon_coefficient * (CHARGE_MATRIX**2 * difference).trace()
        - theta["theta_photon_two"].xreplace(chiral_substitution))
    generic_axial = symbolic_hermitian_matrix("a")
    up_coefficient, down_coefficient = component_coefficients
    components_form = (THETA_PHOTON[2] + up_coefficient * generic_axial[0, 0]
                       + down_coefficient * (generic_axial[1, 1] + generic_axial[2, 2]))
    components_residual = sp.expand(
        theta["theta_photon_two"].xreplace(d3_with_axial(generic_axial)) - components_form)
    return {"theta_residual": theta_residual,
            "theta_photon_residual": photon_residual,
            "theta_photon_components_residual": components_residual}


# ---------------------------------------------------------------------------
# D5: real quark masses (S_0 = m_q, P_0 = 0, everything else generic)
# ---------------------------------------------------------------------------

def real_mass_vector_axial_forms(mass, scalar, pseudoscalar, vector, axial,
                                 bracket=commutator, anti_bracket=anticommutator,
                                 cross_term=None):
    """hat S_2, hat P_2 for S_0 = m_q, P_0 = 0 (vector/axial basis):

    hat S_2 = S_2 - 2i[C_V0, S_1] + 2{C_A0, P_1} - i[C_V1, m_q]
              - [C_V0, [C_V0, m_q]] - {C_A0, {C_A0, m_q}},
    hat P_2 = P_2 - 2i[C_V0, P_1] - 2{C_A0, S_1} - {C_A1, m_q}
              + i([C_V0, {C_A0, m_q}] + {C_A0, [C_V0, m_q]}).
    scalar, pseudoscalar are the lists (_, S_1, S_2), (_, P_1, P_2); entry 0 is not used.
    cross_term(C_V0, C_A0, m_q, bracket, anti_bracket) replaces the expression multiplying
    i in the last line (default: the form written above).
    """
    cross_term = cross_original if cross_term is None else cross_term
    v0, v1, a0, a1 = vector[0], vector[1], axial[0], axial[1]
    dressed_scalar = (scalar[2] - 2 * sp.I * bracket(v0, scalar[1])
                      + 2 * anti_bracket(a0, pseudoscalar[1]) - sp.I * bracket(v1, mass)
                      - bracket(v0, bracket(v0, mass))
                      - anti_bracket(a0, anti_bracket(a0, mass)))
    dressed_pseudoscalar = (pseudoscalar[2] - 2 * sp.I * bracket(v0, pseudoscalar[1])
                            - 2 * anti_bracket(a0, scalar[1]) - anti_bracket(a1, mass)
                            + sp.I * cross_term(v0, a0, mass, bracket, anti_bracket))
    return dressed_scalar, dressed_pseudoscalar


def real_mass_chiral_forms(mass, scalar, pseudoscalar, right, left,
                           bracket=commutator, anti_bracket=anticommutator):
    """Same as real_mass_vector_axial_forms written with C_R, C_L:

    hat S_2 = S_2 - i[C_R0 + C_L0, S_1] + {C_R0 - C_L0, P_1} - (i/2)[C_R1 + C_L1, m_q]
              - (1/4)[C_R0 + C_L0, [C_R0 + C_L0, m_q]] - (1/4){C_R0 - C_L0, {C_R0 - C_L0, m_q}},
    hat P_2 = P_2 - i[C_R0 + C_L0, P_1] - {C_R0 - C_L0, S_1} - (1/2){C_R1 - C_L1, m_q}
              + (i/2)[C_R0^2 - C_L0^2, m_q] - i(C_R0 m_q C_L0 - C_L0 m_q C_R0).
    """
    sum0, difference0 = right[0] + left[0], right[0] - left[0]
    sum1, difference1 = right[1] + left[1], right[1] - left[1]
    quarter = sp.Rational(1, 4)
    dressed_scalar = (scalar[2] - sp.I * bracket(sum0, scalar[1])
                      + anti_bracket(difference0, pseudoscalar[1])
                      - sp.I / 2 * bracket(sum1, mass)
                      - quarter * bracket(sum0, bracket(sum0, mass))
                      - quarter * anti_bracket(difference0, anti_bracket(difference0, mass)))
    dressed_pseudoscalar = (pseudoscalar[2] - sp.I * bracket(sum0, pseudoscalar[1])
                            - anti_bracket(difference0, scalar[1])
                            - anti_bracket(difference1, mass) / 2
                            + sp.I / 2 * bracket(right[0]**2 - left[0]**2, mass)
                            - sp.I * (right[0] * mass * left[0] - left[0] * mass * right[0]))
    return dressed_scalar, dressed_pseudoscalar


# Equivalent forms of the hat P_2 cross term (without the overall factor i).
def cross_original(vector, axial, mass, bracket, anti_bracket):
    """[C_V0, {C_A0, m_q}] + {C_A0, [C_V0, m_q]}."""
    return bracket(vector, anti_bracket(axial, mass)) + anti_bracket(axial, bracket(vector, mass))


def cross_axial_outside(vector, axial, mass, bracket, anti_bracket):
    """{[C_V0, C_A0], m_q} + 2{C_A0, [C_V0, m_q]}
    (from the Leibniz rule [X, {Y, Z}] = {[X, Y], Z} + {Y, [X, Z]})."""
    return (anti_bracket(bracket(vector, axial), mass)
            + 2 * anti_bracket(axial, bracket(vector, mass)))


def cross_vector_outside(vector, axial, mass, bracket, anti_bracket):
    """2[C_V0, {C_A0, m_q}] - {[C_V0, C_A0], m_q} (same Leibniz rule)."""
    return (2 * bracket(vector, anti_bracket(axial, mass))
            - anti_bracket(bracket(vector, axial), mass))


def cross_with_products(vector, axial, mass, bracket, anti_bracket):
    """[{C_V0, C_A0}, m_q] + 2(C_V0 m_q C_A0 - C_A0 m_q C_V0)."""
    return (bracket(anti_bracket(vector, axial), mass)
            + 2 * (vector * mass * axial - axial * mass * vector))


CROSS_TERM_FORMS = {
    "[C_V0,{C_A0,m}] + {C_A0,[C_V0,m]}": cross_original,
    "{[C_V0,C_A0],m} + 2{C_A0,[C_V0,m]}": cross_axial_outside,
    "2[C_V0,{C_A0,m}] - {[C_V0,C_A0],m}": cross_vector_outside,
    "[{C_V0,C_A0},m] + 2(C_V0 m C_A0 - C_A0 m C_V0)": cross_with_products,
}


def real_mass_transport(mass, scalar, pseudoscalar, vector, axial):
    """hat S_2, hat P_2 from the D1 transport with M_0 = m_q, M_n = S_n + i P_n (n = 1, 2),
    C_R = C_V + C_A, C_L = C_V - C_A."""
    full_mass = [mass] + [scalar[n] + sp.I * pseudoscalar[n] for n in (1, 2)]
    right = [v + a for v, a in zip(vector, axial)]
    left = [v - a for v, a in zip(vector, axial)]
    return hermitian_and_antihermitian_parts(
        dressed_mass_derivatives(right, left, full_mass, 2)[2])


def single_nested_brackets(vector, axial, mass):
    """All b1(X, b2(Y, Z)) with (X, Y, Z) a permutation of (C_V0, C_A0, m_q) and b1, b2
    commutator or anticommutator (as sympy expressions, expanded)."""
    brackets = (commutator, anticommutator)
    return [sp.expand(outer(first, inner(second, third)))
            for first, second, third in permutations((vector, axial, mass))
            for outer, inner in product(brackets, repeat=2)]


def derive_real_mass_specialization(forms=None, cross_forms=None):
    """D5 (symbolic): candidates minus transport, with noncommuting Hermitian symbols.

    forms(mass, scalar, pseudoscalar, vector, axial) returns (hat S_2, hat P_2); default
    real_mass_vector_axial_forms. Also checks the chiral forms and each cross-term form, and
    counts monomials to show that no number times a single nested bracket of C_V0, C_A0, m_q
    equals the cross term.
    """
    forms = real_mass_vector_axial_forms if forms is None else forms
    cross_forms = CROSS_TERM_FORMS if cross_forms is None else cross_forms
    mass, vector, axial = QUARK_MASS, VECTOR_CURRENT, AXIAL_CURRENT
    transport_scalar, transport_pseudoscalar = real_mass_transport(
        mass, SCALAR, PSEUDOSCALAR, vector, axial)
    scalar, pseudoscalar = forms(mass, SCALAR, PSEUDOSCALAR, vector, axial)
    right = [v + a for v, a in zip(vector, axial)]
    left = [v - a for v, a in zip(vector, axial)]
    chiral_scalar, chiral_pseudoscalar = real_mass_chiral_forms(mass, SCALAR, PSEUDOSCALAR,
                                                                right, left)
    cross_reference = sp.expand(cross_original(vector[0], axial[0], mass,
                                               commutator, anticommutator))
    cross_residuals = {}
    for name, cross in cross_forms.items():
        cross_residuals[name] = sp.expand(
            cross(vector[0], axial[0], mass, commutator, anticommutator) - cross_reference)
    # A number times one nested bracket has at most four distinct monomials; if the cross
    # term has more, no single nested bracket can equal it.
    cross_monomial_number = len(sp.Add.make_args(cross_reference))
    largest_single_bracket = max(len(sp.Add.make_args(candidate)) for candidate in
                                 single_nested_brackets(vector[0], axial[0], mass))
    return {
        "scalar_residual": sp.expand(sp.expand(scalar) - transport_scalar),
        "pseudoscalar_residual": sp.expand(sp.expand(pseudoscalar) - transport_pseudoscalar),
        "chiral_scalar_residual": sp.expand(sp.expand(chiral_scalar) - transport_scalar),
        "chiral_pseudoscalar_residual": sp.expand(sp.expand(chiral_pseudoscalar)
                                                  - transport_pseudoscalar),
        "cross_residuals": cross_residuals,
        "cross_monomial_number": cross_monomial_number,
        "largest_single_bracket_monomial_number": largest_single_bracket,
        "transport": (transport_scalar, transport_pseudoscalar),
        "display": forms(mass, SCALAR, PSEUDOSCALAR, vector, axial, Commutator,
                         AntiCommutator),
        "display_chiral": real_mass_chiral_forms(mass, SCALAR, PSEUDOSCALAR, RIGHT_CURRENT,
                                                 LEFT_CURRENT, Commutator, AntiCommutator),
        "display_cross": {name: cross(vector[0], axial[0], mass, Commutator, AntiCommutator)
                          for name, cross in cross_forms.items()},
    }


def check_real_mass_on_matrices(seed=DEFAULT_SEED + 6, forms=None, cross_forms=None):
    """D5 (exact 3x3 rationals): m_q real positive diagonal, all other sources random
    Hermitian; candidates and cross-term forms against the matrix transport."""
    forms = real_mass_vector_axial_forms if forms is None else forms
    cross_forms = CROSS_TERM_FORMS if cross_forms is None else cross_forms
    generator = random.Random(seed)
    size = FLAVOR_NUMBER
    mass = sp.diag(*[sp.Rational(generator.randint(1, 9), generator.randint(1, 5))
                     for _ in range(size)])
    scalar = [None] + [random_rational_hermitian(generator, size) for _ in range(2)]
    pseudoscalar = [None] + [random_rational_hermitian(generator, size) for _ in range(2)]
    vector = [random_rational_hermitian(generator, size) for _ in range(2)]
    axial = [random_rational_hermitian(generator, size) for _ in range(2)]
    transport_scalar, transport_pseudoscalar = real_mass_transport(mass, scalar, pseudoscalar,
                                                                   vector, axial)
    candidate_scalar, candidate_pseudoscalar = forms(mass, scalar, pseudoscalar, vector, axial)
    right = [v + a for v, a in zip(vector, axial)]
    left = [v - a for v, a in zip(vector, axial)]
    chiral_scalar, chiral_pseudoscalar = real_mass_chiral_forms(mass, scalar, pseudoscalar,
                                                                right, left)
    reference = cross_original(vector[0], axial[0], mass, commutator, anticommutator)
    return {
        "noncommuting": not is_zero(commutator(vector[0], axial[0]))
                        and not is_zero(commutator(vector[0], mass)),
        "scalar_residual": expand_any(candidate_scalar - transport_scalar),
        "pseudoscalar_residual": expand_any(candidate_pseudoscalar - transport_pseudoscalar),
        "chiral_scalar_residual": expand_any(chiral_scalar - transport_scalar),
        "chiral_pseudoscalar_residual": expand_any(chiral_pseudoscalar - transport_pseudoscalar),
        "cross_residuals": {name: expand_any(cross(vector[0], axial[0], mass, commutator,
                                                   anticommutator) - reference)
                            for name, cross in cross_forms.items()},
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    results = []
    print("Derivative convention: F_n = d^nF/dx^n at 0, F(x) = sum_n F_n x^n/n!, x = a/f_a.\n")

    transport = derive_transport()
    for chirality, values in transport.items():
        show(f"D1 W_{chirality},1", values["transport"][1])
        show(f"D1 W_{chirality},2", values["transport"][2])
        results.append(report(f"D1 W_{chirality}: ODE residual at x^0, x^1 vanishes; "
                              f"W^dagger W = 1 + O(x^3); W_1 = iC_0, W_2 = iC_1 - C_0^2",
                              all(r == 0 for r in values["ode_residuals"]
                                  + values["unitarity_residuals"]
                                  + values["closed_form_residuals"])))

    agreement = check_existing_transport_agreement()
    show("D2 hat M_2 (closed form)",
         closed_form_dressed_mass_two(MASS, RIGHT_CURRENT, LEFT_CURRENT))
    results.append(report("D2 transport hat M_1 = M_1 - i(C_R0 M_0 - M_0 C_L0)",
                          agreement["transport_one_residual"] == 0))
    results.append(report("D2 agreement with the existing transport recursion "
                          "(invariant_mass_coefficients, order 2): transport and closed form",
                          agreement["existing_minus_transport"] == 0
                          and agreement["existing_minus_closed_form"] == 0))

    groupings = check_groupings()
    for name, row in groupings["table"].items():
        verified = row["scalar_residual"] == 0 and row["pseudoscalar_residual"] == 0
        print(f"D2 grouping '{name}': terms (top-level, brackets) hat S_2 {row['scalar_count']}, "
              f"hat P_2 {row['pseudoscalar_count']}")
        results.append(report(f"D2 grouping '{name}' equals the transport hat S_2, hat P_2",
                              verified))
    for label, key in (("chiral", "shortest_chiral"), ("vector/axial", "shortest_vector_axial")):
        name = groupings[key]
        print(f"\nShortest verified {label} grouping: {name}")
        if name is not None:
            show("hat S_2", groupings["table"][name]["display"][0])
            show("hat P_2", groupings["table"][name]["display"][1])
    results.append(report("D2 shortest chiral grouping is 'chiral, nested'",
                          groupings["shortest_chiral"] == "chiral, nested"))

    theta = derive_theta_two()
    show("D3 hat Theta_2", theta["theta_two"])
    show("D3 hat Theta_{gamma,2}", theta["theta_photon_two"])
    results.append(report("D3 hat Theta_2 = Theta_2 + tr C_A1; "
                          "hat Theta_gamma,2 = Theta_gamma,2 + 2N_c tr(Q^2 C_A1)",
                          theta["theta_residual"] == 0 and theta["theta_photon_residual"] == 0))

    alda = derive_alda_specialization()
    show("D4 hat S_2 (Alda parametrization)", alda["display"][0])
    show("D4 hat P_2 (Alda parametrization)", alda["display"][1])
    show("D4 hat S_2 with c^VS, c^AP", alda["display_invariant_couplings"][0])
    show("D4 hat P_2 with c^VS, c^AP", alda["display_invariant_couplings"][1])
    print(f"D4 hat S_2 depends on c^S, c^P at fixed c^VS, c^AP: {alda['scalar_depends_on_raw']}")
    print(f"D4 hat P_2 depends on c^S, c^P at fixed c^VS, c^AP: "
          f"{alda['pseudoscalar_depends_on_raw']}\n")
    results.append(report("D4 compact forms equal the chosen grouping with Alda sources "
                          "(also in c^VS, c^AP)",
                          alda["scalar_residual"] == 0 and alda["pseudoscalar_residual"] == 0
                          and alda["scalar_invariant_residual"] == 0
                          and alda["pseudoscalar_invariant_residual"] == 0))

    for label, charge_conserving in (("general g", False), ("[g, Q] = 0", True)):
        invariance = check_invariance(charge_conserving=charge_conserving)
        results.append(report(f"C1 ({label}) Cayley g unitary through x^3, [H_1,H_2] != 0, "
                              f"tr H_1 != 0",
                              invariance["group_is_unitary"]
                              and invariance["generators_noncommuting"]
                              and invariance["generator_traces_nonzero"]))
        results.append(report(f"C1 ({label}) raw S_2, Theta_2, Theta_gamma,2 change",
                              invariance["raw_scalar_changes"] and invariance["raw_theta_changes"]
                              and invariance["raw_theta_photon_changes"]))
        for key, residual in invariance["residuals"].items():
            if key == "theta_photon_two" and not charge_conserving:
                continue
            results.append(report(f"C1 ({label}) {key} invariant", is_zero(residual)))
        if charge_conserving:
            results.append(report("C1 ([g, Q] = 0) block-determinant anomaly = "
                                  "N_c tr[Q^2 (G_R - G_L)] through x^2",
                                  all(r == 0 for r in invariance["photon_trace_form_residual"])))

    constant = check_constant_connections()
    results.append(report("C2 constant connections: exp(-ixC_R) M exp(ixC_L) gives the closed "
                          "hat M_2 and the chosen hat S_2, hat P_2",
                          all(is_zero(value) for value in constant.values())))

    numerical = check_numerical_transport()
    print(f"C3 relative deviation of fitted 2nd derivative: "
          f"{mpmath.nstr(numerical['relative_deviation'], 5)} "
          f"(tolerance {mpmath.nstr(numerical['tolerance'], 3)})")
    results.append(report("C3 mpmath transport integration reproduces hat M_2",
                          numerical["passed"]))

    for name, passed in check_groupings_on_matrices().items():
        results.append(report(f"C4 grouping '{name}' on exact random 3x3 matrices", passed))

    regression = check_first_order_regression()
    results.append(report("C5 hat S_1 = -i[c^VS, m_q], hat P_1 = -{c^AP, m_q}",
                          regression["scalar_one_residual"] == 0
                          and regression["pseudoscalar_one_residual"] == 0))
    results.append(report("C5 hat Theta_1 = c_G + tr c^A, hat Theta_gamma,1 = c_gamma + "
                          "2N_c tr(Q^2 c^A)",
                          regression["theta_one_residual"] == 0
                          and regression["theta_photon_one_residual"] == 0))

    vacuum = check_vacuum_angle()
    results.append(report("C6 hat Theta + arg det hat M / 2 = Theta + arg det M / 2 "
                          "at x^1 and x^2",
                          all(r == 0 for r in vacuum["residuals"])
                          and vacuum["second_order_phase_changes"]))

    hermiticity = check_hermiticity()
    results.append(report("C7 chosen hat S_2, hat P_2 Hermitian (symbolic and exact matrices)",
                          hermiticity["symbolic_scalar"] == 0
                          and hermiticity["symbolic_pseudoscalar"] == 0
                          and hermiticity["numeric_scalar"]
                          and hermiticity["numeric_pseudoscalar"]))

    chiral_traces = check_chiral_trace_forms()
    results.append(report("C8 Theta_2 + tr(C_R1 - C_L1)/2 equals D3 with C_A1 = (C_R1 - C_L1)/2",
                          chiral_traces["theta_residual"] == 0))
    results.append(report("C8 Theta_gamma,2 + N_c tr[Q^2 (C_R1 - C_L1)] equals D3",
                          chiral_traces["theta_photon_residual"] == 0))
    results.append(report("C8 D3 hat Theta_gamma,2 = Theta_gamma,2 + (8/3)(C_A1)_uu + "
                          "(2/3)[(C_A1)_dd + (C_A1)_ss] (N_c = 3, generic C_A1)",
                          chiral_traces["theta_photon_components_residual"] == 0))

    real_mass = derive_real_mass_specialization()
    show("D5 hat S_2 (S_0 = m_q, P_0 = 0)", real_mass["display"][0])
    show("D5 hat P_2 (S_0 = m_q, P_0 = 0)", real_mass["display"][1])
    show("D5 hat S_2 with C_R, C_L", real_mass["display_chiral"][0])
    show("D5 hat P_2 with C_R, C_L", real_mass["display_chiral"][1])
    for label, key, transport_index in (("hat S_2", "scalar_residual", 0),
                                        ("hat P_2", "pseudoscalar_residual", 1),
                                        ("hat S_2 (C_R, C_L)", "chiral_scalar_residual", 0),
                                        ("hat P_2 (C_R, C_L)", "chiral_pseudoscalar_residual",
                                         1)):
        passed = report(f"D5 real-mass {label} equals the transport with S_0 = m_q, P_0 = 0",
                        real_mass[key] == 0)
        if not passed:
            show("  residual (candidate - transport)", real_mass[key])
            show("  correct form (transport)", real_mass["transport"][transport_index])
        results.append(passed)
    for name, residual in real_mass["cross_residuals"].items():
        show(f"D5 cross term '{name}'", real_mass["display_cross"][name])
        results.append(report(f"D5 cross term '{name}' equals [C_V0,{{C_A0,m}}] + "
                              f"{{C_A0,[C_V0,m]}} (any m_q)", residual == 0))
    print(f"D5 cross term: {real_mass['cross_monomial_number']} monomials; a single nested "
          f"bracket has at most {real_mass['largest_single_bracket_monomial_number']}")
    results.append(report("D5 no number times a single nested bracket equals the cross term",
                          real_mass["cross_monomial_number"]
                          > real_mass["largest_single_bracket_monomial_number"]))
    real_mass_matrices = check_real_mass_on_matrices()
    results.append(report("D5 matrix data: [C_V0, C_A0] != 0 and [C_V0, m_q] != 0",
                          real_mass_matrices["noncommuting"]))
    for key in ("scalar_residual", "pseudoscalar_residual", "chiral_scalar_residual",
                "chiral_pseudoscalar_residual"):
        passed = report(f"D5 {key} vanishes on exact 3x3 matrices (m_q real positive diagonal)",
                        is_zero(real_mass_matrices[key]))
        if not passed:
            print(f"  residual: {real_mass_matrices[key]}")
        results.append(passed)
    for name, residual in real_mass_matrices["cross_residuals"].items():
        results.append(report(f"D5 cross term '{name}' on exact 3x3 matrices", is_zero(residual)))

    print(f"\n{sum(results)} PASS, {len(results) - sum(results)} FAIL")
    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
