"""Scalar and pseudoscalar quark-mass sources of an ALP from a Peccei-Quinn phase.

Run (from the repository root):
    .venv/bin/python scripts/derive_alp_pq_phase_mass_sources.py
    .venv/bin/python scripts/derive_alp_pq_phase_mass_sources.py --max-order 8

Conventions and units
---------------------
Natural units. Quark flavours q = (u, d, s). The quark mass matrix
m_q = diag(m_u, m_d, m_s) is real, positive and has mass dimension one; the
Peccei-Quinn (PQ) charge matrix N = diag(n_u, n_d, n_s) is real and
dimensionless (integer in the models used, kept symbolic here). The ALP
variable x = a/f_a is real and dimensionless.

Memory convention (memory/alp/shared/alp-general-chiral-redefinition-invariants.md):
the partonic Lagrangian contains
    -bar(q_R) X q_L - bar(q_L) X^dagger q_R = -bar(q) (S - i P gamma_5) q,
with X(x) = S(x) + i P(x), S and P Hermitian, X(0) = m_q. Right- and
left-handed derivative-current coefficients C_R, C_L and the topological
coefficient Theta are zero / ALP-independent in the basis used here, where
the whole phase is carried by the mass term.

Taylor-coefficient convention: F(x) = sum_k x^k / k! F_k, i.e.
F_k = d^k F / dx^k at x = 0.

Inputs
------
* memory/extracts/DiLuzio:2026poh/model-lagrangian-and-field-normalization.md
  and angular-mode-couplings.md: the mass term of DiLuzio:2026poh,
  Eqs. (3), (7), -m_q (sqrt2 phi)^{n_q} bar(q_L) q_R + h.c., with
  phi = (1 + sigma)/sqrt2 exp(i a) (Eq. (9)). At sigma = 0 this is
  -m_q exp(i n_q a) bar(q_L) q_R + h.c.
* memory/extracts/DiLuzio:2026poh/angular-mode-parametrization.md: a is the
  dimensionless phase, a = a_can / f_phi.
* memory/alp/shared/alp-chiral-source-transformations.md: source
  transformation laws M' = g_R M g_L^dagger,
  C_R' = g_R C_R g_R^dagger + i g_R d_x g_R^dagger, transport
  d_x W = i C W with W(0) = 1, dressed source hat(M) = W_R^dagger M W_L.
* memory/extracts/Alda:2026nrd/appendix-a-expansion.md (through
  memory/alp/shared/alp-dressed-mass-first-order-invariants.md): first-order
  parametrisation S_1 = i [c^S, m_q], P_1 = {c^P, m_q}.

Assumptions (stated, not derived)
---------------------------------
* The radial mode is set to zero, sigma = 0.
* x is identified with the paper's dimensionless phase a, i.e.
  x = a_can / f_phi. This equals a_can / f_a only if f_a is taken equal to
  f_phi. With the paper's light-axion definition f_a = f_phi / sum_q n_q
  (Eq. (53)), every N below is to be read as N / sum_q n_q.

Results
-------
X(x) = m_q exp(-i N x), so S(x) = m_q cos(N x), P(x) = -m_q sin(N x), and
S_{2j} = (-1)^j m_q N^{2j},   S_{2j+1} = 0,
P_{2j} = 0,                   P_{2j+1} = -(-1)^j m_q N^{2j+1}.
First order: c^S = 0 (any diagonal c^S is equivalent), c^P = -N/2.

All checks raise AssertionError on failure; the script exits with status 1.
"""

import argparse
import sys

import sympy as sp

ALP_FIELD_RATIO = sp.Symbol("x", real=True)
GAMMA_FIVE = sp.Symbol("gamma_5")
FLAVOURS = ("u", "d", "s")

QUARK_MASSES = sp.symbols("m_u m_d m_s", positive=True)
PQ_CHARGES = sp.symbols("n_u n_d n_s", real=True)
REDEFINITION_PHASES = sp.symbols("c_u c_d c_s", real=True)

DEFAULT_MAX_ORDER = 6


def quark_mass_matrix(masses=QUARK_MASSES):
    """Return m_q = diag(m_u, m_d, m_s)."""
    return sp.diag(*masses)


def pq_charge_matrix(charges=PQ_CHARGES):
    """Return N = diag(n_u, n_d, n_s)."""
    return sp.diag(*charges)


def mass_source(masses=QUARK_MASSES, charges=PQ_CHARGES):
    """Return X(x) = m_q exp(-i N x) for diagonal m_q and N."""
    return sp.diag(*[mass * sp.exp(-sp.I * charge * ALP_FIELD_RATIO)
                     for mass, charge in zip(masses, charges)])


def hermitian_parts(matrix):
    """Split X = S + i P into Hermitian S and P (x and entries real)."""
    scalar = (matrix + matrix.H) / 2
    pseudoscalar = (matrix - matrix.H) / (2 * sp.I)
    return scalar.applyfunc(sp.simplify), pseudoscalar.applyfunc(sp.simplify)


def taylor_coefficients(matrix, max_order):
    """Return [F_0, ..., F_K] with F_k = d^k F/dx^k at x = 0."""
    coefficients = []
    derivative = matrix
    for _ in range(max_order + 1):
        coefficients.append(
            derivative.subs(ALP_FIELD_RATIO, 0).applyfunc(sp.simplify))
        derivative = derivative.diff(ALP_FIELD_RATIO)
    return coefficients


def source_taylor_coefficients(max_order=DEFAULT_MAX_ORDER,
                               masses=QUARK_MASSES, charges=PQ_CHARGES):
    """Return the lists S_k and P_k, k = 0..max_order."""
    scalar, pseudoscalar = hermitian_parts(mass_source(masses, charges))
    return (taylor_coefficients(scalar, max_order),
            taylor_coefficients(pseudoscalar, max_order))


def closed_form_scalar(order, masses=QUARK_MASSES, charges=PQ_CHARGES):
    """S_k = (-1)^j m_q N^{2j} for k = 2j, zero for odd k."""
    if order % 2:
        return sp.zeros(3, 3)
    half_order = order // 2
    return (-1) ** half_order * quark_mass_matrix(masses) \
        * pq_charge_matrix(charges) ** order


def closed_form_pseudoscalar(order, masses=QUARK_MASSES, charges=PQ_CHARGES):
    """P_k = -(-1)^j m_q N^{2j+1} for k = 2j + 1, zero for even k."""
    if order % 2 == 0:
        return sp.zeros(3, 3)
    half_order = (order - 1) // 2
    return -(-1) ** half_order * quark_mass_matrix(masses) \
        * pq_charge_matrix(charges) ** order


def check_closed_forms(max_order=DEFAULT_MAX_ORDER, masses=QUARK_MASSES,
                       charges=PQ_CHARGES, scalar_coefficients=None,
                       pseudoscalar_coefficients=None):
    """Assert S_k, P_k equal the closed forms for k = 0..max_order."""
    if scalar_coefficients is None:
        scalar_coefficients, pseudoscalar_coefficients = \
            source_taylor_coefficients(max_order, masses, charges)
    for order in range(max_order + 1):
        scalar_residual = (scalar_coefficients[order]
                           - closed_form_scalar(order, masses, charges)
                           ).applyfunc(sp.simplify)
        pseudoscalar_residual = (pseudoscalar_coefficients[order]
                                 - closed_form_pseudoscalar(order, masses,
                                                            charges)
                                 ).applyfunc(sp.simplify)
        assert scalar_residual == sp.zeros(3, 3), f"S_{order} closed form"
        assert pseudoscalar_residual == sp.zeros(3, 3), \
            f"P_{order} closed form"
    return True


def dirac_coefficients_from_chiral(right_left_matrix, left_right_matrix):
    """Coefficients of bar(q) q and bar(q) gamma_5 q of a chiral mass term.

    The input is the Lagrangian
        -bar(q_R) A q_L - bar(q_L) B q_R,
    with A = right_left_matrix, B = left_right_matrix (flavour diagonal).
    Using bar(q_R) q_L = bar(q) (1 - gamma_5)/2 q and
    bar(q_L) q_R = bar(q) (1 + gamma_5)/2 q, the Lagrangian is
    bar(q) (D_1 + D_5 gamma_5) q; the function returns (D_1, D_5).
    gamma_5 is represented by a commuting symbol, which is exact for these
    flavour-diagonal bilinears.
    """
    projector_left = (1 - GAMMA_FIVE) / 2
    projector_right = (1 + GAMMA_FIVE) / 2
    dirac_matrix = (-right_left_matrix * projector_left
                    - left_right_matrix * projector_right)
    dirac_matrix = dirac_matrix.applyfunc(sp.expand)
    unit_part = dirac_matrix.applyfunc(lambda entry: entry.coeff(GAMMA_FIVE, 0))
    gamma_five_part = dirac_matrix.applyfunc(
        lambda entry: entry.coeff(GAMMA_FIVE, 1))
    return unit_part, gamma_five_part


def check_sign_mapping(masses=QUARK_MASSES, charges=PQ_CHARGES):
    """Assert that the paper's mass term equals -bar(q)(S - i P gamma_5) q.

    Paper (DiLuzio:2026poh Eqs. (3), (7), sigma = 0):
        -m exp(i N x) bar(q_L) q_R + h.c.
      = -m [cos(N x) bar(q) q + i sin(N x) bar(q) gamma_5 q].
    Memory: -bar(q) (S - i P gamma_5) q, with X = m exp(-i N x).
    Returns the four coefficient matrices.
    """
    mass_matrix = quark_mass_matrix(masses)
    phase = sp.diag(*[sp.exp(sp.I * charge * ALP_FIELD_RATIO)
                      for charge in charges])
    # Paper: bar(q_L) q_R carries m e^{iNx}; its h.c. bar(q_R) q_L carries
    # m e^{-iNx}.
    paper_unit, paper_gamma_five = dirac_coefficients_from_chiral(
        mass_matrix * phase.H, mass_matrix * phase)
    paper_unit_expected = -mass_matrix * sp.diag(
        *[sp.cos(charge * ALP_FIELD_RATIO) for charge in charges])
    paper_gamma_five_expected = -sp.I * mass_matrix * sp.diag(
        *[sp.sin(charge * ALP_FIELD_RATIO) for charge in charges])
    for computed, expected, label in (
            (paper_unit, paper_unit_expected, "paper bar(q) q"),
            (paper_gamma_five, paper_gamma_five_expected,
             "paper bar(q) gamma_5 q")):
        residual = (computed - expected).applyfunc(
            lambda entry: sp.simplify(entry.rewrite(sp.cos)))
        assert residual == sp.zeros(3, 3), label

    source = mass_source(masses, charges)
    scalar, pseudoscalar = hermitian_parts(source)
    memory_unit, memory_gamma_five = dirac_coefficients_from_chiral(
        source, source.H)
    # -bar(q)(S - i P gamma_5) q has coefficients (-S, +i P).
    for computed, expected, label in (
            (memory_unit, -scalar, "memory bar(q) q = -S"),
            (memory_gamma_five, sp.I * pseudoscalar,
             "memory bar(q) gamma_5 q = i P")):
        residual = (computed - expected).applyfunc(
            lambda entry: sp.simplify(entry.rewrite(sp.cos)))
        assert residual == sp.zeros(3, 3), label
    for computed, expected, label in (
            (memory_unit, paper_unit_expected, "scalar paper = memory"),
            (memory_gamma_five, paper_gamma_five_expected,
             "pseudoscalar paper = memory")):
        residual = (computed - expected).applyfunc(
            lambda entry: sp.simplify(entry.rewrite(sp.cos)))
        assert residual == sp.zeros(3, 3), label
    return {
        "paper_unit": paper_unit_expected,
        "paper_gamma_five": paper_gamma_five_expected,
        "scalar": scalar,
        "pseudoscalar": pseudoscalar,
    }


def first_order_couplings(charges=PQ_CHARGES):
    """Return (c^S, c^P) = (0, -N/2) for the first-order parametrisation."""
    return sp.zeros(3, 3), -pq_charge_matrix(charges) / 2


def check_first_order_matching(scalar_first, pseudoscalar_first,
                               masses=QUARK_MASSES, charges=PQ_CHARGES):
    """Assert S_1 = i[c^S, m_q] and P_1 = {c^P, m_q} with c^S=0, c^P=-N/2.

    Also asserts the explicit values S_1 = 0, P_1 = -N m_q.
    """
    mass_matrix = quark_mass_matrix(masses)
    scalar_coupling, pseudoscalar_coupling = first_order_couplings(charges)
    scalar_parametrised = sp.I * (scalar_coupling * mass_matrix
                                  - mass_matrix * scalar_coupling)
    pseudoscalar_parametrised = (pseudoscalar_coupling * mass_matrix
                                 + mass_matrix * pseudoscalar_coupling)
    zero = sp.zeros(3, 3)
    assert (scalar_first - scalar_parametrised).applyfunc(sp.simplify) \
        == zero, "S_1 = i[c^S, m_q]"
    assert (pseudoscalar_first - pseudoscalar_parametrised).applyfunc(
        sp.simplify) == zero, "P_1 = {c^P, m_q}"
    assert scalar_first.applyfunc(sp.simplify) == zero, "S_1 = 0"
    assert (pseudoscalar_first + pq_charge_matrix(charges) * mass_matrix
            ).applyfunc(sp.simplify) == zero, "P_1 = -N m_q"
    return pseudoscalar_coupling


def leibniz_product(left_coefficients, right_coefficients, order):
    """Order-k Taylor coefficient of a product F G (order preserved)."""
    return sum((sp.binomial(order, index) * left_coefficients[index]
                * right_coefficients[order - index]
                for index in range(order + 1)), sp.zeros(3, 3))


def transport_coefficients(connection_coefficients, max_order):
    """Taylor coefficients of W solving d_x W = i C W, W(0) = 1.

    Recursion: W_{k+1} = i sum_j binom(k, j) C_j W_{k-j}.
    """
    transport = [sp.eye(3)]
    for order in range(max_order):
        transport.append((sp.I * leibniz_product(
            connection_coefficients, transport, order)).applyfunc(sp.expand))
    return transport


def check_redefinition_invariance(max_order=DEFAULT_MAX_ORDER,
                                  redefinition_generator=None,
                                  masses=QUARK_MASSES, charges=PQ_CHARGES):
    """Assert the dressed source is unchanged by g_R = exp(i c x), g_L = 1.

    Steps, all at the level of Taylor coefficients up to max_order:
    1. g_R,k = (i c)^k, (g_R^dagger)_k = (-i c)^k.
    2. M'_k from M' = g_R M (Leibniz rule; g_L = 1).
    3. C_R'_k from C_R' = g_R C_R g_R^dagger + i g_R d_x g_R^dagger with
       C_R = 0, using (d_x g_R^dagger)_k = (g_R^dagger)_{k+1}.
    4. W_R'_k from the transport recursion; W_L' = 1.
    5. hat(M)_k = sum binom W_R'^dagger_j M'_{k-j}, compared with M_k.
    The default generator c = diag(c_u, c_d, c_s); any Hermitian matrix,
    including a non-diagonal one, may be passed.
    Returns the list of C_R'_k.
    """
    if redefinition_generator is None:
        redefinition_generator = sp.diag(*REDEFINITION_PHASES)
    mass_coefficients = taylor_coefficients(mass_source(masses, charges),
                                            max_order + 1)
    rotation = [(sp.I * redefinition_generator) ** order
                for order in range(max_order + 2)]
    rotation_conjugate = [(-sp.I * redefinition_generator) ** order
                          for order in range(max_order + 2)]
    rotation_conjugate_derivative = rotation_conjugate[1:]
    transformed_mass = [leibniz_product(rotation, mass_coefficients, order)
                        for order in range(max_order + 1)]
    transformed_connection = [
        (sp.I * leibniz_product(rotation, rotation_conjugate_derivative,
                                order)).applyfunc(sp.expand)
        for order in range(max_order + 1)]
    assert transformed_connection[0] == redefinition_generator.applyfunc(
        sp.expand), "C_R'_0 = c"
    for order in range(1, max_order + 1):
        assert transformed_connection[order] == sp.zeros(3, 3), \
            f"C_R'_{order} = 0"
    transport_right = transport_coefficients(transformed_connection,
                                             max_order)
    transport_right_conjugate = [matrix.H for matrix in transport_right]
    for order in range(max_order + 1):
        dressed = leibniz_product(transport_right_conjugate, transformed_mass,
                                  order)
        residual = (dressed - mass_coefficients[order]).applyfunc(sp.simplify)
        assert residual == sp.zeros(3, 3), f"hat(M)_{order} = M_{order}"
    return transformed_connection


def check_numerical_examples(max_order=DEFAULT_MAX_ORDER):
    """Closed forms with exact rational masses and n = (1,0,0), (2,-1,1)."""
    masses = (sp.Rational(11, 5), sp.Rational(47, 10), sp.Integer(93))
    for charges in ((1, 0, 0), (2, -1, 1)):
        charges = tuple(sp.Integer(value) for value in charges)
        check_closed_forms(max_order, masses, charges)
        check_sign_mapping(masses, charges)
        scalar_coefficients, pseudoscalar_coefficients = \
            source_taylor_coefficients(1, masses, charges)
        check_first_order_matching(scalar_coefficients[1],
                                   pseudoscalar_coefficients[1],
                                   masses, charges)
    return True


def print_summary(scalar_coefficients, pseudoscalar_coefficients,
                  use_latex=False):
    """Print diagonal entries (u, d, s) of S_k and P_k."""
    formatter = sp.latex if use_latex else sp.sstr
    print("Taylor coefficients F_k = d^k F/dx^k at x=0; diagonal (u, d, s)")
    for order, (scalar, pseudoscalar) in enumerate(
            zip(scalar_coefficients, pseudoscalar_coefficients)):
        scalar_entries = [formatter(sp.factor(scalar[index, index]))
                          for index in range(3)]
        pseudoscalar_entries = [formatter(sp.factor(pseudoscalar[index, index]))
                                for index in range(3)]
        print(f"S_{order} = diag({', '.join(scalar_entries)})")
        print(f"P_{order} = diag({', '.join(pseudoscalar_entries)})")


def main(arguments=None):
    """Run all derivations and checks; return the process exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-order", type=int, default=DEFAULT_MAX_ORDER)
    parser.add_argument("--latex", action="store_true")
    options = parser.parse_args(arguments)
    max_order = options.max_order

    scalar_coefficients, pseudoscalar_coefficients = \
        source_taylor_coefficients(max_order)
    print_summary(scalar_coefficients, pseudoscalar_coefficients,
                  options.latex)
    checks = (
        ("closed forms S_k, P_k",
         lambda: check_closed_forms(max_order,
                                    scalar_coefficients=scalar_coefficients,
                                    pseudoscalar_coefficients=
                                    pseudoscalar_coefficients)),
        ("sign mapping paper -> memory", check_sign_mapping),
        ("first-order matching c^S = 0, c^P = -N/2",
         lambda: check_first_order_matching(scalar_coefficients[1],
                                            pseudoscalar_coefficients[1])),
        ("invariance under g_R = exp(i c x), diagonal c",
         lambda: check_redefinition_invariance(max_order)),
        ("invariance under g_R = exp(i c x), non-diagonal rational c",
         lambda: check_redefinition_invariance(
             max_order, sp.Matrix([[1, 2 - sp.I, 0],
                                   [2 + sp.I, sp.Rational(-1, 3), sp.I],
                                   [0, -sp.I, 2]]))),
        ("rational spot checks n=(1,0,0), (2,-1,1)",
         lambda: check_numerical_examples(max_order)),
    )
    failures = 0
    for label, check in checks:
        try:
            check()
            print(f"PASS  {label}")
        except AssertionError as error:
            failures += 1
            print(f"FAIL  {label}: {error}")
    _, pseudoscalar_coupling = first_order_couplings()
    print(f"c^P = {sp.sstr(pseudoscalar_coupling)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
