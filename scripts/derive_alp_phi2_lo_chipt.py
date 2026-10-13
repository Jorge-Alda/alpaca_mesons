"""Derive the one-ALP, two-pseudoscalar c_q^V and c_q^S terms at LO.

Invocation:
    .venv/bin/python scripts/derive_alp_phi2_lo_chipt.py

Conventions follow Appendix A of Alda:2026nrd. Natural units are used: fields,
f_a, f_pi, B_0, and m_q have their usual mass dimensions.  The bookkeeping
symbols count ALP and Phi powers and are removed by coefficient extraction.
No WZW term or field heavier than the pseudoscalar nonet is included.
"""

from itertools import combinations_with_replacement

import sympy as sp
try:
    from scripts.symbolic_utils import anticommutator, commutator, map_power, matrix_unit, scalar_coefficient, show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import anticommutator, commutator, map_power, matrix_unit, scalar_coefficient, show


# ---------------------------------------------------------------------------
# Matrix-level derivation from the Appendix-A CCWZ series
# ---------------------------------------------------------------------------
n = sp.Symbol("N_f", integer=True, positive=True)
f_pi, f_a, B_0 = sp.symbols("f_pi f_a B_0", positive=True)
rho, epsilon = sp.symbols("rho epsilon")
a, da = sp.symbols("a partial_mu_a")

Phi = sp.MatrixSymbol("Phi", n, n)
dPhi = sp.MatrixSymbol("partial^mu_Phi", n, n)
cV = sp.MatrixSymbol("c_q^V", n, n)
cS = sp.MatrixSymbol("c_q^S", n, n)
mq = sp.MatrixSymbol("m_q", n, n)

phi = rho * Phi / (2 * f_pi)
dphi = rho * dPhi / (2 * f_pi)

# The O(Phi) pieces are the only terms needed for the kinetic a Phi^2 term.
u_meson = -2 * dphi
u_alp_vector = -2 * sp.I * epsilon * da / f_a * commutator(cV, phi)
kinetic_matrix = sp.expand(
    f_pi**2
    / 4
    * (u_meson * u_alp_vector + u_alp_vector * u_meson)
)
kinetic_raw = sp.Trace(
    scalar_coefficient(scalar_coefficient(kinetic_matrix, epsilon, 1), rho, 2)
)
kinetic_reduced = (
    sp.I * da / (2 * f_a) * sp.Trace(dPhi * commutator(cV, Phi))
)

# The ALP part of chi_S is already O(a).  At O(Phi^2),
# 2 cos(A_phi)[chi_S] contributes -A_phi^2[chi_S].
chi_scalar_alp = 2 * sp.I * B_0 * epsilon * a / f_a * commutator(cS, mq)
chi_plus_scalar = 2 * (
    chi_scalar_alp
    - map_power(anticommutator, chi_scalar_alp, phi, 2) / sp.factorial(2)
)
breaking_matrix = sp.expand(f_pi**2 * chi_plus_scalar / 4)
breaking_raw = sp.Trace(
    scalar_coefficient(scalar_coefficient(breaking_matrix, epsilon, 1), rho, 2)
)
breaking_reduced = (
    -sp.I * B_0 * a / (2 * f_a) * sp.Trace(Phi**2 * commutator(cS, mq))
)

show("Dimensionless pseudoscalar field phi", Phi / (2 * f_pi))
show("Mesonic vielbein needed at O(Phi)", u_meson.subs(rho, 1))
show("Vector-source one-ALP vielbein needed at O(Phi)", u_alp_vector.subs({rho: 1, epsilon: 1}))
show("Raw kinetic coefficient at one ALP and Phi^2", kinetic_raw)
show("Kinetic result after trace cyclicity", kinetic_reduced)
show("ALP-dependent scalar source chi_S", chi_scalar_alp.subs(epsilon, 1))
show("Raw symmetry-breaking coefficient at one ALP and Phi^2", breaking_raw)
show("Symmetry-breaking result after trace cyclicity", breaking_reduced)
show("Total c_q^V plus c_q^S result", kinetic_reduced + breaking_reduced)


# ---------------------------------------------------------------------------
# Generic two-field coefficients
# ---------------------------------------------------------------------------
phi_i, phi_j = sp.symbols("phi_i phi_j")
dphi_i, dphi_j = sp.symbols("partial^mu_phi_i partial^mu_phi_j")
lambda_i = sp.MatrixSymbol("lambda_i", n, n)
lambda_j = sp.MatrixSymbol("lambda_j", n, n)

F_V_ij = sp.Trace(cV * commutator(lambda_i, lambda_j))
F_S_ij = sp.Trace(anticommutator(lambda_i, lambda_j) * commutator(cS, mq))
component_vector = (
    sp.I
    * da
    / (2 * f_a)
    * F_V_ij
    * (phi_i * dphi_j - phi_j * dphi_i)
)
component_scalar = -sp.I * B_0 * a / (2 * f_a) * F_S_ij * phi_i * phi_j

show("Distinct-field c_q^V component interaction", component_vector)
show("Distinct-field c_q^S component interaction", component_scalar)


# ---------------------------------------------------------------------------
# Physical light-channel classification in the repository generator convention
# ---------------------------------------------------------------------------
def E(row: int, column: int) -> sp.ImmutableMatrix:
    return matrix_unit(row, column)


theta_P = sp.Symbol("theta_P", real=True)
eta8 = sp.diag(1, 1, -2) / sp.sqrt(3)
eta0 = sp.sqrt(sp.Rational(2, 3)) * sp.eye(3)
physical_generators = {
    "pi+": sp.sqrt(2) * E(0, 1),
    "pi-": sp.sqrt(2) * E(1, 0),
    "pi0": sp.diag(1, -1, 0),
    "K+": sp.sqrt(2) * E(0, 2),
    "K-": sp.sqrt(2) * E(2, 0),
    "K0": sp.sqrt(2) * E(1, 2),
    "K0bar": sp.sqrt(2) * E(2, 1),
    "eta": sp.cos(theta_P) * eta8 - sp.sin(theta_P) * eta0,
    "eta_prime": sp.sin(theta_P) * eta8 + sp.cos(theta_P) * eta0,
}
charges = {
    "pi+": 1,
    "pi-": -1,
    "pi0": 0,
    "K+": 1,
    "K-": -1,
    "K0": 0,
    "K0bar": 0,
    "eta": 0,
    "eta_prime": 0,
}

mu, md, ms = sp.symbols("m_u m_d m_s", real=True)
v_u, v_d, v_s, v_sd = sp.symbols("c_u^V c_d^V c_s^V c_sd^V")
s_sd = sp.Symbol("c_sd^S")
mq3 = sp.diag(mu, md, ms)
cV_diagonal = sp.diag(v_u, v_d, v_s)
# E_sd closes a flavour trace with an S=+1 final state; its Hermitian-conjugate
# source gives the charge-conjugate S=-1 channels.
cV_sd = v_sd * E(2, 1)
cS_sd = s_sd * E(2, 1)


def vector_flavour(source, first, second):
    return sp.trigsimp(sp.trace(source * commutator(first, second)))


def scalar_flavour(source, first, second):
    return sp.trigsimp(
        sp.trace(anticommutator(first, second) * commutator(source, mq3))
    )


neutral_pairs = [
    (first_name, second_name)
    for first_name, second_name in combinations_with_replacement(physical_generators, 2)
    if charges[first_name] + charges[second_name] == 0
]


def nonzero_channels(coefficient_function, source):
    result = {}
    for first_name, second_name in neutral_pairs:
        coefficient = sp.simplify(
            coefficient_function(
                source,
                physical_generators[first_name],
                physical_generators[second_name],
            )
        )
        if coefficient != 0:
            result[(first_name, second_name)] = coefficient
    return result


vector_diagonal_channels = nonzero_channels(vector_flavour, cV_diagonal)
vector_sd_channels = nonzero_channels(vector_flavour, cV_sd)
scalar_diagonal_channels = nonzero_channels(scalar_flavour, sp.diag(*sp.symbols("c_u^S c_d^S c_s^S")))
scalar_sd_channels = nonzero_channels(scalar_flavour, cS_sd)

print("Physical Q=0 channels with nonzero flavour coefficients:")
for heading, channels in (
    ("diagonal c_q^V", vector_diagonal_channels),
    ("E_sd c_q^V (S=+1; conjugate source gives S=-1)", vector_sd_channels),
    ("diagonal c_q^S", scalar_diagonal_channels),
    ("E_sd c_q^S (S=+1; conjugate source gives S=-1)", scalar_sd_channels),
):
    print(f"  {heading}:")
    if not channels:
        print("    none")
    for channel, coefficient in channels.items():
        print(f"    {channel[0]} {channel[1]}: {sp.latex(coefficient)}")
print()

# On shell, the vector-source momentum contraction reduces to a mass-squared
# difference.  This makes the diagonal-source particle-antiparticle channels
# vanish in the exact charge-conjugate mass limit; their off-shell vertex is
# nevertheless nonzero.
m_a, m_i, m_j = sp.symbols("m_a m_i m_j", positive=True)
p_a_dot_p_i = (m_a**2 + m_i**2 - m_j**2) / 2
p_a_dot_p_j = (m_a**2 + m_j**2 - m_i**2) / 2
on_shell_vector_contraction = sp.simplify(p_a_dot_p_j - p_a_dot_p_i)
show("On-shell p_a dot (p_j-p_i)", on_shell_vector_contraction)


# ---------------------------------------------------------------------------
# Deterministic verification checks
# ---------------------------------------------------------------------------
Phi_test = sp.Matrix([[1, 2], [3, 5]])
dPhi_test = sp.Matrix([[7, 11], [13, 17]])
cV_test = sp.Matrix([[19, 23], [29, 31]])
cS_test = sp.Matrix([[37, 41], [43, 47]])
mq_test = sp.diag(53, 59)

assert sp.simplify(
    kinetic_raw.doit().subs({Phi: Phi_test, dPhi: dPhi_test, cV: cV_test})
    - kinetic_reduced.doit().subs({Phi: Phi_test, dPhi: dPhi_test, cV: cV_test})
) == 0
assert sp.simplify(
    breaking_raw.doit().subs({Phi: Phi_test, cS: cS_test, mq: mq_test})
    - breaking_reduced.doit().subs({Phi: Phi_test, cS: cS_test, mq: mq_test})
) == 0
assert scalar_diagonal_channels == {}
assert set(vector_diagonal_channels) == {
    ("pi+", "pi-"),
    ("K+", "K-"),
    ("K0", "K0bar"),
}
assert set(vector_sd_channels) == {
    ("pi-", "K+"),
    ("pi0", "K0"),
    ("K0", "eta"),
    ("K0", "eta_prime"),
}
assert set(scalar_sd_channels) == set(vector_sd_channels)
assert sp.simplify(on_shell_vector_contraction - (m_j**2 - m_i**2)) == 0

print("Verification: all symbolic trace, channel-classification, and on-shell checks passed.")
