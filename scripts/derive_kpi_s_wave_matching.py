"""Derive the Alda-consistent a -> K+ pi0 K- scalar production kernel.

Invocation:
    .venv/bin/python scripts/derive_kpi_s_wave_matching.py

The calculation uses natural units, the mostly-minus metric, the isospin
limit, diagonal hatted ALP-quark couplings, and the meson-generator
normalization of Alda:2026nrd. The amplitude convention is S = 1 + i T.
"""

from itertools import permutations

import sympy as sp
try:
    from scripts.symbolic_utils import commutator, show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import commutator, show


def kinetic_tensor(
    first: sp.MatrixBase, second: sp.MatrixBase, differentiated: sp.MatrixBase
) -> sp.Expr:
    matrix = commutator(first, commutator(second, differentiated)) + commutator(
        second, commutator(first, differentiated)
    )
    return sp.expand(sp.trace(cA * matrix))


f_a, f_pi, B_0 = sp.symbols("f_a f_pi B_0", positive=True)
m_a, m_K, m_pi, m_hat, m_s = sp.symbols(
    "m_a m_K m_pi m_hat m_s", positive=True
)
s, t, u, z = sp.symbols("s t u z", real=True)
A_u, A_d, A_s = sp.symbols("A_u A_d A_s", real=True)
P_u, P_d, P_s = sp.symbols("P_u P_d P_s", real=True)

sqrt2 = sp.sqrt(2)
K_plus = sp.zeros(3)
K_plus[0, 2] = sqrt2
K_minus = sp.zeros(3)
K_minus[2, 0] = sqrt2
pi_zero = sp.diag(1, -1, 0)
cA = sp.diag(A_u, A_d, A_s)
cP = sp.diag(P_u, P_d, P_s)
mq = sp.diag(m_hat, m_hat, m_s)
P_matrix = cP * mq + mq * cP

# The field ordering is K+(p1), pi0(p2), K-(p3). The semicolon identifies
# the differentiated field in the generic three-distinct-meson tensors.
T_Kp_pi_Km = kinetic_tensor(K_plus, pi_zero, K_minus)
T_Kp_Km_pi = kinetic_tensor(K_plus, K_minus, pi_zero)
T_pi_Km_Kp = kinetic_tensor(pi_zero, K_minus, K_plus)

generators = (K_plus, pi_zero, K_minus)
ordered_sum = sp.zeros(3)
for ordering in permutations(generators):
    ordered_sum += ordering[0] * ordering[1] * ordering[2]
T_P = sp.expand(sp.trace(ordered_sum * P_matrix))

show("lambda_K+", K_plus)
show("lambda_pi0", pi_zero)
show("lambda_K-", K_minus)
show("T^A_(K+ pi0; K-)", T_Kp_pi_Km)
show("T^A_(K+ K-; pi0)", T_Kp_Km_pi)
show("T^A_(pi0 K-; K+)", T_pi_Km_Kp)
show("Jacobi sum", sp.simplify(T_Kp_pi_Km + T_Kp_Km_pi + T_pi_Km_Kp))
show("Sum of the six ordered generator products", ordered_sum)
show("T^P_(K+ pi0 K-)", T_P)

K_coefficient = 1 / (6 * f_a * f_pi)
M_coefficient = -B_0 / (6 * f_a * f_pi)

# Generic three-distinct amplitude, with s=(p1+p2)^2,
# t=(p2+p3)^2, and u=(p1+p3)^2.
amplitude_raw = (
    K_coefficient
    / 2
    * (
        T_Kp_pi_Km * (m_a**2 + m_K**2 - s)
        + T_Kp_Km_pi * (m_a**2 + m_pi**2 - u)
        + T_pi_Km_Kp * (m_a**2 + m_K**2 - t)
    )
    + M_coefficient * T_P
)

Delta_A = sp.Symbol("Delta_A", real=True)
amplitude_compact = sp.factor(amplitude_raw.subs(A_u, Delta_A + A_s))
amplitude_compact = sp.collect(amplitude_compact, (Delta_A, B_0))
expected_amplitude = (
    Delta_A * (s + t - 2 * u + 2 * m_pi**2 - 2 * m_K**2)
    - 4 * B_0 * (2 * P_u * m_hat + P_s * m_s)
) / (6 * f_a * f_pi)

show("Raw LO invariant amplitude", amplitude_raw)
show("LO invariant amplitude with Delta_A = A_u - A_s", amplitude_compact)
show("Compact LO invariant amplitude", expected_amplitude)

mandelstam_u = m_a**2 + 2 * m_K**2 + m_pi**2 - s - t
amplitude_st = sp.factor(expected_amplitude.subs(u, mandelstam_u))
band_kernel = (
    Delta_A * (3 * s - m_a**2 - 3 * m_K**2)
    - 2 * B_0 * (2 * P_u * m_hat + P_s * m_s)
) / (6 * f_a * f_pi)
band_kernel_t = band_kernel.subs(s, t)

show("LO amplitude after eliminating u", amplitude_st)
show("Single K pi band production kernel P_LO(s)", band_kernel)
show("Crossed-band production kernel P_LO(t)", band_kernel_t)
show("P_LO(s) + P_LO(t)", sp.factor(band_kernel + band_kernel_t))

# Project the complete contact amplitude onto the s-channel S wave. This is
# distinct from the single-band kernel because the crossed t-channel term also
# contributes to the s-channel projection.
kallen = lambda x, y, zz: x**2 + y**2 + zz**2 - 2 * x * y - 2 * x * zz - 2 * y * zz
q = sp.sqrt(kallen(s, m_K**2, m_pi**2)) / (2 * sp.sqrt(s))
k = sp.sqrt(kallen(m_a**2, s, m_K**2)) / (2 * sp.sqrt(s))
E1 = (s + m_K**2 - m_pi**2) / (2 * sp.sqrt(s))
E2 = (s + m_pi**2 - m_K**2) / (2 * sp.sqrt(s))
E3 = (m_a**2 - s - m_K**2) / (2 * sp.sqrt(s))
u_sz = 2 * m_K**2 + 2 * E1 * E3 - 2 * q * k * z
t_sz = m_K**2 + m_pi**2 + 2 * E2 * E3 + 2 * q * k * z
amplitude_sz = sp.expand(expected_amplitude.subs({t: t_sz, u: u_sz}))
amplitude_s_wave = sp.factor(sp.integrate(amplitude_sz, (z, -1, 1)) / 2)
amplitude_non_s_wave = sp.factor(amplitude_sz - amplitude_s_wave)

show("Complete amplitude in the K+ pi0 pair frame", amplitude_sz)
show("Complete-amplitude K+ pi0 S-wave projection", amplitude_s_wave)
show("Residual angular term", amplitude_non_s_wave)

# Aloni's per-band momentum polynomial is quadratic in s, whereas the
# chiral production kernel is linear. Their ratio therefore cannot be a
# universal constant.
aloni_momentum_polynomial = (
    (m_a**2 + m_K**2 - s) * (s - m_K**2 - m_pi**2) / 4
)
forced_prefactor = sp.factor(band_kernel / aloni_momentum_polynomial)
show("Aloni per-band momentum polynomial", aloni_momentum_polynomial)
show("Forced coefficient before an S_Kpi normalization is specified", forced_prefactor)

# Diagnostics.
assert sp.simplify(T_Kp_pi_Km - 2 * (A_s - A_u)) == 0
assert sp.simplify(T_Kp_Km_pi - 4 * (A_u - A_s)) == 0
assert sp.simplify(T_pi_Km_Kp - 2 * (A_s - A_u)) == 0
assert sp.simplify(T_Kp_pi_Km + T_Kp_Km_pi + T_pi_Km_Kp) == 0
assert ordered_sum == sp.diag(4, 0, 2)
assert sp.simplify(T_P - (8 * P_u * m_hat + 4 * P_s * m_s)) == 0
assert sp.simplify(amplitude_compact - expected_amplitude) == 0
assert sp.simplify(amplitude_st - band_kernel - band_kernel_t) == 0
assert sp.simplify(amplitude_non_s_wave.subs(z, 0)) == 0
assert sp.simplify(sp.integrate(amplitude_non_s_wave, (z, -1, 1))) == 0
