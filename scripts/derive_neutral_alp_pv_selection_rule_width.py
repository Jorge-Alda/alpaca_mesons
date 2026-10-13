#!/usr/bin/env python3
"""Derive neutral diagonal a->PV flavor zeros and derivative-vertex widths.

Invocation:
    .venv/bin/python scripts/derive_neutral_alp_pv_selection_rule_width.py

Natural units and the mostly-minus metric are used.  The calculation keeps a
generic 3x3 ALP spurion and does not impose its Hermiticity or diagonality.
"""

import sympy as sp
try:
    from scripts.symbolic_utils import show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import show



def minkowski_dot(left: sp.Matrix, right: sp.Matrix, metric: sp.Matrix) -> sp.Expr:
    """Return the mostly-minus scalar product of two contravariant vectors."""
    return sp.expand((left.T * metric * right)[0])


# ---------------------------------------------------------------------------
# 1. Conventions, units, and normalized diagonal U(3) generators
# ---------------------------------------------------------------------------
GeV = sp.Symbol("GeV", positive=True)
f_pi = sp.Rational(93, 1000) * GeV
F = f_pi  # LO convention U=exp(i sqrt(2) Phi/F), as recorded in memory.

theta_P, theta_V = sp.symbols("theta_P theta_V", real=True)
pi0 = sp.diag(1, -1, 0) / sp.sqrt(2)
eta8 = sp.diag(1, 1, -2) / sp.sqrt(6)
eta0 = sp.eye(3) / sp.sqrt(3)

# Pseudoscalar convention: eta=cos(theta_P) eta8-sin(theta_P) eta0.
eta = sp.cos(theta_P) * eta8 - sp.sin(theta_P) * eta0
eta_prime = sp.sin(theta_P) * eta8 + sp.cos(theta_P) * eta0

# Vector convention: omega=cos(theta_V) eta0+sin(theta_V) eta8.  The phi
# phase below is chosen so ideal phi=diag(0,0,1), matching memory/rcht-basics.
omega = sp.cos(theta_V) * eta0 + sp.sin(theta_V) * eta8
phi = sp.sin(theta_V) * eta0 - sp.cos(theta_V) * eta8
omega_ideal = sp.diag(1, 1, 0) / sp.sqrt(2)
phi_ideal = sp.diag(0, 0, 1)
ideal_vector_mixing = {
    sp.cos(theta_V): sp.sqrt(sp.Rational(2, 3)),
    sp.sin(theta_V): 1 / sp.sqrt(3),
}

show("Explicit pion decay constant f_pi [GeV]", f_pi / GeV)
show("Normalized pi0 generator", pi0)
show("Normalized eta8 generator", eta8)
show("Normalized eta0 generator", eta0)
show("Generic eta(theta_P) generator", eta)
show("Generic eta-prime(theta_P) generator", eta_prime)
show("Generic omega(theta_V) generator", omega)
show("Generic phi(theta_V) generator", phi)
show("Ideal omega generator", omega_ideal)
show("Ideal phi generator (chosen positive strange phase)", phi_ideal)

base_generators = (pi0, eta8, eta0)
gram = sp.Matrix([[sp.trace(a * b) for b in base_generators] for a in base_generators])
eta_gram = sp.simplify(
    sp.Matrix([[sp.trace(a * b) for b in (eta, eta_prime)] for a in (eta, eta_prime)])
)
vector_gram = sp.simplify(
    sp.Matrix([[sp.trace(a * b) for b in (omega, phi)] for a in (omega, phi)])
)
show("Gram matrix Tr(T_i T_j) for pi0, eta8, eta0", gram)
show("Generic eta/eta-prime Gram matrix", eta_gram)
show("Generic omega/phi Gram matrix", vector_gram)
show("Omega at ideal vector mixing", sp.simplify(omega.subs(ideal_vector_mixing)))
show("Phi at ideal vector mixing", sp.simplify(phi.subs(ideal_vector_mixing)))


# ---------------------------------------------------------------------------
# 2. Generic diagonal flavor selection rules
# ---------------------------------------------------------------------------
p_u, p_d, p_s = sp.symbols("p_u p_d p_s")
v_u, v_d, v_s = sp.symbols("v_u v_d v_s")
t_u, t_d, t_s = sp.symbols("t_u t_d t_s")
spurion_entries = sp.symbols(
    "c_uu c_ud c_us c_du c_dd c_ds c_su c_sd c_ss"
)
C_A = sp.Matrix(3, 3, spurion_entries)
P_diag = sp.diag(p_u, p_d, p_s)
V_diag = sp.diag(v_u, v_d, v_s)
T_gamma_diag = sp.diag(t_u, t_d, t_s)
Q = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))

commutator_PC = P_diag * C_A - C_A * P_diag
tensor_vector_trace = sp.expand(sp.trace(V_diag * commutator_PC))
mixing_photon_trace = sp.expand(sp.trace(T_gamma_diag * commutator_PC))
direct_trace = sp.expand(sp.trace((Q * C_A - C_A * Q) * P_diag))

show("Generic ALP spurion C_A", C_A)
show("Generic tensor-origin Tr(V [P,C_A])", tensor_vector_trace)
show("Generic mixing-induced Tr(T_gamma [P,C_A])", mixing_photon_trace)
show("Generic direct-current Tr([Q,C_A] P)", direct_trace)

# These are the flavor factors in the existing RChT derivative interactions.
G_V, f_a = sp.symbols("G_V f_a", nonzero=True, real=True)
tau_PW, tau_Pgamma = sp.symbols("tau_PW tau_Pgamma")
coefficient_PW_before_trace = -2 * sp.I * G_V * tau_PW / (F * f_a)
coefficient_Pgamma_before_trace = -2 * sp.I * G_V * tau_Pgamma / (F * f_a)
formal_C_PW = -2 * sp.I * G_V * tensor_vector_trace / (F * f_a)
formal_C_Pgamma = -2 * sp.I * G_V * mixing_photon_trace / (F * f_a)
show("Formal C_PW before imposing tau_PW=Tr(V[P,C_A])=0", coefficient_PW_before_trace)
show(
    "Formal C_Pgamma before imposing tau_Pgamma=Tr(T_gamma[P,C_A])=0",
    coefficient_Pgamma_before_trace,
)
show("Coefficient C_PW after the generic diagonal trace", formal_C_PW)
show("Coefficient C_Pgamma after the generic diagonal trace", formal_C_Pgamma)


# ---------------------------------------------------------------------------
# 3. Lorentz amplitude from C_X X_{mu nu} partial^mu P partial^nu a
# ---------------------------------------------------------------------------
metric = sp.diag(1, -1, -1, -1)
q_components = sp.symbols("q0 q1 q2 q3")
k_components = sp.symbols("k0 k1 k2 k3")
p_components = sp.symbols("p0 p1 p2 p3")
eps_components = sp.symbols("eps0 eps1 eps2 eps3")
q = sp.Matrix(q_components)
k = sp.Matrix(k_components)
p = sp.Matrix(p_components)
eps = sp.Matrix(eps_components)
q_cov, eps_cov = metric * q, metric * eps
field_strength = q_cov * eps_cov.T - eps_cov * q_cov.T
lorentz_contraction = sp.expand((k.T * field_strength * p)[0])
invariant_contraction = sp.expand(
    minkowski_dot(q, k, metric) * minkowski_dot(eps, p, metric)
    - minkowski_dot(eps, k, metric) * minkowski_dot(q, p, metric)
)
ward_contraction = sp.expand(
    lorentz_contraction.subs(dict(zip(eps_components, q_components)), simultaneous=True)
)

show("Field-strength momentum contraction F_{mu nu} k^mu p^nu", lorentz_contraction)
show("Equivalent invariant contraction", invariant_contraction)
show("Ward replacement epsilon -> q", ward_contraction)

C_X = sp.Symbol("C_X")
m_a, m_P, m_X = sp.symbols("m_a m_P m_X", positive=True, real=True)
epsilon_dot_k = sp.Symbol("epsilon^* dot k")
q_dot_k = (m_a**2 - m_P**2 - m_X**2) / 2
q_dot_p = (m_a**2 - m_P**2 + m_X**2) / 2

# Momentum conservation p=k+q and epsilon.q=0 imply epsilon.p=epsilon.k.
amplitude_before_reduction = C_X * (q_dot_k - q_dot_p) * epsilon_dot_k
amplitude_massive = sp.factor(amplitude_before_reduction)
amplitude_photon = sp.simplify(amplitude_massive.subs(m_X, 0))
show("On-shell amplitude before simplifying q.k-q.p", amplitude_before_reduction)
show("Massive-vector amplitude M(a->P X)", amplitude_massive)
show("Real-photon amplitude M(a->P gamma), q^2=0", amplitude_photon)


# ---------------------------------------------------------------------------
# 4. Polarization sum and two-body width for a massive vector
# ---------------------------------------------------------------------------
lambda_K = sp.expand(
    m_a**4 + m_P**4 + m_X**4
    - 2 * m_a**2 * m_P**2
    - 2 * m_a**2 * m_X**2
    - 2 * m_P**2 * m_X**2
)
polarization_sum_eps_k = sp.factor(-m_P**2 + q_dot_k**2 / m_X**2)
polarization_sum_expected = lambda_K / (4 * m_X**2)
C_X_abs_sq = sp.Abs(C_X) ** 2
amplitude_squared_raw = sp.factor(C_X_abs_sq * m_X**4 * polarization_sum_eps_k)
amplitude_squared = C_X_abs_sq * m_X**2 * lambda_K / 4
three_momentum = sp.sqrt(lambda_K) / (2 * m_a)
width_raw = sp.simplify(three_momentum * amplitude_squared_raw / (8 * sp.pi * m_a**2))
width = C_X_abs_sq * m_X**2 * lambda_K ** sp.Rational(3, 2) / (64 * sp.pi * m_a**3)
lambda_notation = sp.Function("lambda")(m_a**2, m_P**2, m_X**2)
width_compact = C_X_abs_sq * m_X**2 * lambda_notation ** sp.Rational(3, 2) / (
    64 * sp.pi * m_a**3
)

show("Kallen lambda(m_a^2,m_P^2,m_X^2)", lambda_K)
show("Massive-vector polarization sum Sum|epsilon.k|^2", polarization_sum_eps_k)
show("Polarization-summed |M|^2", amplitude_squared)
show("Two-body momentum magnitude", three_momentum)
show("Width Gamma(a->P X) in compact Kallen notation for m_X>0", width_compact)


# ---------------------------------------------------------------------------
# 5. State-level selection rule and exact PDG central-mass thresholds
# ---------------------------------------------------------------------------
# J, P, and C are encoded as exact integers.  The photon needs only J and C
# here because its channels are excluded by the independently checked on-shell
# Ward identity rather than by the massive-vector parity argument.
state_quantum_numbers = {
    "a": {"J": 0, "P": -1, "C": 1},
    "pi0": {"J": 0, "P": -1, "C": 1},
    "eta": {"J": 0, "P": -1, "C": 1},
    "eta_prime": {"J": 0, "P": -1, "C": 1},
    "omega": {"J": 1, "P": -1, "C": -1},
    "phi": {"J": 1, "P": -1, "C": -1},
    "gamma": {"J": 1, "C": -1},
}

J_initial = state_quantum_numbers["a"]["J"]
final_total_spin = state_quantum_numbers["omega"]["J"]
allowed_orbital_L = tuple(
    L_trial
    for L_trial in range(4)
    if abs(L_trial - final_total_spin) <= J_initial <= L_trial + final_total_spin
)
L_PV = allowed_orbital_L[0]
P_initial = state_quantum_numbers["a"]["P"]
C_initial = state_quantum_numbers["a"]["C"]
P_final_PV = (
    state_quantum_numbers["pi0"]["P"]
    * state_quantum_numbers["omega"]["P"]
    * (-1) ** L_PV
)
C_final_PV = (
    state_quantum_numbers["pi0"]["C"]
    * state_quantum_numbers["omega"]["C"]
)
CP_initial = P_initial * C_initial
CP_final_PV = P_final_PV * C_final_PV

show("Angular-momentum solution L for 0 -> 0 + 1", sp.Tuple(*allowed_orbital_L))
show("Initial and final parity (P_a, P_PV)", sp.Tuple(P_initial, P_final_PV))
show("Initial and final charge conjugation (C_a, C_PV)", sp.Tuple(C_initial, C_final_PV))
show("Initial and final CP ((CP)_a, (CP)_PV)", sp.Tuple(CP_initial, CP_final_PV))

# Exact decimal PDG central values converted to rational GeV quantities.
pdg_masses = {
    "pi0": sp.Rational("0.1349768") * GeV,
    "eta": sp.Rational("0.547862") * GeV,
    "eta_prime": sp.Rational("0.95778") * GeV,
    "omega": sp.Rational("0.78266") * GeV,
    "phi": sp.Rational("1.019461") * GeV,
    "gamma": sp.Integer(0) * GeV,
}
for state_name in ("pi0", "eta", "eta_prime", "omega", "phi", "gamma"):
    show(f"PDG central mass m_{state_name} [GeV]", pdg_masses[state_name] / GeV)

thresholds = {
    ("pi0", "gamma"): pdg_masses["pi0"] + pdg_masses["gamma"],
    ("eta", "gamma"): pdg_masses["eta"] + pdg_masses["gamma"],
    ("eta_prime", "gamma"): pdg_masses["eta_prime"] + pdg_masses["gamma"],
    ("pi0", "omega"): pdg_masses["pi0"] + pdg_masses["omega"],
    ("eta", "omega"): pdg_masses["eta"] + pdg_masses["omega"],
    ("eta_prime", "omega"): pdg_masses["eta_prime"] + pdg_masses["omega"],
    ("pi0", "phi"): pdg_masses["pi0"] + pdg_masses["phi"],
    ("eta", "phi"): pdg_masses["eta"] + pdg_masses["phi"],
    ("eta_prime", "phi"): pdg_masses["eta_prime"] + pdg_masses["phi"],
}
memory_table_thresholds = {
    ("pi0", "gamma"): sp.Rational("0.1349768") * GeV,
    ("eta", "gamma"): sp.Rational("0.547862") * GeV,
    ("eta_prime", "gamma"): sp.Rational("0.95778") * GeV,
    ("pi0", "omega"): sp.Rational("0.9176368") * GeV,
    ("eta", "omega"): sp.Rational("1.330522") * GeV,
    ("eta_prime", "omega"): sp.Rational("1.74044") * GeV,
    ("pi0", "phi"): sp.Rational("1.1544378") * GeV,
    ("eta", "phi"): sp.Rational("1.567323") * GeV,
    ("eta_prime", "phi"): sp.Rational("1.977241") * GeV,
}
for (pseudoscalar_name, vector_name), threshold in thresholds.items():
    show(
        f"Threshold m_a for {pseudoscalar_name} {vector_name} [GeV]",
        threshold / GeV,
    )


# ---------------------------------------------------------------------------
# 6. Diagnostics
# ---------------------------------------------------------------------------
assert gram == sp.eye(3)
assert sp.trigsimp(eta_gram - sp.eye(2)) == sp.zeros(2)
assert sp.trigsimp(vector_gram - sp.eye(2)) == sp.zeros(2)
assert sp.simplify(omega.subs(ideal_vector_mixing) - omega_ideal) == sp.zeros(3)
assert sp.simplify(phi.subs(ideal_vector_mixing) - phi_ideal) == sp.zeros(3)
assert tensor_vector_trace == 0
assert mixing_photon_trace == 0
assert direct_trace == 0
for P_concrete in (pi0, eta8, eta0, eta, eta_prime):
    for V_concrete in (omega_ideal, phi_ideal, omega, phi):
        assert sp.simplify(sp.trace(V_concrete * (P_concrete * C_A - C_A * P_concrete))) == 0
    assert sp.simplify(sp.trace((Q * C_A - C_A * Q) * P_concrete)) == 0
assert sp.expand(lorentz_contraction - invariant_contraction) == 0
assert ward_contraction == 0
assert amplitude_massive == -C_X * m_X**2 * epsilon_dot_k
assert amplitude_photon == 0
assert sp.simplify(polarization_sum_eps_k - polarization_sum_expected) == 0
assert sp.simplify(amplitude_squared_raw - amplitude_squared) == 0
assert sp.simplify(width_raw - width) == 0
assert allowed_orbital_L == (1,)
assert P_final_PV == P_initial == -1
assert C_final_PV == -1 != C_initial
assert CP_initial == -1
assert CP_final_PV == 1 != CP_initial
for pseudoscalar_name in ("pi0", "eta", "eta_prime"):
    assert state_quantum_numbers[pseudoscalar_name]["J"] == 0
    assert state_quantum_numbers[pseudoscalar_name]["P"] == -1
    assert state_quantum_numbers[pseudoscalar_name]["C"] == 1
for vector_name in ("omega", "phi"):
    assert state_quantum_numbers[vector_name]["J"] == 1
    assert state_quantum_numbers[vector_name]["P"] == -1
    assert state_quantum_numbers[vector_name]["C"] == -1
    for pseudoscalar_name in ("pi0", "eta", "eta_prime"):
        channel_parity = (
            state_quantum_numbers[pseudoscalar_name]["P"]
            * state_quantum_numbers[vector_name]["P"]
            * (-1) ** L_PV
        )
        channel_C = (
            state_quantum_numbers[pseudoscalar_name]["C"]
            * state_quantum_numbers[vector_name]["C"]
        )
        assert channel_parity == P_initial
        assert channel_C != C_initial
        assert channel_parity * channel_C != CP_initial
assert thresholds == memory_table_thresholds

print("Assumptions and evidence status:")
print("- LITERATURE INPUT [ParticleDataGroup:2024cfk]: meson/photon quantum numbers and quoted PDG central masses.")
print("- DERIVED: L=1; parity is allowed; C and CP forbid all six neutral P omega/phi channels.")
print("- DERIVED: all displayed flavor zeros, the Ward identity, amplitude, width, and nine thresholds.")
print("- INPUT: LO RChT derivative operator and U=exp(i sqrt(2) Phi/F) convention from repository memory.")
print("- CONVENTION: mostly-minus metric; all external momenta are physical and p=k+q.")
print("- CONVENTION: the amplitude is stated up to an irrelevant overall S-matrix phase.")
print("- ASSUMPTION: P and V/T_gamma are diagonal; C_A is an otherwise arbitrary 3x3 spurion.")
print("- ASSUMPTION: the massive-width formula uses m_a>m_P+m_X and m_X>0.")
print("- NO APPROXIMATION: generic eta and omega mixing angles are retained in the flavor proof.")
print("- The real-photon zero follows independently from q^2=0, epsilon.q=0, and the Ward identity.")
print("Units: natural units; f_pi=93/1000 GeV; [C_X]=GeV^(-2); masses and widths are in GeV.")
print("Diagnostic assertions: all passed.")
