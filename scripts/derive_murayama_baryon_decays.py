"""ALP decays into octet baryon pairs, a -> B_1 bar(B_2), in the only-up model.

Run:   .venv/bin/python scripts/derive_murayama_baryon_decays.py
Tests: .venv/bin/python scripts/derive_murayama_baryon_decays.py --test
       (equivalently .venv/bin/python -m unittest scripts.test_derive_murayama_baryon_decays)

Model: only-up phase model of DiLuzio:2026poh with (n_u, n_d, n_s) = (1, 0, 0),
f_a = sqrt(2) m_u / kappa_u, radial field frozen (sigma = 0); inputs from
scripts/implementation/only_up_inputs.py and the physical ALP pole mass from
scripts/implementation/only_up_decays.spectrum.  Conventions of Alda:2026nrd:
Phi = sum_k P_k lambda_k with <lambda_k lambda_l> = 2 delta_kl, U = u^2,
u = exp(i Phi/(2 f_pi)), u_mu = i(u^dag d_mu u - u d_mu u^dag), meson
kinetic term (1/4)<d Phi d Phi>, ALP kinetic term (1/2)(d a)^2.  In this model
c_q^{AP} = diag(1/2, 0, 0), c_GA = 0, c_q^{VS} = 0 and
C = -(1/2){chi_0, c_q^{AP}} = B_0 hat P_01 with hat P_01 = -m_u diag(1, 0, 0).

Basis: the basis without ALP-meson mixing of the draft, U = g_R U' g_R,
g_R = exp(i a B_0 M^{-1}_{m_a^2}(hat P_01)/f_a), followed by the linear completion
a' = a'' - ..., a_phys = sqrt(N(m_a^2)) a''.  The baryon vertex is that of
Alda:2026nrd Sec. V at the physical pole mass, divided by sqrt(N(m_a^2)).

Baryon Lagrangian: LO octet Lagrangian with D_B, F_B, S_B and the singlet
coefficient tilde S_B (Alda:2026nrd Eqs. baryoncoeff, newsinglet).  Baryon
form factors are set to one (pointlike LO vertex).  Section V numerical inputs:
D_B = 0.815, F_B = 0.434, S_B = -0.417, tilde S_B = 0.025 (central value of the
naive estimate, uncertainty 0.031; conservative uncertainty 0.062).
Baryon masses: PDG 2024 (ParticleDataGroup:2024cfk), isospin-split values.

Units: masses, f_a, f_pi and widths in GeV; squared masses in GeV^2; the
couplings c^V, c^A are dimensionless; |M|^2 is dimensionless.
"""

import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import brentq

_IMPLEMENTATION = str(Path(__file__).resolve().parent / "implementation")
if _IMPLEMENTATION not in sys.path:
    sys.path.insert(0, _IMPLEMENTATION)
from only_up_inputs import only_up_inputs  # noqa: E402
from only_up_decays import (  # noqa: E402
    BARYON_MASSES, D_BARYON, F_BARYON, S_BARYON,
    TOPOLOGICAL_SINGLET_BARYON, TOPOLOGICAL_SINGLET_BARYON_CONSERVATIVE_UNCERTAINTY,
    TOPOLOGICAL_SINGLET_BARYON_NAIVE_UNCERTAINTY, alp_baryon_sources, baryon_pair_widths,
    inverse_mass_operator, kallen, neutral_flavour_conserving_channels, octet_generators,
    section_five_couplings, spectrum, spin_summed_squared_amplitude)
from only_up_decays import baryon_pair_width_from_couplings as decay_width  # noqa: E402

# Numerical inputs (D_B, F_B, S_B, tilde S_B, PDG 2024 masses), octet generators,
# M_s^{-1}, Section V couplings and the width formula live in
# scripts/implementation/only_up_decays.py; this script keeps the SymPy derivations
# and the cross-checks in other bases.

KAPPA_GRID = np.logspace(-3, 0, 500)


# --------------------------------------------------------------------------
# Lagrangian couplings (independent check of the Section V traces)
# --------------------------------------------------------------------------

def lagrangian_couplings(first, second, vielbein, connection, topological_gradient,
                         topological_singlet=TOPOLOGICAL_SINGLET_BARYON):
    """Coefficients of (d_mu a/f_a) bar B_1 gamma^mu (c^V + c^A gamma_5) B_2 read off the Lagrangian.

    vielbein, connection, topological_gradient: coefficients of d_mu a/f_a in u_mu,
    Gamma_mu and d_mu X_TC.  Uses <bar B gamma^mu i nabla_mu B> with nabla B = dB + [Gamma, B],
    (D_B/2)<bar B gamma gamma_5 {u, B}>, (F_B/2)<bar B gamma gamma_5 [u, B]>,
    (S_B/2)<bar B gamma gamma_5 B><u>, (tilde S_B/2)<bar B gamma gamma_5 B> d X_TC.
    """
    generators = octet_generators()
    bar_first, second_matrix = generators[first].conj().T, generators[second]
    vector = 1j * np.trace(bar_first @ (connection @ second_matrix - second_matrix @ connection))
    axial = (D_BARYON / 2 * np.trace(bar_first @ (vielbein @ second_matrix + second_matrix @ vielbein))
             + F_BARYON / 2 * np.trace(bar_first @ (vielbein @ second_matrix - second_matrix @ vielbein))
             + S_BARYON / 2 * np.trace(bar_first @ second_matrix) * np.trace(vielbein)
             + topological_singlet / 2 * np.trace(bar_first @ second_matrix) * topological_gradient)
    return complex(vector), complex(axial)


# --------------------------------------------------------------------------
# Only-up model in the basis without ALP-meson mixing
# --------------------------------------------------------------------------

def only_up_alp_data(kappa_u):
    """Physical ALP mass, f_a, M_a^{-1}(C) at the physical pole, and N(m_a^2) for the only-up model."""
    inputs = only_up_inputs(kappa_u)
    alp_mass, axial_source, normalization = alp_baryon_sources(inputs)
    return {"kappa_u": kappa_u, "m_a": alp_mass, "f_a": inputs["f_a"], "f_pi": inputs["f_pi"],
            "B_0": inputs["B_0"], "quark_masses": np.diag(inputs["hat_S"][0]),
            "m_1_squared": inputs["m_0_squared"], "axial_source": axial_source,
            "normalization": normalization, "vector_source": np.zeros((3, 3)),
            "gluon_axial_coefficient": 0.}


def physical_couplings(data, first, second, topological_singlet=TOPOLOGICAL_SINGLET_BARYON):
    """(c^V, c^A) for the canonically normalized ALP: Section V couplings divided by sqrt(N(m_a^2))."""
    vector, axial = section_five_couplings(first, second, data["axial_source"], data["vector_source"],
                                           data["gluon_axial_coefficient"], topological_singlet)
    root = math.sqrt(data["normalization"])
    return vector / root, axial / root


def original_basis_neutral_pole(kappa_u):
    """ALP-like eigenvector of the neutral quadratic form in the original polar fields (pi_3, pi_8, eta_0, a).

    Canonical kinetic terms; returns (m_a^2, meson matrix Phi per unit a_phys, a component).
    """
    inputs = only_up_inputs(kappa_u)
    B_0, f_pi, f_a = inputs["B_0"], inputs["f_pi"], inputs["f_a"]
    generators = [np.diag([1., -1., 0.]), np.diag([1., 1., -2.]) / math.sqrt(3),
                  np.eye(3) * math.sqrt(2. / 3.)]
    # only_up_inputs stores derivatives in a/f_a; the Taylor coefficient hat S_02 is half the second derivative.
    scalar, pseudoscalar, scalar_second = inputs["hat_S"][0], inputs["hat_P"][1], inputs["hat_S"][2] / 2
    mass_squared = np.zeros((4, 4))
    for k, first in enumerate(generators):
        for l, second in enumerate(generators):
            mass_squared[k, l] = (B_0 / 2 * np.trace(scalar @ (first @ second + second @ first))
                                  + inputs["m_0_squared"] / 6 * np.trace(first) * np.trace(second))
        mass_squared[k, 3] = mass_squared[3, k] = -B_0 * f_pi / f_a * np.trace(pseudoscalar @ first)
    mass_squared[3, 3] = -2 * B_0 * f_pi**2 / f_a**2 * np.trace(scalar_second)
    eigenvalues, eigenvectors = np.linalg.eigh(mass_squared)
    target = float(np.real(spectrum(inputs)["a"]))**2
    index = int(np.argmin(abs(eigenvalues - target)))
    vector = eigenvectors[:, index] * np.sign(eigenvectors[3, index])
    meson_matrix = sum(vector[k] * generators[k] for k in range(3))
    return eigenvalues[index], meson_matrix, vector[3]


# --------------------------------------------------------------------------
# Amplitude and width
# --------------------------------------------------------------------------

def channel_widths(kappa_u, topological_singlet=TOPOLOGICAL_SINGLET_BARYON):
    """Widths in GeV of all neutral flavour-conserving channels at kappa_u (only_up_decays.baryon_pair_widths)."""
    return baryon_pair_widths(only_up_inputs(kappa_u), topological_singlet)


def threshold_kappa(first, second):
    """kappa_u at which m_a = m_B1 + m_B2, or None if not reached for kappa_u <= 1."""
    def difference(log_kappa):
        return only_up_alp_data(10**log_kappa)["m_a"] - BARYON_MASSES[first] - BARYON_MASSES[second]
    if difference(0.) <= 0:
        return None
    return 10**brentq(difference, -1., 0., xtol=1e-12)


# --------------------------------------------------------------------------
# Explicit spinors (numerical spin sum)
# --------------------------------------------------------------------------

GAMMA_ZERO = np.diag([1, 1, -1, -1]).astype(complex)
PAULI = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]], dtype=complex)]
GAMMA = [GAMMA_ZERO] + [np.block([[np.zeros((2, 2)), sigma], [-sigma, np.zeros((2, 2))]]) for sigma in PAULI]
GAMMA_FIVE = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]]).astype(complex)


def slash(momentum):
    return GAMMA[0] * momentum[0] - sum(GAMMA[i] * momentum[i] for i in range(1, 4))


def dirac_spinors(momentum, mass):
    """(u(p,s), v(p,s)) for s = 1, 2 in the Dirac representation, normalized to bar u u = 2m."""
    energy, three = momentum[0], np.asarray(momentum[1:])
    sigma_p = sum(PAULI[i] * three[i] for i in range(3))
    root = math.sqrt(energy + mass)
    particles, antiparticles = [], []
    for chi in (np.array([1, 0], dtype=complex), np.array([0, 1], dtype=complex)):
        particles.append(np.concatenate([root * chi, sigma_p @ chi / root]))
        antiparticles.append(np.concatenate([sigma_p @ chi / root, root * chi]))
    return particles, antiparticles


def numerical_spin_sum(vector, axial, first_mass, second_mass, alp_mass, f_a, direction=(0.3, -0.5, 0.81)):
    """sum_spins |bar u_1 qslash (c^V + c^A gamma_5) v_2|^2 / f_a^2 with explicit spinors in the ALP rest frame."""
    momentum = math.sqrt(kallen(alp_mass**2, first_mass**2, second_mass**2)) / (2 * alp_mass)
    unit = np.asarray(direction) / np.linalg.norm(direction)
    first = np.concatenate([[math.sqrt(momentum**2 + first_mass**2)], momentum * unit])
    second = np.concatenate([[math.sqrt(momentum**2 + second_mass**2)], -momentum * unit])
    total = first + second
    vertex = slash(total) @ (vector * np.eye(4) + axial * GAMMA_FIVE)
    particles, _ = dirac_spinors(first, first_mass)
    _, antiparticles = dirac_spinors(second, second_mass)
    result = 0.
    for spinor_one in particles:
        for spinor_two in antiparticles:
            amplitude = spinor_one.conj() @ GAMMA_ZERO @ vertex @ spinor_two / f_a
            result += abs(amplitude)**2
    return result


# --------------------------------------------------------------------------
# Symbolic derivations
# --------------------------------------------------------------------------

def derive_source_identification():
    """C = -(1/2){chi_0, c^{AP}} - (m_1^2/3) c_GA 1 equals B_0 hat P_01 for the only-up sources."""
    B_0, up, down, strange, singlet = sp.symbols("B_0 m_u m_d m_s m_1^2", positive=True)
    chi_zero = 2 * B_0 * sp.diag(up, down, strange)
    axial_minus_pseudoscalar = sp.diag(sp.Rational(1, 2), 0, 0)
    gluon_axial = 0
    invariant = -(chi_zero * axial_minus_pseudoscalar + axial_minus_pseudoscalar * chi_zero) / 2 \
        - singlet / 3 * gluon_axial * sp.eye(3)
    hat_p01 = -up * sp.diag(1, 0, 0)
    return invariant, B_0 * hat_p01


def derive_building_blocks():
    """u_mu, Gamma_mu and X_TC linear in a from u = g_R = exp(i a A/f_a), A = B_0 M^{-1}(hat P_01) diagonal."""
    x = sp.Symbol("x", real=True)
    f_a, f_pi = sp.symbols("f_a f_pi", positive=True)
    alp = sp.Function("a", real=True)(x)
    entries = sp.symbols("A_u A_d A_s", real=True)
    shift = sp.diag(*entries)
    rotation = sp.diag(*[sp.exp(sp.I * alp * entry / f_a) for entry in entries])
    rotation_adjoint = sp.diag(*[sp.exp(-sp.I * alp * entry / f_a) for entry in entries])
    vielbein = sp.simplify(sp.I * (rotation_adjoint * rotation.diff(x) - rotation * rotation_adjoint.diff(x)))
    connection = sp.simplify((rotation_adjoint * rotation.diff(x) + rotation * rotation_adjoint.diff(x)) / 2)
    meson_matrix = 2 * f_pi / f_a * shift * alp  # linear part of Phi = Phi' + 2 f_pi a A/f_a at Phi' = 0
    topological = -meson_matrix.trace() / f_pi
    return {"vielbein": vielbein, "connection": connection, "topological_gradient": sp.diff(topological, x),
            "expected_vielbein": -2 * alp.diff(x) / f_a * shift,
            "expected_topological_gradient": -2 * alp.diff(x) / f_a * shift.trace()}


def derive_kinetic_normalization():
    """Coefficient of (1/2)(d a)^2 after Phi = Phi' + a L, from (1/4)<dPhi dPhi> + (1/2)(da)^2 at Phi' = 0."""
    alp_derivative = sp.Symbol("da", real=True)
    f_a, f_pi, B_0 = sp.symbols("f_a f_pi B_0", positive=True)
    entries = sp.symbols("A_u A_d A_s", real=True)  # B_0 M^{-1}(hat P_01) = M_a^{-1}(C)
    leg = 2 * f_pi / f_a * sp.diag(*entries)
    kinetic = sp.Rational(1, 4) * ((leg * alp_derivative)**2).trace() + alp_derivative**2 / 2
    normalization = sp.expand(2 * kinetic / alp_derivative**2)
    expected = 1 + 2 * f_pi**2 / f_a**2 * (sp.diag(*entries)**2).trace()
    return normalization, expected


def derive_only_up_diagonal_inverse():
    """Closed form of X = M_s^{-1}(C) for diagonal C = -B_0 m_u diag(1,0,0)."""
    B_0, singlet, s = sp.symbols("B_0 m_1^2 s", positive=True)
    masses = sp.symbols("m_u m_d m_s", positive=True)
    unknowns = sp.symbols("X_u X_d X_s")
    invariant = [-B_0 * masses[0], 0, 0]
    trace = sum(unknowns)
    equations = [2 * B_0 * masses[q] * unknowns[q] + singlet / 3 * trace - s * unknowns[q] - invariant[q]
                 for q in range(3)]
    solution = sp.solve(equations, unknowns, dict=True)[0]
    denominators = [2 * B_0 * masses[q] - s for q in range(3)]
    trace_closed = (-B_0 * masses[0] / denominators[0]) / (1 + singlet / 3 * sum(1 / d for d in denominators))
    closed = [(invariant[q] - singlet / 3 * trace_closed) / denominators[q] for q in range(3)]
    return solution, unknowns, closed


def derive_spin_sum():
    """Spin sum tr[(p1slash+m1) Gamma (p2slash-m2) bar Gamma]/f_a^2 with explicit Dirac matrices."""
    m1, m2, alp_mass, f_a = sp.symbols("m_1 m_2 m_a f_a", positive=True)
    vector, axial = sp.symbols("c_V c_A")
    gamma_zero = sp.diag(1, 1, -1, -1)
    pauli = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    zero = sp.zeros(2)
    gammas = [gamma_zero] + [sp.Matrix(sp.BlockMatrix([[zero, sigma], [-sigma, zero]])) for sigma in pauli]
    gamma_five = sp.Matrix(sp.BlockMatrix([[zero, sp.eye(2)], [sp.eye(2), zero]]))
    energy_one = (alp_mass**2 + m1**2 - m2**2) / (2 * alp_mass)
    energy_two = (alp_mass**2 - m1**2 + m2**2) / (2 * alp_mass)
    momentum = sp.Symbol("k", positive=True)
    first = [energy_one, 0, 0, momentum]
    second = [energy_two, 0, 0, -momentum]
    total = [first[i] + second[i] for i in range(4)]

    def feynman_slash(p):
        return gammas[0] * p[0] - sum((gammas[i] * p[i] for i in range(1, 4)), sp.zeros(4))

    vertex = feynman_slash(total) * (vector * sp.eye(4) + axial * gamma_five)
    vertex_bar = gamma_zero * vertex.H * gamma_zero
    trace = ((feynman_slash(first) + m1 * sp.eye(4)) * vertex
             * (feynman_slash(second) - m2 * sp.eye(4)) * vertex_bar).trace() / f_a**2
    momentum_squared = kallen(alp_mass**2, m1**2, m2**2) / (4 * alp_mass**2)
    trace = sp.simplify(sp.expand(trace).subs(momentum**2, momentum_squared))
    closed = 2 / f_a**2 * (axial * sp.conjugate(axial) * (m1 + m2)**2 * (alp_mass**2 - (m1 - m2)**2)
                           + vector * sp.conjugate(vector) * (m1 - m2)**2 * (alp_mass**2 - (m1 + m2)**2))
    width = sp.sqrt(kallen(alp_mass**2, m1**2, m2**2)) / (16 * sp.pi * alp_mass**3) * closed
    return {"trace": trace, "closed": closed, "width": width,
            "symbols": (m1, m2, alp_mass, f_a, vector, axial)}


def main():
    invariant, expected = derive_source_identification()
    print("Only-up invariant source C:", invariant, " B_0 hat P_01:", expected)
    blocks = derive_building_blocks()
    print("u_mu|_a =", blocks["vielbein"])
    print("Gamma_mu|_a =", blocks["connection"])
    print("d_mu X_TC|_a =", blocks["topological_gradient"])
    normalization, _ = derive_kinetic_normalization()
    print("K'_aa = N(m_a^2) =", normalization)
    solution, unknowns, _ = derive_only_up_diagonal_inverse()
    for unknown in unknowns:
        print(f"{unknown} =", sp.factor(solution[unknown]))
    spin = derive_spin_sum()
    print("sum |M|^2 =", sp.factor(spin["closed"]))
    print("Gamma =", spin["width"])
    print()
    print("Inputs: D_B=%.3f F_B=%.3f S_B=%.3f tilde S_B=%.3f (+-%.3f naive); pointlike form factors."
          % (D_BARYON, F_BARYON, S_BARYON, TOPOLOGICAL_SINGLET_BARYON,
             TOPOLOGICAL_SINGLET_BARYON_NAIVE_UNCERTAINTY))
    channels = neutral_flavour_conserving_channels()
    print("Neutral flavour-conserving channels:", channels)
    open_any = []
    for first, second in channels:
        kappa = threshold_kappa(first, second)
        status = f"opens at kappa_u = {kappa:.6f}" if kappa is not None else "closed for kappa_u <= 1"
        print(f"  {first} anti-{second}: threshold {BARYON_MASSES[first] + BARYON_MASSES[second]:.6f} GeV, {status}")
        if kappa is not None:
            open_any.append((first, second))
    grid_open = {pair: 0 for pair in open_any}
    for kappa in KAPPA_GRID:
        widths = channel_widths(kappa)
        for pair in open_any:
            grid_open[pair] += widths[pair] > 0
    print("Number of kappa_u grid points (logspace(-3,0,500)) with open channel:", grid_open)
    print()
    for kappa in (0.85, 0.9, 0.95, 1.0):
        data = only_up_alp_data(kappa)
        print(f"kappa_u={kappa}: m_a={data['m_a']:.6f} GeV, f_a={data['f_a']:.6e} GeV, "
              f"N(m_a^2)={data['normalization']:.6f}, M_a^-1(C)={np.real(np.diag(data['axial_source']))}")
        central = channel_widths(kappa)
        low = channel_widths(kappa, TOPOLOGICAL_SINGLET_BARYON - TOPOLOGICAL_SINGLET_BARYON_NAIVE_UNCERTAINTY)
        high = channel_widths(kappa, TOPOLOGICAL_SINGLET_BARYON + TOPOLOGICAL_SINGLET_BARYON_NAIVE_UNCERTAINTY)
        for pair in open_any:
            vector, axial = physical_couplings(data, *pair)
            if central[pair] > 0:
                print(f"   {pair[0]} anti-{pair[1]}: c^V={vector.real:+.5f} c^A={axial.real:+.5f} "
                      f"Gamma={central[pair]:.5e} GeV  [tilde S_B band: {min(low[pair], high[pair]):.5e}, "
                      f"{max(low[pair], high[pair]):.5e}]")


if __name__ == "__main__":
    if "--test" in sys.argv:
        import unittest
        root = str(Path(__file__).resolve().parent.parent)
        if root not in sys.path:
            sys.path.insert(0, root)
        suite = unittest.defaultTestLoader.loadTestsFromName("scripts.test_derive_murayama_baryon_decays")
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        sys.exit(not result.wasSuccessful())
    main()
