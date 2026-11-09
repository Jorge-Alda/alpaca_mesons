#!/usr/bin/env python3
"""Isospin-breaking VPP operators that can split the charged and neutral K* couplings.

Invocation:
    .venv/bin/python scripts/derive_isospin_breaking_vpp_operators.py

Derivation only: no refit, no change of the baseline likelihood. Units: GeV for
masses, widths, G_V and F; GeV^2 for the vacuum spurion chi_+ = 4 B_0 m_q.
Conventions of Alda:2026nrd and of scripts/fit_vector_gv.py: U = u^2 =
exp(i sqrt2 Phi/F), u_mu = -sqrt2 d_mu Phi/F + O(Phi^3), the LO interaction
i G_V/sqrt2 <V_{mu nu} u^mu u^nu>, the squared flavor trace
tr(T_V [T_P1, T_P2])^2 of `flavor_factors`, and the width of
`width_coefficient`. The flavor matrices of the decaying vector and of the
two pseudoscalars are those of `flavor_factors` (published pi0 diagonals of
Alda:2026nrd, nonideal omega-phi angle of the provisional fit).

Operators added to the LO term (coefficient names proposed, awaiting user
approval; in the code: `coefficient_anticommutator`, `coefficient_sandwich`,
`coefficient_vector_singlet`):
    i lambda_1 <V_{mu nu} {u^mu u^nu, chi_+}>
  + i lambda_2 <V_{mu nu} u^mu chi_+ u^nu>
  + i lambda_3 <V_{mu nu}> <u^mu u^nu chi_+>.
At two pseudoscalars and constant quark masses the VPP vertex keeps the LO
momentum structure, with the flavor trace replaced by
    T = tr(V[P1,P2]) + (sqrt2/G_V) [ lambda_1 tr(V{[P1,P2],X})
        + lambda_2 tr(V(P1 X P2 - P2 X P1)) + lambda_3 tr(V) tr([P1,P2] X) ],
with X = chi_+ (vacuum) = 4 diag(B0 m_u, B0 m_d, B0 m_s).
The script works with the ratios r_i = sqrt2 lambda_i / G_V (GeV^-2).

Steps (all executed with assertions):
  1. Enumerate C-even, two-pseudoscalar flavor structures linear in a diagonal
     spurion (all single and multiple traces); compute their ranks on random
     3x3 and 5x5 matrices; check that the chi_- operator
     <V[u^mu, nabla^nu chi_-]> reduces to the single-trace pair.
  2. Symbolic channel traces; at r_i = 0 they reproduce `flavor_factors`.
  3. Exact charged/neutral K* width ratio as a function of r_1 or r_2 alone;
     solve for the PDG-matching ratio kappa_c/kappa_n.
  4. Predicted phi -> K+K-/K0K0bar, rho+/rho0 and omega -> pi pi effects,
     and the SU(3)-breaking side effects of each solution.
  5. Electromagnetic single-trace structures with two charge matrices Q:
     channel pattern and required size compared with alpha/pi (estimate).

Output summary (printed; values of the 2026-11-09 run):
  * 5 independent C-even structures linear in a diagonal spurion (2 single
    traces), equal ranks for 3x3 and 5x5 matrices; the chi_- operator is
    redundant at two pseudoscalars.
  * r_i = 0 reproduces flavor_factors to 1.1e-16.
  * linear isospin-limit splits: K* (charged-neutral) = (r1 - r2/3)(X_u-X_d),
    phi (K+K- minus K0K0bar) = (r2 - r3)(X_u-X_d), rho+ minus rho0 = 0.
  * required Gamma_c/Gamma_n relative to r=0: 1.094483.
    r1 alone: r1 = -0.6951 GeV^-2, r1 X_s = -0.656 (K* coupling/rho = 0.34).
    r2 alone: r2 = +3.7362 GeV^-2, r2 X_s = +3.53; phi rate ratio -16.3%
    (preferred -4.4 +- 2.0%); rho, omega -> K Kbar factor 3.85 relative to K*.
  * electromagnetic single traces (4 independent) need 58 to 178 alpha/pi.
"""
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_vector_gv as gv_fit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MEASUREMENTS_PATH = ROOT / 'data/vector_sector/gv_measurements.json'
FINE_STRUCTURE_CONSTANT = 1/137.035999084
# Required ratio quoted in joint-vector-kstar-kpi-saturation-freedom.md (Step 2).
RECORDED_KAPPA_CHARGED, RECORDED_KAPPA_NEUTRAL = 1.147715, 1.048637
# Endpoint phi predictions from joint-vector-phi-kaon-residual-diagnostic.md.
RECORDED_PHI_PREDICTIONS = (0.5084953, 0.3355123)

E = gv_fit.matrix_unit


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def alda_inputs():
    """Return the Alda:2026nrd input dictionary (B0 m_q, pi0 diagonals, masses)."""
    return json.loads(gv_fit.ALDA_PATH.read_text())


def vacuum_spurion_values():
    """Return 4 B0 m_q for u, d, s in GeV^2 (chi_+ at zero pseudoscalar fields)."""
    return [4*value for value in alda_inputs()['B0mq_gev2']]


def measurements():
    """Return the PDG rows and covariance used by the conditional G_V fit."""
    document = json.loads(MEASUREMENTS_PATH.read_text())
    rows = {row['channel']: row for row in document['rows']}
    return rows, np.array(document['covariance'])


# ---------------------------------------------------------------------------
# Step 1: enumeration of flavor structures
# ---------------------------------------------------------------------------
def monomial_value(cycles, matrices):
    """Evaluate a product of traces; each cycle is a tuple of letters."""
    value = 1
    for cycle in cycles:
        product = np.eye(matrices[cycle[0]].shape[0], dtype=complex)
        for letter in cycle:
            product = product @ matrices[letter]
        value *= np.trace(product)
    return value


def charge_conjugate(cycles):
    """C image of a trace monomial: every trace reversed, overall sign -1 from V.

    V -> -V^T, u_mu -> u_mu^T, chi_+ -> chi_+^T (Q -> Q^T); the Lorentz index
    carried by each u stays attached to its pseudoscalar, so the reversed
    monomial keeps the labels P1 (index mu) and P2 (index nu).
    """
    return -1, tuple(tuple(reversed(cycle)) for cycle in cycles)


def antisymmetrized_c_even_value(cycles, matrices):
    """Value of (1 + C) applied to the monomial, antisymmetrized in P1 <-> P2."""
    sign, image = charge_conjugate(cycles)
    swapped = dict(matrices, P1=matrices['P2'], P2=matrices['P1'])
    def antisymmetric(monomial):
        return monomial_value(monomial, matrices) - monomial_value(monomial, swapped)
    return antisymmetric(cycles) + sign*antisymmetric(image)


def trace_monomials(letters):
    """All products of traces of the given letters (one per permutation)."""
    seen, result = set(), []
    for permutation in itertools.permutations(range(len(letters))):
        remaining, cycles = set(range(len(letters))), []
        while remaining:
            start = min(remaining)
            cycle, position = [], start
            while position in remaining:
                remaining.discard(position)
                cycle.append(letters[position])
                position = permutation[position]
            # canonical cyclic rotation
            rotations = [tuple(cycle[k:]+cycle[:k]) for k in range(len(cycle))]
            cycles.append(min(rotations))
        key = tuple(sorted(cycles))
        if key not in seen:
            seen.add(key)
            result.append(key)
    return result


def structure_rank(letters, dimension, sample_count=40, single_trace_only=False, seed=20261109):
    """Rank of the C-even antisymmetrized structures on random matrices.

    The spurion letters (S or Q) are diagonal, as the vacuum chi_+ and Q are.
    """
    generator = np.random.default_rng(seed)
    monomials = [m for m in trace_monomials(letters)
                 if not single_trace_only or len(m) == 1]
    rows = []
    for _ in range(sample_count):
        matrices = {name: generator.normal(size=(dimension, dimension))
                    + 1j*generator.normal(size=(dimension, dimension)) for name in ('V', 'P1', 'P2')}
        diagonal = np.diag(generator.normal(size=dimension))
        for name in set(letters) - {'V', 'P1', 'P2'}:
            matrices[name] = diagonal
        rows.append([antisymmetrized_c_even_value(m, matrices) for m in monomials])
    return int(np.linalg.matrix_rank(np.array(rows), tol=1e-9)), monomials


def independent_structures(letters, dimension=5, seed=7):
    """Greedy basis of C-even antisymmetrized monomials, fewest traces first."""
    generator = np.random.default_rng(seed)
    monomials = sorted(trace_monomials(letters), key=len)
    samples = []
    for _ in range(40):
        matrices = {n: generator.normal(size=(dimension, dimension)) for n in ('V', 'P1', 'P2')}
        diagonal = np.diag(generator.normal(size=dimension))
        for name in set(letters) - {'V', 'P1', 'P2'}:
            matrices[name] = diagonal
        samples.append([antisymmetrized_c_even_value(m, matrices).real for m in monomials])
    samples = np.array(samples)
    chosen = []
    for index in range(len(monomials)):
        if np.linalg.matrix_rank(samples[:, chosen + [index]], tol=1e-9) == len(chosen) + 1:
            chosen.append(index)
    return [monomials[index] for index in chosen]


def readable(monomial):
    """Product of traces as text, each trace rotated to start with V if present."""
    parts = []
    for cycle in monomial:
        if 'V' in cycle:
            cycle = cycle[cycle.index('V'):] + cycle[:cycle.index('V')]
        parts.append('<' + ' '.join(cycle) + '>')
    return ''.join(parts)


def check_chi_minus_reduction():
    """[P1,{P2,S}] - [P2,{P1,S}] = {[P1,P2],S} + 2(P1 S P2 - P2 S P1).

    At two pseudoscalars u^mu ~ d^mu Phi and nabla^nu chi_- ~ {d^nu Phi, M}, so
    <V_{mu nu}[u^mu, nabla^nu chi_-]> antisymmetrized in the two fields is a
    fixed combination of the two single-trace structures. Exact with SymPy.
    """
    first, second = sp.Matrix(3, 3, sp.symbols('a0:9')), sp.Matrix(3, 3, sp.symbols('b0:9'))
    spurion = sp.diag(*sp.symbols('s0:3'))
    def commutator(x, y):
        return x*y-y*x
    def anticommutator(x, y):
        return x*y+y*x
    left = commutator(first, anticommutator(second, spurion)) - commutator(second, anticommutator(first, spurion))
    right = anticommutator(commutator(first, second), spurion) + 2*(first*spurion*second - second*spurion*first)
    assert sp.expand(left-right) == sp.zeros(3)
    return True


# ---------------------------------------------------------------------------
# Step 2: symbolic channel traces
# ---------------------------------------------------------------------------
SYMBOLS = dict(
    r1=sp.Symbol('r_1', real=True), r2=sp.Symbol('r_2', real=True), r3=sp.Symbol('r_3', real=True),
    Xu=sp.Symbol('X_u', positive=True), Xd=sp.Symbol('X_d', positive=True),
    Xs=sp.Symbol('X_s', positive=True), theta=sp.Symbol('theta', real=True),
    a=sp.Symbol('a', real=True), b=sp.Symbol('b', real=True), c=sp.Symbol('c', real=True))


def channel_matrices(pion_zero_diagonal=None, mixing_angle=None):
    """Return {channel: (V, P1, P2)} with the matrices of `flavor_factors`."""
    a, b, c, theta = (SYMBOLS[k] for k in ('a', 'b', 'c', 'theta'))
    if pion_zero_diagonal is not None:
        a, b, c = pion_zero_diagonal
    if mixing_angle is not None:
        theta = mixing_angle
    pion_zero = sp.diag(a, b, c)/sp.sqrt(2)
    omega_ideal, phi_ideal = sp.diag(1, 1, 0)/sp.sqrt(2), sp.diag(0, 0, 1)
    phi = sp.sin(theta)*omega_ideal + sp.cos(theta)*phi_ideal
    omega = sp.cos(theta)*omega_ideal - sp.sin(theta)*phi_ideal
    rho_zero = sp.diag(1, -1, 0)/sp.sqrt(2)
    return {
        'kstar_charged_K+_pi0': (E(2, 0), E(0, 2), pion_zero),
        'kstar_charged_K0_pi+': (E(2, 0), E(1, 2), E(0, 1)),
        'kstar_neutral_K0_pi0': (E(2, 1), E(1, 2), pion_zero),
        'kstar_neutral_K+_pi-': (E(2, 1), E(0, 2), E(1, 0)),
        'phi_K+_K-': (phi, E(0, 2), E(2, 0)),
        'phi_K0_K0bar': (phi, E(1, 2), E(2, 1)),
        'rho0_pi+_pi-': (rho_zero, E(0, 1), E(1, 0)),
        'rho+_pi+_pi0': (E(1, 0), E(0, 1), pion_zero),
        'omega_pi+_pi-': (omega, E(0, 1), E(1, 0)),
        'phi_pi+_pi-': (phi, E(0, 1), E(1, 0)),
        'rho0_K+_K-': (rho_zero, E(0, 2), E(2, 0)),
        'omega_K+_K-': (omega, E(0, 2), E(2, 0)),
        'rho0_K0_K0bar': (rho_zero, E(1, 2), E(2, 1)),
        'omega_K0_K0bar': (omega, E(1, 2), E(2, 1)),
    }


def structure_traces(vector, first, second, spurion):
    """LO trace and the three mass-spurion traces (before the r_i factors)."""
    commutator = first*second - second*first
    return dict(
        leading=sp.trace(vector*commutator),
        anticommutator=sp.trace(vector*(commutator*spurion + spurion*commutator)),
        sandwich=sp.trace(vector*(first*spurion*second - second*spurion*first)),
        vector_singlet=sp.trace(vector)*sp.trace(commutator*spurion))


def effective_traces(pion_zero_diagonal=None, mixing_angle=None):
    """Effective trace T(r_1, r_2, r_3) for every channel (dimensionless)."""
    spurion = sp.diag(SYMBOLS['Xu'], SYMBOLS['Xd'], SYMBOLS['Xs'])
    result = {}
    for name, (vector, first, second) in channel_matrices(pion_zero_diagonal, mixing_angle).items():
        parts = structure_traces(vector, first, second, spurion)
        result[name] = sp.expand(parts['leading'] + SYMBOLS['r1']*parts['anticommutator']
                                 + SYMBOLS['r2']*parts['sandwich'] + SYMBOLS['r3']*parts['vector_singlet'])
    return result


def isospin_limit_relative_shifts():
    """Linear relative vertex shifts dT/T per unit r_i, isospin-limit pi0, ideal mixing."""
    traces = effective_traces((1, -1, 0), 0)
    shifts = {}
    for name, trace in traces.items():
        leading = trace.subs({SYMBOLS['r1']: 0, SYMBOLS['r2']: 0, SYMBOLS['r3']: 0})
        if leading == 0:
            continue
        shifts[name] = {key: sp.simplify(sp.diff(trace, SYMBOLS[key])/leading) for key in ('r1', 'r2', 'r3')}
    return shifts


def check_zero_coefficient_reduction():
    """At r_i = 0 the squared traces equal `flavor_factors` (published pi0, fit angle)."""
    _, _, mixing_angle = gv_fit.fixed_parameters()
    published = alda_inputs()['published_lambda_diagonals']['pi0']
    reference = gv_fit.flavor_factors(mixing_angle)
    traces = effective_traces([sp.Float(x, 30) for x in published], mixing_angle)
    zero = {SYMBOLS['r1']: 0, SYMBOLS['r2']: 0, SYMBOLS['r3']: 0}
    differences = {}
    for name in reference:
        value = float(sp.N(traces[name].subs(zero)**2, 30))
        differences[name] = abs(value - reference[name])
        assert differences[name] < 1e-13, (name, value, reference[name])
    # exact isospin-limit values of the saturation record
    exact = effective_traces((1, -1, 0), 0)
    assert exact['kstar_charged_K+_pi0'].subs(zero)**2 == sp.Rational(1, 2)
    assert exact['kstar_neutral_K+_pi-'].subs(zero)**2 == 1
    return differences


# ---------------------------------------------------------------------------
# Step 3: K* charged/neutral width ratio and the coefficient solve
# ---------------------------------------------------------------------------
KSTAR_MODES = {
    'charged': [('kstar_charged_K+_pi0', 'K_charged', 'pi0'), ('kstar_charged_K0_pi+', 'K_neutral', 'pi_charged')],
    'neutral': [('kstar_neutral_K0_pi0', 'K_neutral', 'pi0'), ('kstar_neutral_K+_pi-', 'K_charged', 'pi_charged')],
}


def kstar_width_ratio_function(spurion_values=None):
    """Return R(r_1, r_2, r_3) = [Gamma_c(r)/Gamma_c(0)] / [Gamma_n(r)/Gamma_n(0)].

    Widths: width_coefficient with the squared effective trace; PDG parent
    masses of the G_V measurement rows; the common tensor mass cancels.
    """
    _, mixing_angle = None, gv_fit.fixed_parameters()[2]
    inputs = alda_inputs()
    masses = inputs['physical_masses_gev']
    rows, _ = measurements()
    parent = {'charged': rows['kstar_charged_Kpi']['physical_vector_mass_GeV'],
              'neutral': rows['kstar_neutral_Kpi']['physical_vector_mass_GeV']}
    traces = effective_traces([sp.Float(x, 30) for x in inputs['published_lambda_diagonals']['pi0']], mixing_angle)
    spurion_values = spurion_values or vacuum_spurion_values()
    numbers = dict(zip((SYMBOLS['Xu'], SYMBOLS['Xd'], SYMBOLS['Xs']), spurion_values))
    zero = {SYMBOLS['r1']: 0, SYMBOLS['r2']: 0, SYMBOLS['r3']: 0}
    widths = {}
    for charge, modes in KSTAR_MODES.items():
        total = 0
        for channel, first, second in modes:
            phase_space = gv_fit.width_coefficient(parent[charge], 1.0, masses[first], masses[second], 1.0)
            total += phase_space*traces[channel].subs(numbers)**2
        widths[charge] = sp.expand(total)
    ratio = (widths['charged']/widths['charged'].subs(zero))/(widths['neutral']/widths['neutral'].subs(zero))
    return ratio, widths


def required_kstar_ratio():
    """kappa_c/kappa_n from PDG widths and the saturated (r=0) width ratio."""
    rows, _ = measurements()
    _, widths = kstar_width_ratio_function()
    zero = {SYMBOLS['r1']: 0, SYMBOLS['r2']: 0, SYMBOLS['r3']: 0}
    saturated = float(widths['charged'].subs(zero)/widths['neutral'].subs(zero))
    measured = rows['kstar_charged_Kpi']['value']/rows['kstar_neutral_Kpi']['value']
    value = measured/saturated
    assert abs(value - RECORDED_KAPPA_CHARGED/RECORDED_KAPPA_NEUTRAL) < 2e-6
    return value


def solve_single_coefficient(name, target):
    """Exact real roots of R(r)=target with only r_name nonzero; return root nearest zero."""
    ratio, _ = kstar_width_ratio_function()
    others = {SYMBOLS[k]: 0 for k in ('r1', 'r2', 'r3') if k != name}
    numerator, denominator = sp.fraction(sp.together(ratio.subs(others) - target))
    roots = [complex(root) for root in sp.Poly(sp.expand(numerator), SYMBOLS[name]).nroots(n=30)]
    real_roots = sorted((root.real for root in roots if abs(root.imag) < 1e-20), key=abs)
    return real_roots[0], real_roots


# ---------------------------------------------------------------------------
# Step 4: predictions for phi, rho, omega at a given coefficient
# ---------------------------------------------------------------------------
def relative_couplings(coefficients):
    """T(r)/T(0) per channel with published pi0, fit angle and Alda spurions.

    For channels that vanish at r=0 (omega, phi -> pi+ pi-) return T(r)
    divided by the rho0 -> pi+ pi- LO trace.
    """
    traces = numeric_traces()
    values = {SYMBOLS[k]: coefficients.get(k, 0) for k in ('r1', 'r2', 'r3')}
    zero = {SYMBOLS['r1']: 0, SYMBOLS['r2']: 0, SYMBOLS['r3']: 0}
    rho_leading = float(traces['rho0_pi+_pi-'].subs(zero))
    result = {}
    for name, trace in traces.items():
        leading = float(trace.subs(zero))
        value = float(trace.subs(values))
        result[name] = value/leading if abs(leading) > 1e-12 else value/rho_leading
    return result


def phi_ratio_comparison(phi_rate_ratio_change):
    """Compare a relative change of Gamma(K+K-)/Gamma(K0K0bar) with the measured one."""
    rows, covariance = measurements()
    charged, neutral = rows['phi_K+_K-']['value'], rows['phi_KL_KS']['value']
    variance = covariance[2, 2]/charged**2 + covariance[3, 3]/neutral**2 - 2*covariance[2, 3]/(charged*neutral)
    predicted_ratio = RECORDED_PHI_PREDICTIONS[0]/RECORDED_PHI_PREDICTIONS[1]
    preferred_change = (charged/neutral)/predicted_ratio - 1
    return dict(preferred_change=preferred_change, relative_sigma=math.sqrt(variance),
                predicted_change=phi_rate_ratio_change,
                pull_after=(predicted_ratio*(1+phi_rate_ratio_change) - charged/neutral)/(charged/neutral*math.sqrt(variance)))


def check_phi_leading_ratio():
    """The LO phi ratio is pure phase space and reproduces the endpoint predictions."""
    masses = alda_inputs()['physical_masses_gev']
    rows, _ = measurements()
    parent = rows['phi_K+_K-']['physical_vector_mass_GeV']
    charged = gv_fit.width_coefficient(parent, 1, masses['K_charged'], masses['K_charged'], 1)
    neutral = gv_fit.width_coefficient(parent, 1, masses['K_neutral'], masses['K_neutral'], 1)
    recorded = RECORDED_PHI_PREDICTIONS[0]/RECORDED_PHI_PREDICTIONS[1]
    assert abs(charged/neutral/recorded - 1) < 1e-6
    return charged/neutral


def linear_isospin_limit_splits():
    """Linear K* and phi charge splits in the isospin limit (SymPy).

    K* split: rate-weighted dT/T of the charged minus the neutral K*, with the
    isospin-limit weights 1/3 (K pi0) and 2/3 (charged pion), i.e. equal phase
    space. phi split: dT/T of K+K- minus K0K0bar. rho split: rho+ minus rho0.
    """
    shifts = isospin_limit_relative_shifts()
    def total(name):
        return sum(SYMBOLS[k]*shifts[name][k] for k in ('r1', 'r2', 'r3'))
    third = sp.Rational(1, 3)
    kstar = sp.factor(third*total('kstar_charged_K+_pi0') + 2*third*total('kstar_charged_K0_pi+')
                      - third*total('kstar_neutral_K0_pi0') - 2*third*total('kstar_neutral_K+_pi-'))
    phi = sp.factor(total('phi_K+_K-') - total('phi_K0_K0bar'))
    rho = sp.simplify(total('rho+_pi+_pi0') - total('rho0_pi+_pi-'))
    r1, r2, Xu, Xd = (SYMBOLS[k] for k in ('r1', 'r2', 'Xu', 'Xd'))
    assert sp.simplify(kstar - (r1 - r2/3)*(Xu - Xd)) == 0
    assert sp.simplify(phi - (r2 - SYMBOLS['r3'])*(Xu - Xd)) == 0
    assert rho == 0
    return dict(kstar=kstar, phi=phi, rho=rho)


def numeric_traces():
    """Effective traces with published pi0, fit angle and Alda spurions; r_i symbolic."""
    _, _, mixing_angle = gv_fit.fixed_parameters()
    inputs = alda_inputs()
    traces = effective_traces([sp.Float(x, 30) for x in inputs['published_lambda_diagonals']['pi0']], mixing_angle)
    numbers = dict(zip((SYMBOLS['Xu'], SYMBOLS['Xd'], SYMBOLS['Xs']), vacuum_spurion_values()))
    return {name: trace.subs(numbers) for name, trace in traces.items()}


def phi_rate_ratio_change_expression():
    """[Gamma(K+K-)/Gamma(K0K0bar)](r) / (same at r=0) - 1, symbolic in r_i."""
    traces = numeric_traces()
    zero = {SYMBOLS['r1']: 0, SYMBOLS['r2']: 0, SYMBOLS['r3']: 0}
    charged, neutral = traces['phi_K+_K-'], traces['phi_K0_K0bar']
    return (charged/neutral)**2/(charged.subs(zero)/neutral.subs(zero))**2 - 1


def solve_combined(free_names, target, phi_change, start):
    """Solve R_K*(r) = target and phi rate-ratio change = phi_change for two r_i.

    The third coefficient is zero. Uses sympy.nsolve from the given start.
    """
    ratio, _ = kstar_width_ratio_function()
    others = {SYMBOLS[k]: 0 for k in ('r1', 'r2', 'r3') if k not in free_names}
    equations = [ratio.subs(others) - target, phi_rate_ratio_change_expression().subs(others) - phi_change]
    unknowns = [SYMBOLS[k] for k in free_names]
    solution = sp.nsolve(equations, unknowns, start, tol=1e-25, prec=30)
    return {k: float(v) for k, v in zip(free_names, solution)}


def scenario_summary(coefficients):
    """Charge splits and SU(3)-breaking side effects at given r_i."""
    relative = relative_couplings(coefficients)
    ratio, _ = kstar_width_ratio_function()
    kstar_ratio = float(ratio.subs({SYMBOLS[k]: coefficients.get(k, 0) for k in ('r1', 'r2', 'r3')}))
    phi_change = (relative['phi_K+_K-']/relative['phi_K0_K0bar'])**2 - 1
    return dict(
        kstar_width_ratio=kstar_ratio,
        phi_rate_ratio_change=phi_change,
        phi=phi_ratio_comparison(phi_change),
        rho_plus_over_rho_zero_coupling=relative['rho+_pi+_pi0']/relative['rho0_pi+_pi-'],
        kstar_charged_over_rho_coupling=relative['kstar_charged_K0_pi+']/relative['rho0_pi+_pi-'],
        phi_over_rho_coupling=relative['phi_K0_K0bar']/relative['rho0_pi+_pi-'],
        omega_pipi_over_rho_leading=relative['omega_pi+_pi-'],
        phi_pipi_over_rho_leading=relative['phi_pi+_pi-'],
        phi_over_kstar_coupling=relative['phi_K0_K0bar']/relative['kstar_charged_K0_pi+'],
        light_vector_to_kaons_over_kstar={name: relative[name]/relative['kstar_charged_K0_pi+'] for name in
                                          ('rho0_K+_K-', 'omega_K+_K-', 'rho0_K0_K0bar', 'omega_K0_K0bar')},
        modes=relative)


# ---------------------------------------------------------------------------
# Step 5: electromagnetic single-trace structures with two charge matrices
# ---------------------------------------------------------------------------
def electromagnetic_structures():
    """Independent C-even single traces of {V,P1,P2,Q,Q} and their channel shifts.

    Returns, for each representative single-trace word, the linear relative
    shift dT/T per unit coefficient in the K* modes, the charged-minus-neutral
    rate-weighted K* split, and the phi K+K- minus K0K0bar split (isospin-limit
    pi0, ideal mixing). Q = diag(2/3, -1/3, -1/3).
    """
    charge = np.diag([2/3, -1/3, -1/3])
    pion_zero = np.diag([1, -1, 0])/math.sqrt(2)
    def unit(i, j):
        m = np.zeros((3, 3))
        m[i, j] = 1
        return m
    channels = {
        'kstar_charged_K+_pi0': (unit(2, 0), unit(0, 2), pion_zero),
        'kstar_charged_K0_pi+': (unit(2, 0), unit(1, 2), unit(0, 1)),
        'kstar_neutral_K0_pi0': (unit(2, 1), unit(1, 2), pion_zero),
        'kstar_neutral_K+_pi-': (unit(2, 1), unit(0, 2), unit(1, 0)),
        'phi_K+_K-': (unit(2, 2), unit(0, 2), unit(2, 0)),
        'phi_K0_K0bar': (unit(2, 2), unit(1, 2), unit(2, 1)),
        'rho0_pi+_pi-': (np.diag([1, -1, 0])/math.sqrt(2), unit(0, 1), unit(1, 0)),
        'rho+_pi+_pi0': (unit(1, 0), unit(0, 1), pion_zero),
    }
    rank, monomials = structure_rank(('V', 'P1', 'P2', 'Q', 'Q'), 5, single_trace_only=True)
    # choose a representative set of words with linearly independent values
    generator = np.random.default_rng(7)
    samples = []
    for _ in range(40):
        mats = {n: generator.normal(size=(5, 5)) for n in ('V', 'P1', 'P2')}
        mats['Q'] = np.diag(generator.normal(size=5))
        samples.append([antisymmetrized_c_even_value(m, mats).real for m in monomials])
    samples = np.array(samples)
    chosen = []
    for index in range(len(monomials)):
        trial = chosen + [index]
        if np.linalg.matrix_rank(samples[:, trial], tol=1e-9) == len(trial):
            chosen = trial
    assert len(chosen) == rank
    leading = {name: antisymmetrized_c_even_value(((('V', 'P1', 'P2'),)), dict(V=v, P1=p, P2=q)).real
               for name, (v, p, q) in channels.items()}
    weights = {}
    masses = alda_inputs()['physical_masses_gev']
    rows, _ = measurements()
    for charge_name, modes in KSTAR_MODES.items():
        parent = rows[f'kstar_{charge_name}_Kpi']['physical_vector_mass_GeV']
        partial = [gv_fit.width_coefficient(parent, 1, masses[f], masses[s], leading[c]**2) for c, f, s in modes]
        weights[charge_name] = {c: w/sum(partial) for (c, _, _), w in zip(modes, partial)}
    result = []
    for index in chosen:
        monomial = monomials[index]
        shifts = {name: antisymmetrized_c_even_value(monomial, dict(V=v, P1=p, P2=q, Q=charge)).real/leading[name]
                  for name, (v, p, q) in channels.items()}
        kstar_split = (sum(w*shifts[c] for c, w in weights['charged'].items())
                       - sum(w*shifts[c] for c, w in weights['neutral'].items()))
        word = monomial[0][monomial[0].index('V'):] + monomial[0][:monomial[0].index('V')]
        result.append(dict(word=' '.join(word), shifts=shifts, kstar_split=kstar_split,
                           phi_split=shifts['phi_K+_K-'] - shifts['phi_K0_K0bar'],
                           rho_split=shifts['rho+_pi+_pi0'] - shifts['rho0_pi+_pi-']))
    return rank, result


# ---------------------------------------------------------------------------
def main():
    print('Step 1: structure ranks (C-even, antisymmetric in the two pseudoscalars)')
    for dimension in (3, 5):
        rank, monomials = structure_rank(('V', 'P1', 'P2', 'S'), dimension)
        single, _ = structure_rank(('V', 'P1', 'P2', 'S'), dimension, single_trace_only=True)
        print(f'  {dimension}x{dimension}: all traces rank {rank} of {len(monomials)} monomials; single traces rank {single}')
    print('  independent mass-spurion structures (each taken with its C image, antisymmetrized):')
    for monomial in independent_structures(('V', 'P1', 'P2', 'S')):
        print('   ', readable(monomial))
    print('  chi_- reduction identity holds:', check_chi_minus_reduction())

    print('Step 2: zero-coefficient reduction, max |difference| =',
          max(check_zero_coefficient_reduction().values()))
    for name, shift in isospin_limit_relative_shifts().items():
        print(f'  {name:24s} dT/T: r1 {shift["r1"]},  r2 {shift["r2"]},  r3 {shift["r3"]}')

    splits = linear_isospin_limit_splits()
    print('  linear isospin-limit splits: K* charged-neutral', splits['kstar'], '; phi K+K- minus K0K0bar',
          splits['phi'], '; rho+ minus rho0', splits['rho'])

    target = required_kstar_ratio()
    xu, xd, xs = vacuum_spurion_values()
    print(f'Step 3: required Gamma_c/Gamma_n relative to r=0: {target:.6f} '
          f'(coupling ratio {math.sqrt(target):.6f})')
    print(f'  X = 4B0 m_q = ({xu:.6f}, {xd:.6f}, {xs:.6f}) GeV^2; (m_d-m_u)/m_s = {(xd-xu)/xs:.5f}; '
          f'alpha/pi = {FINE_STRUCTURE_CONSTANT/math.pi:.5f}')
    for name in ('r1', 'r2'):
        root, roots = solve_single_coefficient(name, target)
        summary = scenario_summary({name: root})
        print(f'  {name} only: r = {root:.6f} GeV^-2 (all real roots {[round(x, 4) for x in roots]})')
        print(f'    r (X_d - X_u) = {root*(xd-xu):+.5f},  r X_s = {root*xs:+.4f}')
        print(f'    check R = {summary["kstar_width_ratio"]:.6f}')
        print(f'    phi K+K-/K0K0bar rate ratio change {summary["phi_rate_ratio_change"]:+.4f}; '
              f'preferred {summary["phi"]["preferred_change"]:+.4f} +- {summary["phi"]["relative_sigma"]:.4f}; '
              f'pull after {summary["phi"]["pull_after"]:+.2f}')
        print(f'    rho+/rho0 coupling {summary["rho_plus_over_rho_zero_coupling"]:.6f}; '
              f'K*(K0pi+)/rho0 {summary["kstar_charged_over_rho_coupling"]:.4f}; '
              f'phi(K0K0bar)/rho0 {summary["phi_over_rho_coupling"]:.4f}; '
              f'omega->pipi/rho LO {summary["omega_pipi_over_rho_leading"]:+.5f}; '
              f'phi->pipi/rho LO {summary["phi_pipi_over_rho_leading"]:+.5f}')
        print(f'    phi(K0K0bar)/K*(K0pi+) {summary["phi_over_kstar_coupling"]:.4f}; (rho0, omega -> KK)/K*(K0pi+): '
              + ', '.join(f'{k} {v:.4f}' for k, v in summary['light_vector_to_kaons_over_kstar'].items()))
        print('    K* mode factors T(r)/T(0): ' + ', '.join(f'{k} {summary["modes"][k]:.5f}' for k in
              ('kstar_charged_K+_pi0', 'kstar_charged_K0_pi+', 'kstar_neutral_K0_pi0', 'kstar_neutral_K+_pi-')))

    phi_preferred = phi_ratio_comparison(0.)['preferred_change']
    print(f'  LO phi K+K-/K0K0bar width ratio (phase space only) {check_phi_leading_ratio():.6f}, equal to the endpoint ratio')
    print('Step 4: two-coefficient solutions with the phi rate ratio at its preferred change')
    r2_alone, _ = solve_single_coefficient('r2', target)
    for free_names, start in ((('r2', 'r3'), (r2_alone, r2_alone)), (('r1', 'r2'), (-0.7, 1.0))):
        solution = solve_combined(free_names, target, phi_preferred, start)
        summary = scenario_summary(solution)
        print(f'  {free_names}: ' + ', '.join(f'{k} = {v:+.6f} GeV^-2 (r X_s = {v*xs:+.4f})' for k, v in solution.items()))
        print(f'    check R = {summary["kstar_width_ratio"]:.6f}, phi change {summary["phi_rate_ratio_change"]:+.4f}; '
              f'K*(K0pi+)/rho0 {summary["kstar_charged_over_rho_coupling"]:.4f}; '
              f'phi(K0K0bar)/rho0 {summary["phi_over_rho_coupling"]:.4f}; '
              f'omega->pipi/rho LO {summary["omega_pipi_over_rho_leading"]:+.5f}')
        print(f'    phi(K0K0bar)/K*(K0pi+) {summary["phi_over_kstar_coupling"]:.4f}; (rho0, omega -> KK)/K*(K0pi+): '
              + ', '.join(f'{k} {v:.4f}' for k, v in summary['light_vector_to_kaons_over_kstar'].items()))

    print('Step 5: electromagnetic single traces with two Q')
    rank, structures = electromagnetic_structures()
    print(f'  rank {rank}')
    linear_target = 0.5*math.log(target)
    for item in structures:
        needed = linear_target/item['kstar_split'] if abs(item['kstar_split']) > 1e-12 else float('inf')
        print(f'  <{item["word"]}>: K* split {item["kstar_split"]:+.4f}, phi split {item["phi_split"]:+.4f}, '
              f'rho+-rho0 {item["rho_split"]:+.4f}; required coefficient {needed:+.3f} '
              f'= {needed/(FINE_STRUCTURE_CONSTANT/math.pi):+.1f} alpha/pi')


if __name__ == '__main__':
    main()
