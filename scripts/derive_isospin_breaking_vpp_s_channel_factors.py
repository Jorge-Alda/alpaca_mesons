#!/usr/bin/env python3
"""Isospin-breaking VPP factors of the s-channel vector -> K Kbar vertices of the BaBar amplitude.

Invocation:
    .venv/bin/python scripts/derive_isospin_breaking_vpp_s_channel_factors.py

Derivation only: no fit, no refit, no change of the baseline likelihood. Units:
GeV^2 for the vacuum spurion X = chi_+ (Phi = 0) = 4 diag(B0 m_u, B0 m_d, B0 m_s);
GeV^-2 for r_i = sqrt2 lambda_i / G_V; the factors are dimensionless.
Conventions of Alda:2026nrd and of scripts/derive_isospin_breaking_vpp_operators.py:
flavor matrices K+ = E(0,2), K0 = E(1,2), K- = E(2,0), K0bar = E(2,1); a
decaying rho+ is represented by E(1,0) and a rho- by E(0,1), as the rho+ -> pi+ pi0
mode of fit_vector_gv.flavor_factors; rho0 = diag(1,-1,0)/sqrt2; light basis
vector L = diag(1,1,0)/sqrt2 and strange S = diag(0,0,1) as in
fit_babar_kkpi_global_vector_sector._babar_blocks_full; physical
omega = cos(theta) L - sin(theta) S, phi = sin(theta) L + cos(theta) S.

Vertex of each mode (operator record, Step 2):
    T(r) = tr(V[P1,P2]) + r1 tr(V{[P1,P2],X}) + r2 tr(V(P1 X P2 - P2 X P1))
           + r3 tr(V) tr([P1,P2] X);
each VPP coupling is multiplied by the mode factor T(r)/T(0). T is linear in V,
so the factor of a light/strange combination is fixed by the light and
strange traces separately.

Where the s-channel vertices sit in the baseline BaBar amplitude:
  * K+K-pi0: s = m^2(K+K-). `rs` is the rho0(s) propagator with the
    rho0 -> K+K- vertex; `lls = ll - sqrt2 ls` is the neutral propagator
    produced on its light end (the pi0 production vertices couple to L only
    in the baseline) and closed on the K+K- vertex, written in the light/strange
    basis: the ratio -sqrt2 is T_S(0)/T_L(0) for K+K- (checked below).
  * K_S K+- pi-+: s = m^2(K_S K+-) is a charged pair, so `rs` is the charged
    rho propagator (isospin-limit width of the baseline) with the
    rho+ -> K+ K0bar (or C-conjugate rho- -> K- K0) vertex; no omega, phi.

Steps (all asserted):
  1. Symbolic factors T(r)/T(0) of rho0 -> K+K-, rho0 -> K0K0bar,
     rho+ -> K+K0bar, rho- -> K-K0, L -> K+K-, S -> K+K-, omega -> K+K-,
     phi -> K+K-.
  2. At r = 0 every factor is exactly 1.
  3. T_S(0)/T_L(0) = -sqrt2 for K+K-, i.e. the baseline lls = ll - sqrt2 ls.
  4. Physical-basis and light/strange-basis assemblies of the s-channel
     neutral term agree for all theta, r:
        sqrt2 sum_n R_{L n} T_n(r) D_n = f_L ll - sqrt2 f_S ls,
     with R_{L n} = tr(V_n L), D_n the omega, phi pole propagators,
     ll = sum_n R_{L n}^2 D_n, ls = sum_n R_{L n} R_{S n} D_n; at r = 0 this
     is the baseline lls. Also sqrt2 T_n(0) = R_{L n} - sqrt2 R_{S n}.
  5. Charge structure: rho+ and rho- factors are equal (C), the charged-rho
     factor is the mean of the rho0 -> K+K- and rho0 -> K0K0bar factors, and
     in the isospin limit X_u = X_d all three rho factors coincide.
  6. Numerical: fitted omega-phi angle from the fitted generators
     sector.generators[1:] at the accepted endpoint basin_probes/even_shared_widths
     (the same matrices the likelihood option uses for the phi rows), and the
     factors at the scan point r = (-0.350, 1.855, 0.729) GeV^-2.

Output summary (2026-11-09 run):
  f(rho0 -> K+K-)    = 1 + 2 X_u r1 + X_s r2
  f(rho0 -> K0K0bar) = 1 + 2 X_d r1 + X_s r2
  f(rho+- -> K Kbar) = 1 + (X_u + X_d) r1 + X_s r2
  f_L(K+K-)          = 1 + 2 X_u r1 + X_s r2 + 2 (X_u - X_s) r3
  f_S(K+K-)          = 1 + 2 X_s r1 + X_u r2 + (X_s - X_u) r3
  fitted angle 0.0952440 rad (fixed provisional angle of the scan: 0.0824135).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_isospin_breaking_vpp_operators as operators  # noqa: E402
import fit_vector_gv as gv_fit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
E = gv_fit.matrix_unit
R1, R2, R3 = (operators.SYMBOLS[k] for k in ('r1', 'r2', 'r3'))
XU, XD, XS = (operators.SYMBOLS[k] for k in ('Xu', 'Xd', 'Xs'))
THETA = operators.SYMBOLS['theta']
COEFFICIENTS = (R1, R2, R3)
SPURIONS = (XU, XD, XS)
LIGHT = sp.diag(1, 1, 0)/sp.sqrt(2)
STRANGE = sp.diag(0, 0, 1)
SCAN_POINT_GEV_MINUS2 = (-.350, 1.855, .729)
ENDPOINT_ID = 'basin_probes/even_shared_widths'


def omega_phi_matrices(angle=THETA):
    """Physical omega and phi flavor matrices at the omega-phi angle (radians)."""
    return (sp.cos(angle)*LIGHT - sp.sin(angle)*STRANGE,
            sp.sin(angle)*LIGHT + sp.cos(angle)*STRANGE)


def s_channel_vertex_matrices(angle=THETA):
    """{mode: (V, P1, P2)} of the s-channel vector -> K Kbar vertices."""
    omega, phi = omega_phi_matrices(angle)
    rho_zero = sp.diag(1, -1, 0)/sp.sqrt(2)
    return {
        'rho0_K+_K-': (rho_zero, E(0, 2), E(2, 0)),
        'rho0_K0_K0bar': (rho_zero, E(1, 2), E(2, 1)),
        'rho+_K+_K0bar': (E(1, 0), E(0, 2), E(2, 1)),
        'rho-_K-_K0': (E(0, 1), E(2, 0), E(1, 2)),
        'light_K+_K-': (LIGHT, E(0, 2), E(2, 0)),
        'strange_K+_K-': (STRANGE, E(0, 2), E(2, 0)),
        'omega_K+_K-': (omega, E(0, 2), E(2, 0)),
        'phi_K+_K-': (phi, E(0, 2), E(2, 0))}


def effective_trace(vector, first, second):
    """T(r) of one mode with symbolic X and r_i (structure traces of the operator script)."""
    parts = operators.structure_traces(vector, first, second, sp.diag(*SPURIONS))
    return sp.expand(parts['leading'] + R1*parts['anticommutator'] + R2*parts['sandwich']
                     + R3*parts['vector_singlet'])


def at_zero(expression):
    """Expression at r1 = r2 = r3 = 0."""
    return expression.subs(dict.fromkeys(COEFFICIENTS, 0))


def symbolic_traces(angle=THETA):
    """{mode: T(r)} for the s-channel modes."""
    return {name: effective_trace(*matrices) for name, matrices in s_channel_vertex_matrices(angle).items()}


def symbolic_factors(angle=THETA):
    """{mode: T(r)/T(0)}, simplified (dimensionless)."""
    return {name: sp.simplify(trace/at_zero(trace)) for name, trace in symbolic_traces(angle).items()}


EXPECTED_FACTORS = {
    'rho0_K+_K-': 1 + 2*XU*R1 + XS*R2,
    'rho0_K0_K0bar': 1 + 2*XD*R1 + XS*R2,
    'rho+_K+_K0bar': 1 + (XU + XD)*R1 + XS*R2,
    'rho-_K-_K0': 1 + (XU + XD)*R1 + XS*R2,
    'light_K+_K-': 1 + 2*XU*R1 + XS*R2 + 2*(XU - XS)*R3,
    'strange_K+_K-': 1 + 2*XS*R1 + XU*R2 + (XS - XU)*R3}


def check_closed_forms(factors):
    """The rho, light and strange factors equal the closed forms of EXPECTED_FACTORS."""
    for name, expected in EXPECTED_FACTORS.items():
        assert sp.simplify(factors[name] - expected) == 0, name
    return True


def check_zero_coefficients(factors):
    """Every factor is exactly 1 at r = 0."""
    for name, factor in factors.items():
        assert sp.simplify(at_zero(factor) - 1) == 0, name
    return True


def check_baseline_light_strange_ratio(traces):
    """T_S(0)/T_L(0) = -sqrt2 for K+K-, the coefficient of ls in the baseline lls."""
    ratio = sp.simplify(at_zero(traces['strange_K+_K-'])/at_zero(traces['light_K+_K-']))
    assert ratio == -sp.sqrt(2)
    return ratio


def check_physical_equals_light_strange(traces, factors):
    """Physical-basis and light/strange assemblies of the s-channel omega/phi term.

    D_omega, D_phi are symbolic pole propagators. Returns the r = 0 limit,
    which must be the baseline ll - sqrt2 ls.
    """
    d_omega, d_phi = sp.symbols('D_omega D_phi')
    omega, phi = omega_phi_matrices()
    rotation = sp.Matrix([[sp.trace(v*basis) for v in (omega, phi)] for basis in (LIGHT, STRANGE)])
    assert sp.simplify(rotation - sp.Matrix([[sp.cos(THETA), sp.sin(THETA)],
                                             [-sp.sin(THETA), sp.cos(THETA)]])) == sp.zeros(2)
    propagators = (d_omega, d_phi)
    light_light = sum(rotation[0, n]**2*propagators[n] for n in range(2))
    light_strange = sum(rotation[0, n]*rotation[1, n]*propagators[n] for n in range(2))
    physical = sp.sqrt(2)*sum(rotation[0, n]*traces[name]*propagators[n]
                              for n, name in enumerate(('omega_K+_K-', 'phi_K+_K-')))
    light_strange_form = (factors['light_K+_K-']*light_light
                          - sp.sqrt(2)*factors['strange_K+_K-']*light_strange)
    assert sp.simplify(sp.expand(physical - light_strange_form)) == 0
    baseline = light_light - sp.sqrt(2)*light_strange
    assert sp.simplify(at_zero(physical) - baseline) == 0
    for n, name in enumerate(('omega_K+_K-', 'phi_K+_K-')):
        assert sp.simplify(sp.sqrt(2)*at_zero(traces[name])
                           - (rotation[0, n] - sp.sqrt(2)*rotation[1, n])) == 0
    return baseline


def check_charge_structure(factors):
    """C equality, charged = neutral mean, isospin-limit equality of the rho factors."""
    assert sp.simplify(factors['rho+_K+_K0bar'] - factors['rho-_K-_K0']) == 0
    mean = (factors['rho0_K+_K-'] + factors['rho0_K0_K0bar'])/2
    assert sp.simplify(factors['rho+_K+_K0bar'] - mean) == 0
    limit = {XD: XU}
    assert sp.simplify((factors['rho+_K+_K0bar'] - factors['rho0_K+_K-']).subs(limit)) == 0
    assert sp.simplify((factors['rho0_K0_K0bar'] - factors['rho0_K+_K-']).subs(limit)) == 0
    return True


def fitted_mixing_angle(values=None):
    """omega-phi angle (rad) of the fitted generators at the accepted endpoint.

    Reads sector.generators[1], [2] from GlobalLikelihood.unpack (reduced
    orders 1, 3, 4: the generators do not depend on quadrature) and asserts
    that both lie in the light/strange plane with omega = cos L - sin S,
    phi = sin L + cos S to 1e-14.
    """
    import fit_babar_kkpi_global_vector_sector as joint_fit
    if values is None:
        comparison = json.loads((ROOT/'data/vector_sector/joint_vector_endpoint_comparison.json').read_text())
        candidate = next(c for c in comparison['candidates'] if c['id'] == ENDPOINT_ID)
        assert candidate['accepted']
        values = np.asarray(candidate['endpoint'], dtype=float)
    model = joint_fit.GlobalLikelihood(1, 3, 4)
    generators = model.unpack(values)[1].generators
    phi = np.real(np.asarray(generators[2]))
    angle = math.atan2(math.sqrt(2)*phi[0, 0], phi[2, 2])
    omega_expected, phi_expected = (np.array(sp.N(m), dtype=float) for m in omega_phi_matrices(angle))
    assert np.max(abs(np.asarray(generators[1]) - omega_expected)) < 1e-14
    assert np.max(abs(np.asarray(generators[2]) - phi_expected)) < 1e-14
    return angle, model.preliminary['omega_phi_mixing_angle_rad']


def numerical_factors(coefficients, angle, spurion_values=None):
    """Factors of all s-channel modes at numerical r (GeV^-2), angle (rad) and X (GeV^2)."""
    spurion_values = spurion_values or operators.vacuum_spurion_values()
    numbers = dict(zip(SPURIONS, spurion_values)) | dict(zip(COEFFICIENTS, coefficients))
    numbers[THETA] = angle
    return {name: float(sp.N(factor.subs(numbers), 30)) for name, factor in symbolic_factors().items()}


def main():
    traces = symbolic_traces()
    factors = symbolic_factors()
    print('Step 1: s-channel vertex factors T(r)/T(0)')
    for name, factor in factors.items():
        print(f'  {name:16s} {sp.simplify(factor)}')
    print('  closed forms (rho, light, strange):', check_closed_forms(factors))
    print('Step 2: all factors equal 1 at r = 0:', check_zero_coefficients(factors))
    print('Step 3: T_S(0)/T_L(0) for K+K- =', check_baseline_light_strange_ratio(traces))
    print('Step 4: physical basis = light/strange basis; r = 0 limit:',
          check_physical_equals_light_strange(traces, factors))
    print('Step 5: C equality, charged rho = neutral mean, isospin-limit equality:',
          check_charge_structure(factors))
    angle, provisional = fitted_mixing_angle()
    print(f'Step 6: fitted omega-phi angle {angle:.7f} rad from sector.generators at {ENDPOINT_ID} '
          f'(provisional fixed angle {provisional:.7f})')
    xu, xd, xs = operators.vacuum_spurion_values()
    print(f'  X = ({xu:.6f}, {xd:.6f}, {xs:.6f}) GeV^2; factors at r = {SCAN_POINT_GEV_MINUS2} GeV^-2:')
    for name, value in numerical_factors(SCAN_POINT_GEV_MINUS2, angle).items():
        print(f'  {name:16s} {value:.6f}')


if __name__ == '__main__':
    main()
