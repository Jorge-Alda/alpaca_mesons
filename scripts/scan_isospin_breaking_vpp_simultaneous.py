#!/usr/bin/env python3
"""Simultaneous variation of r_1, r_2, r_3 in the isospin-breaking VPP model.

Invocation:
    .venv/bin/python scripts/scan_isospin_breaking_vpp_simultaneous.py

Exploratory scan, no refit. Model, conventions and units are those of
scripts/derive_isospin_breaking_vpp_operators.py (imported, not modified):
effective flavor trace
    T = tr(V[P1,P2]) + r1 tr(V{[P1,P2],X}) + r2 tr(V(P1 X P2 - P2 X P1))
        + r3 tr(V) tr([P1,P2] X),
X = chi_+ (vacuum) = 4 B0 diag(m_u, m_d, m_s) in GeV^2, r_i in GeV^-2,
published pi0 diagonals of Alda:2026nrd, nonideal omega-phi angle of the
provisional fit, `width_coefficient` phase space with PDG parent masses.

Targets (from memory/rcht/odd-parity/joint-vector-isospin-breaking-vpp-operators.md
and joint-vector-kstar-kpi-saturation-freedom.md):
  (a) R = [Gamma_c(r)/Gamma_c(0)]/[Gamma_n(r)/Gamma_n(0)] = kappa_c/kappa_n
      = 1.094483 (width-ratio convention; coupling ratio sqrt(R) = 1.046175);
  (b) phi rate-ratio change [Gamma(K+K-)/Gamma(K0K0bar)](r)/(r=0) - 1
      = -4.43 % +- 2.0 % (preferred change and its measured relative error);
  (c) rho+/rho0 coupling ratio relative to r = 0 equal to 1.
Assumption for the chi-square: the error of (a) is the relative error of the
measured PDG width ratio Gamma(K*+-)/Gamma(K*0) (uncorrelated rows of
gv_measurements.json), 1.88 %, times the target; saturated widths carry no error.
Assumption for the exact solution: the phi K+K-/K0K0bar amplitude ratio keeps
its r = 0 sign (positive square-root branch).

Steps:
  1. SymPy: the linearised splits are (r1 - r2/3)(X_u - X_d) and
     (r2 - r3)(X_u - X_d) in the isospin limit with ideal mixing; with the code
     inputs (published pi0, fit angle, physical phase space) compare the
     numerical gradient at r = 0 with these directions. The exact conditions
     are nonlinear (ratios of quadratic forms) and are solved exactly.
  2. Exact one-parameter family: R(r1, r2) = target is a conic in (r1, r2)
     (r3 does not enter R); for each r1 solve the quadratic in r2, then the
     phi condition is linear in r3 and is solved in closed form.
  3. 3D grid r_i in [-5, 5] GeV^-2: chi-square of (a) and (b), couplings.
  4. Repeat with X_u - X_d and X_s varied by +-20 %.

Output summary (values of the 2026-11-09 run, about 10 s):
  * Isospin limit, ideal mixing: splits (r1 - r2/3)(X_u-X_d), (r2 - r3)(X_u-X_d),
    rho+/rho0 trace ratio exactly -1 for all r_i. Code inputs: gradients in
    units of (X_u-X_d) are K* (1, -0.663, 0), phi (-0.124, 1.062, -1.186);
    the r2 K* coefficient doubles through X_s times the s sbar part (-0.0114)
    of the published pi0. Reduction not exact -> exact solve used.
  * Exact family exists for every r1 in [-5, 5] (two conic branches); it
    contains the record points (0, 3.7362, 1.2482) and (-0.8415, -0.7871, 0).
  * Point of the family closest to unit common-mode couplings (Alda X):
    r = (-0.350, +1.855, +0.729) GeV^-2, g(K*)/g(rho) 0.695,
    g(phi)/g(K*) 1.444, g(phi)/g(rho) 1.003 (largest |log| 0.367).
  * Grid (8.1e6 points): min chi2 4e-6; chi2 < 1 point closest to unit
    couplings has largest |log| 0.242 (g(K*)/g(rho) 0.785, g(phi)/g(K*) 1.260).
    |rho+/rho0 - 1| <= 2.0e-3 on the grid, from the published pi0 only (it
    grows to about 2 % on the family where |r2| is about 30).
  * X_d-X_u and X_s by +-20 %: closest point keeps largest |log| 0.33-0.37.
"""
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_isospin_breaking_vpp_operators as derivation  # noqa: E402
import fit_vector_gv as gv_fit  # noqa: E402

R1, R2, R3 = (derivation.SYMBOLS[k] for k in ('r1', 'r2', 'r3'))
XU, XD, XS = (derivation.SYMBOLS[k] for k in ('Xu', 'Xd', 'Xs'))
COEFFICIENT_RANGE = 5.0  # GeV^-2
GRID_STEP = 0.05  # GeV^-2
CHANNELS_USED = ('kstar_charged_K+_pi0', 'kstar_charged_K0_pi+', 'kstar_neutral_K0_pi0',
                 'kstar_neutral_K+_pi-', 'phi_K+_K-', 'phi_K0_K0bar', 'rho0_pi+_pi-', 'rho+_pi+_pi0')


# ---------------------------------------------------------------------------
# Linear trace coefficients and phase space
# ---------------------------------------------------------------------------
def trace_coefficient_functions():
    """Per channel, numerical functions of (X_u, X_d, X_s) giving (T0, dT/dr1, dT/dr2, dT/dr3).

    T is exactly linear in r_i; T0 is the LO trace (independent of X).
    """
    mixing_angle = gv_fit.fixed_parameters()[2]
    published = derivation.alda_inputs()['published_lambda_diagonals']['pi0']
    traces = derivation.effective_traces([sp.Float(x, 30) for x in published], mixing_angle)
    zero = {R1: 0, R2: 0, R3: 0}
    result = {}
    for name in CHANNELS_USED:
        trace = sp.expand(traces[name])
        parts = [trace.subs(zero)] + [sp.diff(trace, symbol) for symbol in (R1, R2, R3)]
        result[name] = [sp.lambdify((XU, XD, XS), part, 'math') for part in parts]
    return result


TRACE_FUNCTIONS = trace_coefficient_functions()


def trace_coefficients(spurion_values):
    """Return {channel: numpy array (T0, A1, A2, A3)} at the given X (GeV^2)."""
    return {name: np.array([float(f(*spurion_values)) for f in functions])
            for name, functions in TRACE_FUNCTIONS.items()}


def phase_space_weights():
    """Width_coefficient per K* mode (unit trace, unit tensor mass) and the phi K+K-/K0K0bar LO ratio."""
    masses = derivation.alda_inputs()['physical_masses_gev']
    rows, _ = derivation.measurements()
    weights = {}
    for charge, modes in derivation.KSTAR_MODES.items():
        parent = rows[f'kstar_{charge}_Kpi']['physical_vector_mass_GeV']
        for channel, first, second in modes:
            weights[channel] = gv_fit.width_coefficient(parent, 1.0, masses[first], masses[second], 1.0)
    parent = rows['phi_K+_K-']['physical_vector_mass_GeV']
    weights['phi_K+_K-'] = gv_fit.width_coefficient(parent, 1.0, masses['K_charged'], masses['K_charged'], 1.0)
    weights['phi_K0_K0bar'] = gv_fit.width_coefficient(parent, 1.0, masses['K_neutral'], masses['K_neutral'], 1.0)
    return weights


PHASE_SPACE = phase_space_weights()


def targets():
    """Return (a) target R and its assumed error, (b) preferred phi change and its error."""
    rows, covariance = derivation.measurements()
    kstar_target = derivation.required_kstar_ratio()
    relative_error = math.sqrt(covariance[0, 0]/rows['kstar_charged_Kpi']['value']**2
                               + covariance[1, 1]/rows['kstar_neutral_Kpi']['value']**2)
    phi = derivation.phi_ratio_comparison(0.0)
    return dict(kstar=kstar_target, kstar_error=kstar_target*relative_error,
                kstar_relative_error=relative_error,
                phi=phi['preferred_change'], phi_error=phi['relative_sigma'])


# ---------------------------------------------------------------------------
# Observables (vectorised in r_i)
# ---------------------------------------------------------------------------
def channel_traces(coefficients, r1, r2, r3):
    """Effective trace T(r) of every channel (arrays broadcast over r_i)."""
    return {name: c[0] + c[1]*r1 + c[2]*r2 + c[3]*r3 for name, c in coefficients.items()}


def observables(coefficients, r1, r2, r3):
    """Charge splits and common-mode couplings at r_i (all relative to r = 0).

    Couplings are signed and relative to their LO values:
      * kstar_mean: K* mode factors T(r)/T(0) averaged with the LO partial
        widths of all four K pi modes as weights;
      * phi_mean: same for phi -> K+K-, K0K0bar;
      * rho: rho0 -> pi+pi- factor.
    The record definitions (K*+ -> K0 pi+, phi -> K0K0bar) are also returned.
    """
    traces = channel_traces(coefficients, r1, r2, r3)
    leading = {name: c[0] for name, c in coefficients.items()}
    factor = {name: traces[name]/leading[name] for name in traces}

    def width(modes):
        return sum(PHASE_SPACE[m]*traces[m]**2 for m in modes)

    def leading_width(modes):
        return sum(PHASE_SPACE[m]*leading[m]**2 for m in modes)
    charged = [m for m, _, _ in derivation.KSTAR_MODES['charged']]
    neutral = [m for m, _, _ in derivation.KSTAR_MODES['neutral']]
    kstar_ratio = (width(charged)/leading_width(charged))/(width(neutral)/leading_width(neutral))
    phi_change = ((traces['phi_K+_K-']/traces['phi_K0_K0bar'])**2
                  / (leading['phi_K+_K-']/leading['phi_K0_K0bar'])**2 - 1)
    rho_split = (traces['rho+_pi+_pi0']/traces['rho0_pi+_pi-'])/(leading['rho+_pi+_pi0']/leading['rho0_pi+_pi-'])

    def weighted_mean(modes):
        weights = {m: PHASE_SPACE[m]*leading[m]**2 for m in modes}
        return sum(w*factor[m] for m, w in weights.items())/sum(weights.values())
    kstar_mean = weighted_mean(charged + neutral)
    phi_mean = weighted_mean(['phi_K+_K-', 'phi_K0_K0bar'])
    rho = factor['rho0_pi+_pi-']
    return dict(
        kstar_ratio=kstar_ratio, phi_change=phi_change, rho_split=rho_split,
        kstar_mean=kstar_mean, phi_mean=phi_mean, rho=rho,
        kstar_over_rho=kstar_mean/rho, phi_over_kstar=phi_mean/kstar_mean, phi_over_rho=phi_mean/rho,
        record_kstar_over_rho=factor['kstar_charged_K0_pi+']/rho,
        record_phi_over_kstar=factor['phi_K0_K0bar']/factor['kstar_charged_K0_pi+'],
        record_phi_over_rho=factor['phi_K0_K0bar']/rho)


def chi_square(values, target):
    """chi^2 of (a) and (b) with the errors of `targets` (assumption for (a))."""
    return (((values['kstar_ratio'] - target['kstar'])/target['kstar_error'])**2
            + ((values['phi_change'] - target['phi'])/target['phi_error'])**2)


# ---------------------------------------------------------------------------
# Step 1: linear reduction
# ---------------------------------------------------------------------------
def symbolic_linear_reduction():
    """SymPy check, isospin-limit pi0 and ideal mixing, equal-phase-space weights.

    Returns the linear K* and phi splits; asserts they equal
    (r1 - r2/3)(X_u - X_d) and (r2 - r3)(X_u - X_d), and that the rho+/rho0
    trace ratio is -1 for all r_i (exact, not only linear).
    """
    splits = derivation.linear_isospin_limit_splits()
    assert sp.simplify(splits['kstar'] - (R1 - R2/3)*(XU - XD)) == 0
    assert sp.simplify(splits['phi'] - (R2 - R3)*(XU - XD)) == 0
    traces = derivation.effective_traces((1, -1, 0), 0)
    assert sp.simplify(traces['rho+_pi+_pi0'] + traces['rho0_pi+_pi-']) == 0
    return splits


def numerical_linear_gradients(spurion_values, step=1e-6):
    """Gradients at r = 0 of (1/2) ln R and (1/2) ln(1 + phi change), code inputs.

    Returned in units of (X_u - X_d), to compare with (1, -1/3, 0) and (0, 1, -1).
    """
    coefficients = trace_coefficients(spurion_values)
    difference = spurion_values[0] - spurion_values[1]
    gradients = {'kstar': [], 'phi': []}
    for index in range(3):
        shift = np.zeros(3)
        shift[index] = step
        plus, minus = observables(coefficients, *shift), observables(coefficients, *(-shift))
        gradients['kstar'].append(0.25*(math.log(plus['kstar_ratio']) - math.log(minus['kstar_ratio']))/step/difference)
        gradients['phi'].append(0.25*(math.log(1 + plus['phi_change']) - math.log(1 + minus['phi_change']))/step/difference)
    return {k: np.array(v) for k, v in gradients.items()}


# ---------------------------------------------------------------------------
# Step 2: exact one-parameter family
# ---------------------------------------------------------------------------
def kstar_conic(coefficients, kstar_target):
    """Quadratic-form coefficients of Gamma_c(r)/Gamma_c0 - target Gamma_n(r)/Gamma_n0 in (r1, r2).

    Returns a function of r1 giving (a, b, c) of a r2^2 + b r2 + c = 0.
    """
    def normalised_quadratic(modes):
        leading_width = sum(PHASE_SPACE[m]*coefficients[m][0]**2 for m in modes)
        matrix = np.zeros((3, 3))  # basis (1, r1, r2)
        for m in modes:
            vector = coefficients[m][:3]
            matrix += PHASE_SPACE[m]*np.outer(vector, vector)
        return matrix/leading_width
    charged = [m for m, _, _ in derivation.KSTAR_MODES['charged']]
    neutral = [m for m, _, _ in derivation.KSTAR_MODES['neutral']]
    form = normalised_quadratic(charged) - kstar_target*normalised_quadratic(neutral)

    def quadratic_in_r2(r1):
        a = form[2, 2]
        b = 2*(form[0, 2] + form[1, 2]*r1)
        c = form[0, 0] + 2*form[0, 1]*r1 + form[1, 1]*r1**2
        return a, b, c
    return quadratic_in_r2


def solve_r3(coefficients, r1, r2, phi_target):
    """Closed-form r3 giving the phi change exactly (positive amplitude-ratio branch)."""
    charged, neutral = coefficients['phi_K+_K-'], coefficients['phi_K0_K0bar']
    ratio = charged[0]/neutral[0]*math.sqrt(1 + phi_target)
    charged_rest = charged[0] + charged[1]*r1 + charged[2]*r2
    neutral_rest = neutral[0] + neutral[1]*r1 + neutral[2]*r2
    return (ratio*neutral_rest - charged_rest)/(charged[3] - ratio*neutral[3])


def exact_family(spurion_values, target, r1_values):
    """Points (r1, r2, r3) satisfying (a) and (b) exactly, for both conic branches.

    Returns {branch sign: array of shape (n, 3)} with NaN where the discriminant is negative.
    """
    coefficients = trace_coefficients(spurion_values)
    quadratic = kstar_conic(coefficients, target['kstar'])
    family = {}
    for sign in (+1, -1):
        points = []
        for r1 in r1_values:
            a, b, c = quadratic(r1)
            discriminant = b*b - 4*a*c
            if discriminant < 0:
                points.append((r1, np.nan, np.nan))
                continue
            r2 = (-b + sign*math.sqrt(discriminant))/(2*a)
            points.append((r1, r2, solve_r3(coefficients, r1, r2, target['phi'])))
        family[sign] = np.array(points)
    return family


def family_table(spurion_values, target, r1_values):
    """Observables on every exact family point with finite r2."""
    coefficients = trace_coefficients(spurion_values)
    rows = []
    for sign, points in exact_family(spurion_values, target, r1_values).items():
        for r1, r2, r3 in points:
            if np.isnan(r2):
                continue
            values = observables(coefficients, r1, r2, r3)
            rows.append(dict(branch=sign, r1=r1, r2=r2, r3=r3, **values))
    return rows


def common_mode_distance(row):
    """Largest |log| of the three common-mode coupling ratios (inf for a sign flip)."""
    ratios = (row['kstar_over_rho'], row['phi_over_kstar'], row['phi_over_rho'])
    return max(abs(math.log(x)) if x > 0 else math.inf for x in ratios)


# ---------------------------------------------------------------------------
# Step 3: 3D grid
# ---------------------------------------------------------------------------
def grid_scan(spurion_values, target, limit=COEFFICIENT_RANGE, step=GRID_STEP):
    """chi^2 and couplings on a cubic grid r_i in [-limit, limit]."""
    axis = np.round(np.arange(-limit, limit + step/2, step), 10)
    r1, r2, r3 = np.meshgrid(axis, axis, axis, indexing='ij')
    values = observables(trace_coefficients(spurion_values), r1, r2, r3)
    chi2 = chi_square(values, target)
    best = np.unravel_index(np.argmin(chi2), chi2.shape)
    inside = chi2 < 1
    summary = dict(points=chi2.size, minimum=float(chi2[best]),
                   best=(float(r1[best]), float(r2[best]), float(r3[best])),
                   count_below_one=int(inside.sum()),
                   max_rho_split_deviation=float(np.max(np.abs(values['rho_split'] - 1))))
    for key in ('kstar_over_rho', 'phi_over_kstar', 'phi_over_rho'):
        selected = values[key][inside]
        summary[key] = (float(selected.min()), float(selected.max())) if selected.size else None
    # point with chi^2 < 1 closest to unit common-mode couplings
    if inside.any():
        distance = np.max(np.abs(np.log(np.abs(np.stack(
            [values[k] for k in ('kstar_over_rho', 'phi_over_kstar', 'phi_over_rho')])))), axis=0)
        sign_ok = (values['kstar_over_rho'] > 0) & (values['phi_over_kstar'] > 0)
        distance = np.where(inside & sign_ok, distance, np.inf)
        closest = np.unravel_index(np.argmin(distance), distance.shape)
        summary['closest_to_unit_couplings'] = dict(
            r=(float(r1[closest]), float(r2[closest]), float(r3[closest])),
            chi2=float(chi2[closest]), max_abs_log=float(distance[closest]),
            **{k: float(values[k][closest]) for k in ('kstar_over_rho', 'phi_over_kstar', 'phi_over_rho')})
    return summary


# ---------------------------------------------------------------------------
def print_family_summary(spurion_values, target, label):
    """Print the exact family summary for one set of X values."""
    r1_values = np.linspace(-COEFFICIENT_RANGE, COEFFICIENT_RANGE, 4001)
    rows = family_table(spurion_values, target, r1_values)
    print(f'  [{label}] X = ({spurion_values[0]:.6f}, {spurion_values[1]:.6f}, {spurion_values[2]:.6f}) GeV^2; '
          f'{len(rows)} exact points with r1 in [-5, 5]')
    within = [r for r in rows if max(abs(r['r2']), abs(r['r3'])) <= COEFFICIENT_RANGE]
    for branch in (+1, -1):
        selected = [r for r in rows if r['branch'] == branch]
        if not selected:
            print(f'    branch {branch:+d}: no real solutions')
            continue
        print(f'    branch {branch:+d}: r2 in [{min(r["r2"] for r in selected):+.3f}, {max(r["r2"] for r in selected):+.3f}], '
              f'r3 in [{min(r["r3"] for r in selected):+.3f}, {max(r["r3"] for r in selected):+.3f}]')
    best = min(within, key=common_mode_distance) if within else None
    if best is None:
        print('    no exact point with all |r_i| <= 5')
        return None
    for key in ('kstar_over_rho', 'phi_over_kstar', 'phi_over_rho', 'rho_split'):
        values = [r[key] for r in within]
        print(f'    {key} over exact points with |r_i|<=5: [{min(values):+.4f}, {max(values):+.4f}]')
    print(f'    exact point closest to unit common-mode couplings (branch {best["branch"]:+d}): '
          f'r = ({best["r1"]:+.4f}, {best["r2"]:+.4f}, {best["r3"]:+.4f}) GeV^-2, '
          f'r X_s = ({best["r1"]*spurion_values[2]:+.3f}, {best["r2"]*spurion_values[2]:+.3f}, {best["r3"]*spurion_values[2]:+.3f})')
    print(f'      R = {best["kstar_ratio"]:.6f}, phi change {best["phi_change"]:+.5f}, rho+/rho0 {best["rho_split"]:.6f}; '
          f'g(K*)/g(rho) {best["kstar_over_rho"]:.4f}, g(phi)/g(K*) {best["phi_over_kstar"]:.4f}, '
          f'g(phi)/g(rho) {best["phi_over_rho"]:.4f}; K* mean {best["kstar_mean"]:.4f}, '
          f'phi mean {best["phi_mean"]:.4f}, rho {best["rho"]:.4f}; max|log| {common_mode_distance(best):.3f}')
    print(f'      record definitions: g(K*0pi+)/g(rho) {best["record_kstar_over_rho"]:.4f}, '
          f'g(phi K0K0bar)/g(K*) {best["record_phi_over_kstar"]:.4f}, g(phi K0K0bar)/g(rho) {best["record_phi_over_rho"]:.4f}')
    return best


def main():
    target = targets()
    alda = derivation.vacuum_spurion_values()
    print(f'Targets: (a) R = {target["kstar"]:.6f} +- {target["kstar_error"]:.4f} '
          f'(assumed: relative PDG width-ratio error {target["kstar_relative_error"]:.4f}); '
          f'(b) phi change {target["phi"]:+.5f} +- {target["phi_error"]:.5f}; (c) rho+/rho0 = 1')

    print('Step 1: linear reduction')
    splits = symbolic_linear_reduction()
    print(f'  isospin limit, ideal mixing (SymPy, exact): K* split {splits["kstar"]}, phi split {splits["phi"]}, '
          f'rho+/rho0 trace ratio -1 for all r_i')
    gradients = numerical_linear_gradients(alda)
    print(f'  code inputs, gradient at r=0 in units of (X_u-X_d): K* {np.round(gradients["kstar"], 5)} '
          f'(isospin limit: [1, -0.33333, 0]); phi {np.round(gradients["phi"], 5)} (isospin limit: [0, 1, -1])')
    print('  exact conditions are ratios of quadratic forms: not functions of (r1-r2/3), (r2-r3) alone '
          '-> solved exactly (conic in r1, r2; linear in r3)')

    coefficients = trace_coefficients(alda)
    for name, point in (('record r1 alone', (-0.695084, 0, 0)), ('record r2 alone', (0, 3.736157, 0)),
                        ('record (r2, r3)', (0, 3.7362, 1.2482)), ('record (r1, r2)', (-0.8415, -0.7871, 0))):
        v = observables(coefficients, *point)
        print(f'  check {name} {point}: R {v["kstar_ratio"]:.6f}, phi change {v["phi_change"]:+.4f}, '
              f'record g(K*)/g(rho) {v["record_kstar_over_rho"]:.4f}, g(phi)/g(K*) {v["record_phi_over_kstar"]:.4f}, '
              f'g(phi)/g(rho) {v["record_phi_over_rho"]:.4f}; mean-coupling g(K*)/g(rho) {v["kstar_over_rho"]:.4f}, '
              f'g(phi)/g(K*) {v["phi_over_kstar"]:.4f}')

    print('Step 2: exact family satisfying (a) and (b)')
    print_family_summary(alda, target, 'Alda X')

    print(f'Step 3: 3D grid, step {GRID_STEP} GeV^-2, |r_i| <= {COEFFICIENT_RANGE}')
    summary = grid_scan(alda, target)
    print(f'  {summary["points"]} points; min chi2 {summary["minimum"]:.3e} at r = {summary["best"]}; '
          f'{summary["count_below_one"]} points with chi2 < 1; max |rho+/rho0 - 1| on grid {summary["max_rho_split_deviation"]:.2e}')
    for key in ('kstar_over_rho', 'phi_over_kstar', 'phi_over_rho'):
        print(f'  {key} range over chi2<1: {summary[key]}')
    print(f'  chi2<1 point closest to unit couplings: {summary.get("closest_to_unit_couplings")}')

    print('Step 4: X variations (X_u - X_d and X_s by +-20 %, X_u + X_d fixed)')
    mean, half_difference, strange = (alda[0] + alda[1])/2, (alda[1] - alda[0])/2, alda[2]
    for difference_scale in (0.8, 1.0, 1.2):
        for strange_scale in (0.8, 1.0, 1.2):
            values = (mean - difference_scale*half_difference, mean + difference_scale*half_difference,
                      strange_scale*strange)
            print_family_summary(values, target, f'X_d-X_u x{difference_scale}, X_s x{strange_scale}')


if __name__ == '__main__':
    main()
