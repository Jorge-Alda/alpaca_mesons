#!/usr/bin/env python3
"""Free the K pi saturation of the K*(892) couplings in the joint vector likelihood.

Invocation:
    .venv/bin/python scripts/derive_kstar_kpi_saturation_freedom.py

Derivation only: no refit and no change of the baseline likelihood
(scripts/fit_babar_kkpi_global_vector_sector.py). Units: GeV for masses,
widths, G_V and F; GeV^2 for invariants. Conventions as in Alda:2026nrd and
the baseline (F = 0.093 GeV, VPP vertex of fit_vector_gv.width_coefficient).

Definitions (notation proposed, awaiting user approval):
    kappa_q = (g_q / g_q^sat)^2,  q in {charged K*, neutral K*},
where g_q is the K*(892)_q -> K pi coupling used in the likelihood and g_q^sat
the coupling implied by the common G_V of the baseline. kappa_q = 1 is the
baseline. Steps:
  1. SymPy: Gamma_q = kappa_q Gamma_q^sat exactly; vertex scales by sqrt(kappa_q).
  2. Numerical: Gamma_q^sat from the SymPy formula reproduces the baseline
     `totals` at the accepted endpoint basin_probes/even_shared_widths.
  3. Numerical: the K* propagators kt, ku enter the BaBar amplitude only in
     G_V-proportional terms (direct, exchange); scaling each K* propagator by
     sqrt(kappa_q) of its charge state equals G_V -> G_V sqrt(kappa_q) on
     exactly those terms; kappa = 1 leaves the amplitude unchanged.
  4. Consequence: changing the width without the vertex is a different model.
Only reduced quadrature orders (1, 3, 4) without energy-bin blocks are used.
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_vector_gv as gv_fit  # noqa: E402
import joint_vector_kstar_charge_options as kstar_options  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASELINE_RUNNING_PROPAGATOR = joint_fit.running_propagator
KSTAR_NAMES = ('Kstarplus', 'Kstar0')  # charged, neutral


# ---------------------------------------------------------------------------
# Step 1: symbolic width and vertex scaling
# ---------------------------------------------------------------------------
def symbolic_width_scaling():
    """Return SymPy expressions and assert Gamma_q = kappa_q Gamma_q^sat."""
    G_V, F, physical_mass, tensor_mass, momentum, trace_squared, kappa = sp.symbols(
        'G_V F M_phys M_T p T2 kappa', positive=True)
    saturated_coupling = G_V*physical_mass**2*sp.sqrt(trace_squared)/(sp.sqrt(2)*F**2*tensor_mass)
    def width(coupling):
        return coupling**2*momentum**3/(6*sp.pi*physical_mass**2)
    saturated_width = width(saturated_coupling)
    code_formula = trace_squared*physical_mass**2*momentum**3/(12*sp.pi*F**4*tensor_mass**2)
    assert sp.simplify(saturated_width/G_V**2-code_formula) == 0
    free_coupling = sp.sqrt(kappa)*saturated_coupling
    assert sp.simplify(width(free_coupling)-kappa*saturated_width) == 0
    # The width fixes only |g_q|: both roots of g^2 = kappa (g^sat)^2.
    vertex_ratio = sp.symbols('r', real=True)
    roots = sp.solve(sp.Eq(width(vertex_ratio*saturated_coupling), kappa*saturated_width), vertex_ratio)
    assert set(roots) == {sp.sqrt(kappa), -sp.sqrt(kappa)}
    return dict(saturated_width=saturated_width, code_formula=code_formula, vertex_roots=roots)


def symbolic_flavor_traces():
    """Exact squared traces tr(V[P1,P2])^2 for the four K pi modes.

    Two pi0 choices: the isospin-limit diag(1,-1,0)/sqrt(2) and the
    Alda:2026nrd published diagonals used by fit_vector_gv.flavor_factors.
    """
    unit = gv_fit.matrix_unit
    kaon_plus, kaon_zero = unit(0, 2), unit(1, 2)
    pion_plus, pion_minus = unit(0, 1), unit(1, 0)
    kstar_minus, antikstar_zero = unit(2, 0), unit(2, 1)
    alda = json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())
    pion_zero_options = dict(
        isospin_limit=sp.diag(1, -1, 0)/sp.sqrt(2),
        alda_published=sp.diag(*[sp.nsimplify(x) for x in alda['published_lambda_diagonals']['pi0']])/sp.sqrt(2))
    result = {}
    for label, pion_zero in pion_zero_options.items():
        def squared(vector, first, second):
            return sp.nsimplify(sp.expand(sp.trace(vector*(first*second-second*first))**2))
        result[label] = {
            'kstar_charged_K+_pi0': squared(kstar_minus, kaon_plus, pion_zero),
            'kstar_charged_K0_pi+': squared(kstar_minus, kaon_zero, pion_plus),
            'kstar_neutral_K0_pi0': squared(antikstar_zero, kaon_zero, pion_zero),
            'kstar_neutral_K+_pi-': squared(antikstar_zero, kaon_plus, pion_minus)}
    return result


# ---------------------------------------------------------------------------
# Step 2: saturated widths at the accepted endpoint
# ---------------------------------------------------------------------------
def accepted_endpoint():
    """Coordinates of the accepted candidate basin_probes/even_shared_widths."""
    comparison = json.loads((ROOT/'data/vector_sector/joint_vector_endpoint_comparison.json').read_text())
    candidate = next(c for c in comparison['candidates'] if c['id'] == 'basin_probes/even_shared_widths')
    assert candidate['accepted']
    return np.asarray(candidate['endpoint'], dtype=float)


def saturated_widths(model, tensor_mass, coupling, symbolic):
    """Gamma_q^sat from the SymPy formula with the baseline masses and traces."""
    G_V, F, physical_mass, tensor, momentum, trace_squared = sp.symbols(
        'G_V F M_phys M_T p T2', positive=True)
    formula = sp.lambdify((G_V, F, physical_mass, tensor, momentum, trace_squared),
                          symbolic['saturated_width'])
    physical = model.alda['physical_masses_gev']
    widths = {}
    for row, name in zip(model.gv_data['rows'][:2], KSTAR_NAMES):
        parent = row['physical_vector_mass_GeV']
        total = 0.
        for first, second, flavor in kstar_options.KSTAR_MODES[name]:
            first_mass, second_mass = physical[first], physical[second]
            p = np.sqrt((parent**2-(first_mass+second_mass)**2)*(parent**2-(first_mass-second_mass)**2))/(2*parent)
            total += formula(coupling, gv_fit.DECAY_CONSTANT_GEV, parent, tensor_mass, p,
                             model.kaon_flavors[flavor])
        widths[name] = total
    return widths


# ---------------------------------------------------------------------------
# Step 3: K* charge-state substitution in the BaBar amplitude
# ---------------------------------------------------------------------------
def babar_amplitudes(model, values, state, kappa_charged=1., kappa_neutral=1.,
                     width_mode='baseline_mean', excited_scaling='ground_only',
                     kstar_widths=None, zero_kstar=False):
    """BaBar amplitudes with the K* propagators of each charge rescaled.

    width_mode: 'baseline_mean' uses the baseline mean of the two (scaled)
    widths in every K* propagator; 'charge_specific' uses the width of the
    propagator's own charge state (as the subclass option
    separate_kstar_charge_states). kstar_widths overrides the (charged,
    neutral) on-shell widths, to hold them fixed while the vertex changes.
    excited_scaling: 'ground_only' multiplies only the K*(892) pole by
    sqrt(kappa_q); 'whole' multiplies the full propagator including the
    beta-weighted K*(1410), K*(1680) terms.
    zero_kstar replaces every K* propagator by 0.
    Returns {channel: complex amplitude array} and the residual vector.
    """
    point, sector, mk, gv, totals = state
    totals = dict(totals)
    totals['Kstarplus'] *= kappa_charged
    totals['Kstar0'] *= kappa_neutral
    if kstar_widths is not None:
        totals['Kstarplus'], totals['Kstar0'] = kstar_widths

    def rescaled(name, invariant, beta, mass, shared_width, reference):
        if name != 'kstar':
            return BASELINE_RUNNING_PROPAGATOR(name, invariant, beta, mass, shared_width, reference)
        if zero_kstar:
            return np.zeros_like(invariant, dtype=complex)
        neutral = invariant is model.phase['kskpi']['t']
        kappa = kappa_neutral if neutral else kappa_charged
        width = shared_width if width_mode == 'baseline_mean' else totals['Kstar0' if neutral else 'Kstarplus']
        full = BASELINE_RUNNING_PROPAGATOR(name, invariant, beta, mass, width, reference)
        if kappa == 1.:
            return full
        if excited_scaling == 'whole':
            return np.sqrt(kappa)*full
        ground = BASELINE_RUNNING_PROPAGATOR(name, invariant, (0., 0.), mass, width, reference)
        return np.sqrt(kappa)*ground+(full-ground)

    joint_fit.running_propagator = rescaled
    model.babar_cache_key = None
    try:
        residual, _ = model._babar_blocks_full(values, (point, sector, mk, gv, totals))
    finally:
        joint_fit.running_propagator = BASELINE_RUNNING_PROPAGATOR
        model.babar_cache_key = None
    return {c: model.current_amplitudes[c].copy() for c in joint_fit.CHANNELS}, residual


def max_relative_difference(first, second):
    """Largest |a-b| over all channels divided by the largest |b|."""
    return max(float(np.max(abs(first[c]-second[c]))/np.max(abs(second[c]))) for c in first)


def relative_size(amplitudes, reference):
    """Largest |amplitude| over all channels divided by the largest |reference|."""
    return max(float(np.max(abs(amplitudes[c]))/np.max(abs(reference[c]))) for c in amplitudes)


def main():
    report = {}
    print('Step 1 (SymPy): width and vertex scaling')
    symbolic = symbolic_width_scaling()
    print('  Gamma^sat =', symbolic['saturated_width'])
    print('  Gamma^sat/G_V^2 equals fit_vector_gv.width_coefficient formula:', True)
    print('  Gamma(sqrt(kappa) g^sat) = kappa Gamma^sat exactly:', True)
    print('  roots of Gamma(r g^sat) = kappa Gamma^sat:', symbolic['vertex_roots'],
          '(width fixes |g_q| only; +sqrt(kappa) keeps the baseline sign)')

    traces = symbolic_flavor_traces()
    for label, values in traces.items():
        charged = values['kstar_charged_K+_pi0']+values['kstar_charged_K0_pi+']
        neutral = values['kstar_neutral_K0_pi0']+values['kstar_neutral_K+_pi-']
        print(f'  squared traces ({label}):', {k: str(v) for k, v in values.items()},
              f'sum charged = {charged}, sum neutral = {neutral}')
    iso = traces['isospin_limit']
    assert iso['kstar_charged_K0_pi+'] == 2*iso['kstar_charged_K+_pi0'] == 1
    assert iso['kstar_neutral_K+_pi-'] == 2*iso['kstar_neutral_K0_pi0'] == 1

    print('\nStep 2: saturated widths at basin_probes/even_shared_widths (orders 1, 3, 4)')
    values = accepted_endpoint()
    model = joint_fit.GlobalLikelihood(1, 3, 4)
    for key, exact in traces['alda_published'].items():
        assert abs(float(exact)-model.kaon_flavors[key]) < 1e-15, key
    state = model.unpack(values)
    point, sector, mk, gv, totals = state
    widths = saturated_widths(model, mk, gv, symbolic)
    for name in KSTAR_NAMES:
        relative = abs(widths[name]/totals[name]-1)
        print(f'  {name}: sympy {widths[name]:.12g} GeV, baseline totals {totals[name]:.12g} GeV, '
              f'relative difference {relative:.2e}')
        assert relative < 1e-13
    pdg = {name: model.odd_data['total_widths'][name]['value_gev'] for name in KSTAR_NAMES}
    kappa_pdg = {name: pdg[name]/totals[name] for name in KSTAR_NAMES}
    print(f'  G_V = {gv:.10g} GeV, m_K* (tensor) = {mk:.10g} GeV')
    print('  PDG total widths:', pdg, ' kappa matching them:', {k: round(v, 6) for k, v in kappa_pdg.items()})

    print('\nStep 3: BaBar amplitude, term-by-term K* charge assignment')
    baseline_residual, _ = model._babar_blocks_full(values, state)
    baseline = {c: model.current_amplitudes[c].copy() for c in joint_fit.CHANNELS}
    unit, unit_residual = babar_amplitudes(model, values, state)
    difference = max_relative_difference(unit, baseline)
    print(f'  kappa = (1, 1), baseline mean width: max relative amplitude difference {difference:.1e}; '
          f'residuals bitwise equal: {np.array_equal(unit_residual, baseline_residual)}')
    assert difference == 0. and np.array_equal(unit_residual, baseline_residual)
    # kappa one unit in the last place above 1 passes through the ground/excited split
    forced_split, _ = babar_amplitudes(model, values, state, np.nextafter(1., 2.), np.nextafter(1., 2.))
    split_difference = max_relative_difference(forced_split, baseline)
    print(f'  kappa = 1 + 1 ulp through the ground/excited split: relative difference {split_difference:.1e}')
    assert split_difference < 1e-13

    # (a) G_V = 0: no K* propagator survives, so no G_V-independent term has a K*.
    zero_state = (point, sector, mk, 0., totals)
    with_kstar, _ = babar_amplitudes(model, values, zero_state)
    without_kstar, _ = babar_amplitudes(model, values, zero_state, zero_kstar=True)
    gv_independent = max_relative_difference(with_kstar, without_kstar)
    print(f'  G_V = 0: K* propagators on vs off, relative difference {gv_independent:.1e}')
    assert gv_independent == 0.

    # (b) With widths fixed, scaling every K* propagator by lambda equals
    # G_V -> lambda G_V on the K* terms only.
    fixed = (totals['Kstarplus'], totals['Kstar0'])
    scale = 1.37
    scaled_propagators, _ = babar_amplitudes(model, values, state, scale**2, scale**2,
                                             excited_scaling='whole', kstar_widths=fixed)
    scaled_coupling_state = (point, sector, mk, scale*gv, totals)
    full_scaled, _ = babar_amplitudes(model, values, scaled_coupling_state)
    rest_scaled, _ = babar_amplitudes(model, values, scaled_coupling_state, zero_kstar=True)
    rest, _ = babar_amplitudes(model, values, state, zero_kstar=True)
    expected = {c: full_scaled[c]-rest_scaled[c]+rest[c] for c in full_scaled}
    identity = max_relative_difference(scaled_propagators, expected)
    print(f'  K* propagators x {scale} == G_V x {scale} on K* terms only: relative difference {identity:.1e}')
    assert identity < 1e-13

    # (c) Charge assignment: with widths fixed the amplitude is affine in
    # (sqrt(kappa_c), sqrt(kappa_n)); kappa_n does not enter K+K-pi0.
    def at(x, y):
        return babar_amplitudes(model, values, state, x**2, y**2, kstar_widths=fixed)[0]
    a00, a10, a01, a11 = at(1., 1.), at(1.5, 1.), at(1., 1.5), at(1.5, 1.5)
    cross = relative_size({c: a11[c]-a10[c]-a01[c]+a00[c] for c in a00}, a00)
    linear = relative_size({c: at(2., 1.)[c]-2*a10[c]+a00[c] for c in a00}, a00)
    neutral_in_kkpi0 = float(np.max(abs(a01['kkpi0']-a00['kkpi0'])))
    neutral_in_kskpi = float(np.max(abs(a01['kskpi']-a00['kskpi']))/np.max(abs(a00['kskpi'])))
    charged_in_kskpi = float(np.max(abs(a10['kskpi']-a00['kskpi']))/np.max(abs(a00['kskpi'])))
    print(f'  no cross term in (sqrt kappa_c, sqrt kappa_n): {cross:.1e}; second difference in sqrt kappa_c: {linear:.1e}')
    print(f'  kappa_n effect on K+K-pi0: {neutral_in_kkpi0:.1e} (absolute); on KsK pi: {neutral_in_kskpi:.2e}; '
          f'kappa_c effect on KsK pi: {charged_in_kskpi:.2e} (relative)')
    assert cross < 1e-13 and linear < 1e-13 and neutral_in_kkpi0 == 0.
    assert neutral_in_kskpi > 1e-3 and charged_in_kskpi > 1e-3

    # (d) Agreement with the subclass charge-specific widths at kappa = 1.
    subclass = kstar_options.KstarChargeLikelihood(1, 3, 4, separate_kstar_charge_states=True,
                                                   free_kstar_mass_difference=True)
    extended = subclass.extend(values)
    extended[-1] = 0.
    subclass_residual = subclass.residual(extended)
    model_state = model.unpack(values)
    _, own_residual = babar_amplitudes(model, values, model_state, width_mode='charge_specific')
    babar_rows = slice(28, 28+len(own_residual))
    subclass_difference = float(np.max(abs(subclass_residual[babar_rows]-own_residual)))
    print(f'  charge-specific widths, kappa = 1, Delta m = 0: max |residual difference| to subclass {subclass_difference:.1e}')
    assert subclass_difference < 1e-10

    print('\nStep 4: width-only versus width-and-vertex at kappa_q = Gamma_q^PDG / Gamma_q^sat')
    kc, kn = kappa_pdg['Kstarplus'], kappa_pdg['Kstar0']
    for excited in ('ground_only', 'whole'):
        _, coupled = babar_amplitudes(model, values, state, kc, kn, width_mode='charge_specific',
                                      excited_scaling=excited)
        report[excited] = coupled
    width_only_amplitudes, width_only = babar_amplitudes(
        model, values, state, width_mode='charge_specific',
        kstar_widths=(kc*totals['Kstarplus'], kn*totals['Kstar0']))
    _, mean_coupled = babar_amplitudes(model, values, state, kc, kn)
    chi2 = {k: float(np.sum(v**2)) for k, v in dict(
        baseline=baseline_residual, width_only=width_only, mean_width_coupled=mean_coupled, **report).items()}
    print('  BaBar chi2 at reduced orders (endpoint values unchanged otherwise):',
          {k: round(v, 4) for k, v in chi2.items()})
    coupled_amplitudes, _ = babar_amplitudes(model, values, state, kc, kn, width_mode='charge_specific')
    print(f'  max relative amplitude difference width-only vs coupled (ground only): '
          f'{max_relative_difference(width_only_amplitudes, coupled_amplitudes):.3e}')


if __name__ == '__main__':
    main()
