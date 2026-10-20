"""Independent symbolic/legacy checks of the numerical bounded SM kernels.

Run .venv/bin/python scripts/verify_odd_bounded_amplitudes.py.
No fit; exact compact artifacts and existing action amplitudes are references.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp
from diagnose_vector_stage4 import module, BASIS, ZERO, E

ROOT = Path(__file__).resolve().parents[1]


def main():
    bounded = module('odd_bounded_amplitudes')
    sectors = module('vector_sector')
    old = module('vector_sector_amplitudes')
    p = module('alda_2026nrd_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree')
    embedding = bounded.original_embedding()
    assert embedding.shape == (241, 11) and np.linalg.matrix_rank(embedding) == 11
    source = np.diag([.11, -.07, .13])
    # One split-pole/source-sign case and two angles expose responses that a
    # universal/common-pole comparison alone would fail to distinguish.
    checks = 0
    for angle in (.073, -.9):
        generators = sectors.physical_generators(angle)
        residues = tuple(float(np.trace(t@source).real) for t in generators)
        sector = sectors.VectorSector((.771, .783, 1.019), generators, residues)
        for name in ('pi0', 'eta', 'etap'):
            x = old.VectorAmplitudeInputs(.093, .770, p[name]['mp'], p[name]['TP'], p[name]['chi'], sector)
            for s, t in ((0., 0.), (-.07, 0.), (-.11, -.03)):
                row = bounded.tff_gradients(x, s, t)
                reference = np.array([sum(old.tff_components(c, x, s, t))
                                      - sum(old.tff_components(ZERO, x, s, t)) for c in BASIS])
                np.testing.assert_allclose(row@embedding, reference, rtol=2e-12, atol=2e-13)
                np.testing.assert_allclose(row, bounded.tff_gradients(x, t, s), rtol=2e-12, atol=2e-13)
                np.testing.assert_allclose(bounded.wzw_offset(x), old.tff_components(ZERO, x, s, t)[0])
                checks += 1
            for i, TV in enumerate(generators):
                row = bounded.radiative_gradients(x, TV, sector.masses[i], .3028)
                reference = np.array([sum(old.radiative_components(c, x, TV, sector.masses[i], .3028))
                                      for c in BASIS]).real
                np.testing.assert_allclose(row@embedding, reference, rtol=2e-12, atol=2e-13)
                checks += 1
        for i, j, mv in ((0, 1, .775), (0, 2, .892), (1, 2, .896)):
            chi = p['pi0']['chi']
            mp = np.sqrt((chi[i, i]+chi[j, j])/2)
            x = old.VectorAmplitudeInputs(.093, .770, mp, E(i, j), chi, sector)
            row = bounded.radiative_gradients(x, E(j, i), mv, .3028)
            reference = np.array([sum(old.radiative_components(c, x, E(j, i), mv, .3028))
                                  for c in BASIS]).real
            np.testing.assert_allclose(row@embedding, reference, rtol=2e-12, atol=2e-13)
            checks += 1
    # Every charged raw column is independently compared to the earlier
    # symbolic artifact at exact-rational input values, not to another
    # implementation of the compact word-evaluation routine.
    saved = json.loads((ROOT/'data/vector_sector/odd_cubic_parent_kernels.json').read_text())
    chi = np.diag([.02, .03, .48])
    sector = sectors.VectorSector.universal(.8, .14)
    w0 = sector.response().real
    substitutions = dict(F=.093, M_V=.770, M_1=.9, e=.3028, s=.9**2,
                         chi_u=.02, chi_d=.03, chi_s=.48,
                         W_0u=w0[0, 0], W_0d=w0[1, 1], W_0s=w0[2, 2])
    raw_checks = 0
    for name, i, j in (('rho_charged_pi', 0, 1), ('Kstar_charged_K', 0, 2),
                       ('Kstar_neutral_K', 1, 2)):
        mp2 = (chi[i, i]+chi[j, j])/2
        x = old.VectorAmplitudeInputs(.093, .770, np.sqrt(mp2), E(i, j), chi, sector)
        row = bounded.radiative_gradients(x, E(j, i), .9, .3028)
        values = []
        for identifier in bounded.parent_ids():
            expr = sp.sympify(saved['charged_channels'][name]['raw_parent_gradients'][identifier])
            replace = {symbol: substitutions.get(str(symbol), mp2 if str(symbol) == 'm_a2' else symbol)
                       for symbol in expr.free_symbols}
            values.append(float(expr.subs(replace)))
        np.testing.assert_allclose(row, values, rtol=2e-12, atol=2e-13)
        raw_checks += len(values)
    for scheme in ('all_derivatives_soft', 'hard_vector_derivatives_promoted'):
        subsets = [set(bounded.retained_columns(scheme, level)) for level in (0, 1, 2)]
        assert subsets[0] <= subsets[1] <= subsets[2] == set(range(241))
    measurements = json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text())
    pi0 = next(item for item in measurements['observations'] if item['channel'] == 'pi0_gamma_gamma')
    pi0_error = sp.sqrt(sp.Rational(52, 1000)**2+sp.Rational(105, 1000)**2)/10**9
    assert pi0['observable'] == 'partial_width' and pi0['unit'] == 'GeV'
    np.testing.assert_allclose(pi0['sigma'], float(pi0_error), rtol=1e-12, atol=0.)
    result = dict(passed=True, original_embedding_rank=11, legacy_channel_checks=checks,
                  independent_symbolic_raw_charged_checks=raw_checks,
                  all_parent_Bose_checks=18*241, parent_filter_nesting=True,
                  direct_pi0_variance_reconstruction=True, pi0_combined_error_latex=sp.latex(pi0_error))
    paths = ['scripts/verify_odd_bounded_amplitudes.py', 'scripts/implementation/odd_bounded_amplitudes.py',
             'data/vector_sector/odd_cubic_parent_kernels.json',
             'data/odd_parity/odd_low_energy_measurements.json']
    result['sha256'] = {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}
    output = ROOT/'data/vector_sector/odd_bounded_amplitude_verification.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(output.relative_to(ROOT))


if __name__ == '__main__':
    main()
