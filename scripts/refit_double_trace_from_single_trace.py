"""Refit the double-trace model from the embedded single-trace solution.

Run: .venv/bin/python scripts/refit_double_trace_from_single_trace.py
At y_V=0 the double-trace source action is exactly the single-trace one.
No sign restriction is imposed on signed source residues.
"""
import hashlib
import json
import math
import numpy as np
import profile_odd_low_energy as P

ROOT = P.ROOT


def main():
    paths = [
        'data/vector_sector/odd_preliminary_fit.json',
        'data/odd_parity/odd_low_energy_measurements.json',
        'data/vector_sector/pdg2024.json',
        'data/vector_sector/stage3_results.json',
        'data/vector_sector/alda_2026nrd_arxiv_v1.json',
        'scripts/refit_double_trace_from_single_trace.py',
        'scripts/profile_odd_low_energy.py',
        'scripts/calibrate_odd_low_energy.py',
    ]
    data_art, data, pdg, stage3, alda = [json.loads((ROOT / p).read_text()) for p in paths[:5]]
    results = data_art['results']
    single = next(r for r in results if r['label'] == 'central_single_trace_negative')
    archived_double = next(r for r in results if r['label'] == 'central_double_trace')
    ass = stage3['assumptions']

    tree = P.C.J.module('alda_2026nrd_neutral_inputs').alda_2026nrd_neutral_inputs(
        mass_scheme='tree', F=P.C.J.F, fpi=P.C.J.F)
    fast = P.C.FastChannels(tree, alda['physical_masses_gev'],
                            np.sqrt(4 * math.pi * ass['alpha_em']))
    single_point = np.asarray(single['even_parameters'], dtype=float)
    embedded_point = np.r_[single_point, 0.]
    single_sector, _ = P.C.J.action_sector(single_point, 'single_trace_negative', ass)
    double_sector, _ = P.C.J.action_sector(embedded_point, 'double_trace', ass)
    np.testing.assert_allclose(double_sector.residues, single_sector.residues,
                               rtol=1e-13, atol=1e-13)

    synthetic = dict(archived_double)
    synthetic['measurement_labels'] = [x['channel'] for x in single['observations']]
    synthetic['even_parameters'] = embedded_point.tolist()
    original_coefficients = np.asarray(list(
        single['reference_backend_representative_coefficients'].values()))
    synthetic['gauge_representative_coefficients'] = (
        P.C.model_maps()['original_11'] @ original_coefficients).tolist()
    synthetic['shared_total_widths_gev'] = single['shared_total_widths_gev']
    likelihood = P.ConditionalLikelihood(synthetic, data, pdg, ass, fast)
    start = likelihood.start.copy()
    initial_residual = likelihood.residual(start)
    n_odd = len(likelihood.obs)
    n_even = 7
    initial_odd_chi2 = float(initial_residual[:n_odd] @ initial_residual[:n_odd])
    initial_even_chi2 = float(initial_residual[n_odd:n_odd + n_even] @
                               initial_residual[n_odd:n_odd + n_even])
    initial_width_chi2 = float(initial_residual[n_odd + n_even:] @
                                initial_residual[n_odd + n_even:])
    initial_chi2 = float(initial_residual @ initial_residual)
    fit = likelihood.fit(start, max_nfev=800)
    summary = likelihood.summary(fit)

    output = dict(
        scope='one conditional double-trace refit initialized at the embedded central negative single-trace result',
        source_sign_policy='signed residues are predictions; no positivity or universal-sign restriction is imposed',
        y_zero_embedding=dict(
            y_V_GeV_inverse=0.,
            single_trace_residues_GeV=list(single_sector.residues),
            double_trace_residues_GeV=list(double_sector.residues),
            residues_equal=True,
            single_trace_fit_chi2=float(single['chi2']),
            double_trace_initial_chi2=initial_chi2,
            initial_chi2_breakdown=dict(
                odd_observables=initial_odd_chi2,
                even_observables=initial_even_chi2,
                shared_total_width_nuisances=initial_width_chi2)),
        observable_definition_comparison=dict(
            agrees_with_single_trace=True,
            odd_observable_count=n_odd,
            odd_observable_labels=likelihood.labels,
            odd_observable_order='included records in data/odd_parity/odd_low_energy_measurements.json order',
            even_observables='rho, Kstar, omega, phi masses and rho, omega, phi leptonic widths',
            shared_nuisances='same total-width nuisance parameters and Gaussian penalties'),
        result=summary,
        archived_double_trace_chi2=float(archived_double['chi2']),
        archived_double_trace_note='Historical result used an additional residue-sign restriction and is not the same fit domain.',
        seeded_result_is_better_than_single=bool(summary['chi2'] <= single['chi2'] + 1e-8),
        global_minimum_proven=False,
        sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
    )
    path = ROOT / 'data/vector_sector/odd_double_trace_seeded_refit.json'
    path.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    print(path.relative_to(ROOT))
    print(f"single-trace fit chi2: {single['chi2']:.12g}")
    print(f"double-trace chi2 at y_V=0: {initial_chi2:.12g}")
    print(f"double-trace refit chi2: {summary['chi2']:.12g}; success={summary['success']}")


if __name__ == '__main__':
    main()
