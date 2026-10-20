"""Re-anchor source-discovered references in the unconstrained likelihood.

Run: .venv/bin/python scripts/reanchor_odd_low_energy_profiles.py
GeV conventions and hypothesis domains are inherited unchanged. Both source
magnitudes vary freely; finite searches do not establish a global minimum.
"""
import hashlib
import json
import math
import numpy as np
from scipy.optimize import least_squares
import profile_odd_low_energy as P
from continue_odd_low_energy_source_profiles import source_constraint_certificate

ROOT = P.ROOT


def physical_source_fit(likelihood, start, max_nfev=1200):
    """Same unconstrained residual, with positive F_l and signed F_s logs."""
    if likelihood.hyp != 'single_trace_negative':
        raise ValueError('implemented only for affected negative single-trace family')
    xl, xs = likelihood.ass['xl_gev2'], likelihood.ass['xs_gev2']
    p = likelihood.point+likelihood.scales*start[:likelihood.np]
    fl, signed_fs = p[4]+2*p[5]*xl, -(p[4]+2*p[5]*xs)
    z = start.copy(); z[4], z[5] = math.log(fl), math.log(signed_fs)
    lower, upper = likelihood.lower.copy(), likelihood.upper.copy()
    lower[4:6], upper[4:6] = -700., 700.

    def expand(z):
        fl, signed_fs = math.exp(z[4]), math.exp(z[5])
        cf = (-signed_fs-fl)/(2*(xs-xl))
        fv = fl-2*cf*xl
        v = z.copy()
        v[4], v[5] = (fv-likelihood.point[4])/likelihood.scales[4], (cf-likelihood.point[5])/likelihood.scales[5]
        return v

    def residual(z):
        try:
            return likelihood.residual(expand(z))
        except (ValueError, OverflowError):
            return np.full(len(likelihood.obs)+14, 1e12)

    fit = least_squares(residual, z, bounds=(lower, upper), max_nfev=max_nfev,
                        x_scale='jac', xtol=1e-10, ftol=1e-10, gtol=1e-7)
    v = expand(fit.x)
    summary = likelihood.summary(type('Fit', (), dict(x=v, success=fit.success,
        nfev=fit.nfev, message=fit.message))())
    cert = source_constraint_certificate(likelihood, v, math.exp(fit.x[5]))
    summary.update(method='free physical source logarithms',
                   free_log_sources=fit.x[4:6].tolist(),
                   source_reconstruction_certificate=cert,
                   numerical_exp_bound_active=bool(np.any(abs(fit.x[4:6])>699.999)),
                   valid=bool(fit.success and summary['chi2']<1e20 and cert['valid'] and
                              np.all(abs(fit.x[4:6])<699.999)))
    return summary


def main():
    paths = ['data/vector_sector/odd_low_energy_profiles.json',
             'data/vector_sector/odd_low_energy_calibration.json',
             'data/odd_parity/odd_low_energy_measurements.json',
             'data/vector_sector/pdg2024.json', 'data/vector_sector/stage3_results.json',
             'data/vector_sector/alda_2026nrd_arxiv_v1.json',
             'data/vector_sector/odd_low_energy_continuation.json',
             'data/vector_sector/odd_low_energy_source_continuation.json',
             'scripts/reanchor_odd_low_energy_profiles.py',
             'scripts/continue_odd_low_energy_source_profiles.py',
             'scripts/profile_odd_low_energy.py', 'scripts/calibrate_odd_low_energy.py']
    profiles, calibration, data, pdg, stage3, alda, continuation, source = [json.loads((ROOT/p).read_text()) for p in paths[:8]]
    ass = stage3['assumptions']
    tree = P.C.J.module('alda_2026nrd_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree', F=P.C.J.F, fpi=P.C.J.F)
    fast = P.C.FastChannels(tree, alda['physical_masses_gev'], np.sqrt(4*math.pi*ass['alpha_em']))
    output = dict(scope='finite unconstrained re-anchoring from accepted source minima',
                  global_nonlinear_minimum_proven=False,
                  sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, results=[])
    path = ROOT/'data/vector_sector/odd_low_energy_profile_reanchoring.json'
    for record in source['results']:
        old = profiles['results'][record['pair_index']]
        previous = next(r for r in continuation['results'] if r['pair_index']==record['pair_index'])
        accepted = [q for q in record['points'] if q['valid'] and q['chi2']<previous['best']['chi2']-1e-4]
        if not accepted:
            continue
        base = next(r for r in calibration['results'] if r['model']=='soft_leading' and
                    r['scenario_index']==old['scenario_index'] and r['control']==old['control'])
        likelihood = P.ConditionalLikelihood(base, data, pdg, ass, fast)
        fits = []
        for point in accepted:
            start = np.asarray(point['coordinates'])
            # The feasible source point is itself an allowed unconstrained
            # reference, even if a subsequent optimizer fails to improve it.
            baseline = dict(point); baseline.pop('candidate_fits', None)
            baseline.update(method='feasible source-discovered incumbent', source_multiplier=point['multiplier'])
            fits.append(baseline)
            try:
                fit = likelihood.fit(start, max_nfev=1200)
                summary = likelihood.summary(fit)
                summary.update(method='original standardized coordinates', source_multiplier=point['multiplier'],
                               valid=bool(fit.success and summary['chi2']<1e20))
                fits.append(summary)
            except (ValueError, np.linalg.LinAlgError) as error:
                fits.append(dict(valid=False, method='original standardized coordinates', reason=str(error)))
            try:
                summary = physical_source_fit(likelihood, start)
                summary['source_multiplier'] = point['multiplier']
                fits.append(summary)
            except (ValueError, np.linalg.LinAlgError, OverflowError) as error:
                fits.append(dict(valid=False, method='free physical source logarithms', reason=str(error)))
        best = min((q for q in fits if q['valid']), key=lambda q:q['chi2'])
        result = dict(pair_index=record['pair_index'], scenario_index=record['scenario_index'], control=record['control'],
                      original_reference_chi2=old['best']['chi2'], previous_continuation_chi2=previous['best']['chi2'],
                      best=best, fits=fits, improvement=previous['best']['chi2']-best['chi2'],
                      additional_improvement_beyond_source=min(q['chi2'] for q in accepted)-best['chi2'],
                      old_source_profile_reference_resolved=False,
                      profiles_require_recomputation=True,
                      maximum_every_evaluation_image_defect=likelihood.maximum_image_defect)
        output['results'].append(result)
        output['complete'] = False
        path.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
        print(f'{path.relative_to(ROOT)}: pair {record["pair_index"]} complete', flush=True)
    output['complete'] = True
    path.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
