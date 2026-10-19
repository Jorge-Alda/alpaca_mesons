"""Continue finite odd-sector searches without changing the likelihood (GeV).

Run: .venv/bin/python scripts/continue_odd_low_energy_profiles.py
Additional numerical starts are a finite search, not a global-minimum proof.
"""
import argparse
import hashlib
import json
import math
import numpy as np
import scipy.linalg as la
from scipy.optimize import least_squares
import profile_odd_low_energy as P

ROOT = P.ROOT


def source_point(likelihood, reference, descriptor, target, starts, max_nfev):
    """Eliminate the same physical source; rebuild odd starts at the target."""
    index = descriptor[1]
    mask = np.arange(len(reference)) != index

    def expand(z):
        v = np.zeros(len(reference))
        v[mask] = z
        low = P.source_quantity(likelihood, v, descriptor)
        v[index] = 1.
        slope = P.source_quantity(likelihood, v, descriptor) - low
        if abs(slope) < 1e-14:
            raise ValueError('source constraint lost sensitivity')
        v[index] = (target - low) / slope
        return v

    def residual(z):
        try:
            return likelihood.residual(expand(z))
        except ValueError:
            return np.full(len(likelihood.obs) + 14, 1e12)

    candidates = []
    for start in starts:
        try:
            constrained = expand(start[mask])
            rows, off, factors = likelihood.blocks(constrained)[3:]
            # At fixed nuisances, use the exact affine amplitude map. Each
            # observed sign gives a deterministic new gauge-coordinate start.
            signs = np.sign(likelihood.predictions(start)[1])
            matrix = np.sqrt(factors)[:, None] * (rows @ likelihood.gauge)
            rhs = signs * np.sqrt(likelihood.y) - np.sqrt(factors) * off
            constrained[likelihood.nn:] = la.lstsq(matrix, rhs, cond=1e-12)[0]
            fit = least_squares(residual, constrained[mask],
                                bounds=(likelihood.lower[mask], likelihood.upper[mask]),
                                max_nfev=max_nfev, x_scale='jac',
                                xtol=1e-10, ftol=1e-10, gtol=1e-7)
            v = expand(fit.x)
            summary = likelihood.summary(type('Fit', (), dict(x=v, success=fit.success,
                nfev=fit.nfev, message=fit.message))())
            defect = abs(P.source_quantity(likelihood, v, descriptor) - target)
            summary.update(valid=bool(fit.success and summary['chi2'] < 1e20 and
                defect < 1e-9 * max(target, 1e-8)), constraint_absolute_defect=defect)
            candidates.append(summary)
        except (ValueError, np.linalg.LinAlgError) as error:
            candidates.append(dict(valid=False, reason=str(error)))
    candidates.sort(key=lambda x: (not x['valid'], x.get('chi2', float('inf'))))
    return candidates[0] if candidates else dict(valid=False, reason='no feasible start')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--count', type=int, default=0)
    parser.add_argument('--max-nfev', type=int, default=800)
    args = parser.parse_args()
    paths = ['data/vector_sector/odd_low_energy_profiles.json',
             'data/vector_sector/odd_low_energy_calibration.json',
             'data/odd_parity/odd_low_energy_measurements.json',
             'data/vector_sector/pdg2024.json', 'data/vector_sector/stage3_results.json',
             'data/vector_sector/alda_arxiv_v1.json',
             'scripts/continue_odd_low_energy_profiles.py',
             'scripts/profile_odd_low_energy.py', 'scripts/calibrate_odd_low_energy.py']
    profiles, calibration, data, pdg, stage3, alda = [json.loads((ROOT / p).read_text()) for p in paths[:6]]
    ass = stage3['assumptions']
    tree = P.C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree', F=P.C.J.F, fpi=P.C.J.F)
    fast = P.C.FastChannels(tree, alda['physical_masses_gev'], np.sqrt(4 * math.pi * ass['alpha_em']))
    selected = [(i, r) for i, r in enumerate(profiles['results']) if
                r['second_batch_continued_improvement'] or
                any(not q['valid'] for q in r.get('one_sided_source_profile', {}).get('points', []))]
    selected = selected[args.start:args.start + (args.count or len(selected))]
    output = dict(scope='finite continuation of flagged sign searches and unresolved source points',
        global_nonlinear_minimum_proven=False, max_nfev=args.max_nfev,
        policy='refit saved best plus next 16 fixed-point signs (ranks 33-48); source target signed-amplitude least-squares starts',
        sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}, results=[])
    suffix = f'_{args.start}_{args.count}' if args.start or args.count else ''
    path = ROOT / f'data/vector_sector/odd_low_energy_continuation{suffix}.json'
    for i, old in selected:
        base = next(r for r in calibration['results'] if r['model'] == 'soft_leading' and
                    r['scenario_index'] == old['scenario_index'] and r['control'] == old['control'])
        likelihood = P.ConditionalLikelihood(base, data, pdg, ass, fast)
        reference = likelihood.rebase(old['best'])
        summaries = [old['best']]
        if old['second_batch_continued_improvement']:
            starts = likelihood.sign_starts(48)[0][32:48]
            for start in [reference] + starts:
                fit = likelihood.fit(start, max_nfev=args.max_nfev)
                summaries.append(likelihood.summary(fit))
        valid = [s for s in summaries if s['success'] and s['chi2'] < 1e20]
        best = min(valid, key=lambda s: s['chi2'])
        record = dict(pair_index=i, scenario_index=old['scenario_index'], control=old['control'],
                      old_chi2=old['best']['chi2'], best=best, search_minima=summaries,
                      improvement=old['best']['chi2'] - best['chi2'], source_points=[])
        source = old.get('one_sided_source_profile')
        if source:
            descriptor = P.closest_source_quantity(likelihood, reference)
            for point in source['points']:
                if point['valid']:
                    continue
                candidates = [reference, likelihood.rebase(best)]
                for previous in source['points']:
                    if previous['valid']:
                        candidates.append(np.asarray(previous['coordinates']))
                q = source_point(likelihood, reference, descriptor,
                                 point['signed_source_target_gev'], candidates, args.max_nfev)
                q.update(multiplier=point['multiplier'], signed_source_target_gev=point['signed_source_target_gev'])
                record['source_points'].append(q)
        record.update(maximum_every_evaluation_image_defect=likelihood.maximum_image_defect,
                      domain_rejections=likelihood.domain_rejections)
        output['results'].append(record)
        output['complete'] = len(output['results']) == len(selected)
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
        print(f'{path.relative_to(ROOT)}: pair {i} complete', flush=True)


if __name__ == '__main__':
    main()
