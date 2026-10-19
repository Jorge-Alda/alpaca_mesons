"""One-sided source profiles in physical source coordinates (GeV).

Run: .venv/bin/python scripts/continue_odd_low_energy_source_profiles.py
The likelihood and strict source hypotheses are inherited unchanged.
"""
from decimal import Decimal, localcontext
import hashlib
import json
import math
import numpy as np
import scipy.linalg as la
from scipy.optimize import least_squares
import profile_odd_low_energy as P

ROOT = P.ROOT


def source_constraint_certificate(likelihood, v, target):
    """Independent decimal reconstruction plus a binary64 forward-error bound.

    For n elementary correctly rounded operations, gamma_n=n*u/(1-n*u).
    The map uses fewer than 32 operations (exp excluded: the represented F_l
    is the free input). Sum absolute affine/source terms to bound cancellation.
    No statistical or physics tolerance is added to this arithmetic bound.
    """
    point, scales = likelihood.point, likelihood.scales
    xs = likelihood.ass['xs_gev2']
    with localcontext() as context:
        context.prec = 80
        D = Decimal.from_float
        p4 = D(float(point[4])) + D(float(scales[4])) * D(float(v[4]))
        p5 = D(float(point[5])) + D(float(scales[5])) * D(float(v[5]))
        recovered = -(p4 + 2 * p5 * D(float(xs)))
        defect = float(abs(recovered - D(float(target))))
    # Include source-coordinate construction, inversion into standardized
    # nuisance coordinates, and their later affine reconstruction.
    magnitude = (abs(point[4]) + abs(scales[4]*v[4]) +
                 2*abs(xs)*(abs(point[5]) + abs(scales[5]*v[5])) + abs(target))
    u = np.finfo(float).eps / 2
    gamma32 = 32*u/(1-32*u)
    bound = gamma32*magnitude
    valid = bool(defect <= bound and target > bound and recovered > 0)
    return dict(valid=valid, independent_decimal_signed_source_gev=str(recovered),
                absolute_defect_gev=defect, forward_error_bound_gev=bound,
                unit_roundoff=u, operation_bound=32,
                strict_positive_source_certified=bool(float(recovered)>bound),
                target_relative_arithmetic_bound=bound/target)


def source_point(likelihood, reference, target, starts, max_nfev=800):
    """Use log(F_l / GeV) freely and fix F_s=-target strictly."""
    if likelihood.hyp != 'single_trace_negative':
        raise ValueError('only the unresolved negative single-trace source family is implemented')
    xl, xs = likelihood.ass['xl_gev2'], likelihood.ass['xs_gev2']
    mask = np.arange(len(reference)) != 4
    free_index = int(np.flatnonzero(np.arange(len(reference))[mask] == 5)[0])

    def expand(z):
        v = np.zeros(len(reference)); v[mask] = z
        fl = math.exp(z[free_index])
        cf = (-target-fl)/(2*(xs-xl))
        fv = fl-2*cf*xl
        v[4] = (fv-likelihood.point[4])/likelihood.scales[4]
        v[5] = (cf-likelihood.point[5])/likelihood.scales[5]
        return v

    def residual(z):
        try:
            return likelihood.residual(expand(z))
        except (ValueError, OverflowError):
            return np.full(len(likelihood.obs)+14, 1e12)

    lower, upper = likelihood.lower[mask].copy(), likelihood.upper[mask].copy()
    # Binary64 exp range is numerical, not a physical prior. Record any hit.
    lower[free_index], upper[free_index] = -700., 700.
    candidates = []
    for start, regenerate_odd in ((s, mode) for s in starts for mode in (False, True)):
        try:
            p = likelihood.point+likelihood.scales*start[:likelihood.np]
            fl = p[4]+2*p[5]*xl
            if fl <= 0:
                raise ValueError('initial F_l is not positive')
            z = start[mask].copy(); z[free_index] = math.log(fl)
            constrained = expand(z)
            rows, off, factors = likelihood.blocks(constrained)[3:]
            signs = np.sign(likelihood.predictions(start)[1])
            matrix = np.sqrt(factors)[:, None]*(rows@likelihood.gauge)
            rhs = signs*np.sqrt(likelihood.y)-np.sqrt(factors)*off
            if regenerate_odd:
                constrained[likelihood.nn:] = la.lstsq(matrix, rhs, cond=1e-12)[0]
            z = constrained[mask]; z[free_index] = math.log(fl)
            fit = least_squares(residual, z, bounds=(lower, upper), max_nfev=max_nfev,
                                x_scale='jac', xtol=1e-10, ftol=1e-10, gtol=1e-7)
            v = expand(fit.x)
            summary = likelihood.summary(type('Fit', (), dict(x=v, success=fit.success,
                nfev=fit.nfev, message=fit.message))())
            certificate = source_constraint_certificate(likelihood, v, target)
            summary.update(source_constraint_certificate=certificate,
                           regenerated_odd_start=regenerate_odd,
                           free_log_source_coordinate=float(fit.x[free_index]),
                           numerical_exp_bound_active=bool(abs(fit.x[free_index])>699.999),
                           valid=bool(fit.success and summary['chi2']<1e20 and
                               certificate['valid'] and abs(fit.x[free_index])<699.999))
            candidates.append(summary)
        except (ValueError, np.linalg.LinAlgError, OverflowError) as error:
            candidates.append(dict(valid=False, reason=str(error)))
    candidates.sort(key=lambda q: (not q['valid'], q.get('chi2', float('inf'))))
    best = dict(candidates[0])
    best['candidate_fits'] = candidates
    return best


def main():
    paths = ['data/vector_sector/odd_low_energy_profiles.json',
             'data/vector_sector/odd_low_energy_calibration.json',
             'data/odd_parity/odd_low_energy_measurements.json',
             'data/vector_sector/pdg2024.json', 'data/vector_sector/stage3_results.json',
             'data/vector_sector/alda_arxiv_v1.json',
             'data/vector_sector/odd_low_energy_continuation.json',
             'scripts/continue_odd_low_energy_source_profiles.py',
             'scripts/profile_odd_low_energy.py', 'scripts/calibrate_odd_low_energy.py']
    profiles, calibration, data, pdg, stage3, alda, continuation = [json.loads((ROOT/p).read_text()) for p in paths[:7]]
    ass = stage3['assumptions']
    tree = P.C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree', F=P.C.J.F, fpi=P.C.J.F)
    fast = P.C.FastChannels(tree, alda['physical_masses_gev'], np.sqrt(4*math.pi*ass['alpha_em']))
    output = dict(scope='unchanged one-sided source targets with physical source coordinates',
                  source_coordinate='free log(F_l/GeV), fixed F_s=-target',
                  sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, results=[])
    path = ROOT/'data/vector_sector/odd_low_energy_source_continuation.json'
    for record in continuation['results']:
        unresolved = [q for q in record['source_points'] if not q['valid']]
        if not unresolved:
            continue
        old = profiles['results'][record['pair_index']]
        base = next(r for r in calibration['results'] if r['model']=='soft_leading' and
                    r['scenario_index']==old['scenario_index'] and r['control']==old['control'])
        likelihood = P.ConditionalLikelihood(base, data, pdg, ass, fast)
        reference = likelihood.rebase(old['best'])
        starts = [reference, likelihood.rebase(record['best'])]
        result = dict(pair_index=record['pair_index'], scenario_index=record['scenario_index'], control=record['control'], points=[])
        for point in unresolved:
            seeds = starts.copy()
            if 'coordinates' in point:
                seeds.append(np.asarray(point['coordinates']))
            if result['points'] and result['points'][-1]['valid']:
                seeds.append(np.asarray(result['points'][-1]['coordinates']))
            fitted = source_point(likelihood, reference, point['signed_source_target_gev'], seeds)
            fitted.update(multiplier=point['multiplier'], signed_source_target_gev=point['signed_source_target_gev'])
            result['points'].append(fitted)
            print(f'{path.relative_to(ROOT)}: pair {record["pair_index"]}, target {point["multiplier"]} complete', flush=True)
        result['maximum_every_evaluation_image_defect'] = likelihood.maximum_image_defect
        output['results'].append(result)
        output['complete'] = False
        path.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    output['complete'] = True
    path.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
