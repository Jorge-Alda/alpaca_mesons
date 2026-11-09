"""Block-exclusion and predefined sensitivity refits of the joint vector fit.

Run (block exclusions, then variations):
  .venv/bin/python scripts/refit_joint_vector_sensitivity.py --study exclusions --output data/vector_sector/joint_vector_block_exclusion.json
  .venv/bin/python scripts/refit_joint_vector_sensitivity.py --study variations --output data/vector_sector/joint_vector_sensitivity_variations.json

Every refit starts from the best accepted endpoint
(basin_probes/even_shared_widths), uses quadrature orders (2,16,16) with
four-bin blocks, scipy trf with max_nfev=5000, ftol=1e-8, xtol=1e-9,
gtol=1e-7 and the fixed optimizer scales of the comparison ledger, and is
re-evaluated at (4,128,128). Variations are additionally audited with the
quadrature ladder and acceptance rule of compare_joint_vector_endpoints.py.
A variation whose endpoint lies at scaled distance > 1 from the best accepted
endpoint is refit again from the other three accepted endpoints.

Units: GeV, nb; residuals and scaled distances are dimensionless. Parameter
shifts are expressed in the baseline (4,128,128) local errors of the best
accepted endpoint, after converting even-sector coordinates to the baseline
reference point (needed for the second inverse-mass branch).
"""
import argparse
import hashlib
import json
import os
import shlex
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_distribution

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assess_joint_vector_fit_adequacy as adequacy
import compare_joint_vector_endpoints as comparison_module
import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_joint_vector_multistart as multi
import joint_vector_likelihood_options as options_module

ROOT = joint_fit.ROOT
BEST = 'basin_probes/even_shared_widths'
ADEQUACY = 'data/vector_sector/joint_vector_fit_adequacy.json'
FIT_ORDERS = (2, 16, 16)
HIGH_ORDERS = (4, 128, 128)
REFINEMENT_ORDERS = (2, 64, 64)
BLOCK_SIZE = 4
OPTIMIZER_SETTINGS = dict(method='trf', max_nfev=5000, ftol=1e-8, xtol=1e-9, gtol=1e-7)
BASIN_SWITCH_DISTANCE = 1.
SINGULAR_RELATIVE_THRESHOLD = 1e-6
SINGULAR_COMPONENT_THRESHOLD = .3

EXCLUSIONS = [dict(name=f'without_{name}', excluded_blocks=(name,),
                   definition=f'Remove the {name} block ({stop-start} rows) from the residual.')
              for name, start, stop in options_module.BLOCKS]

VARIATIONS = [
    dict(name='babar_upper_energy_1p70', options=dict(babar_upper_energy_gev=1.70),
         definition='Drop BaBar bins with upper edge above 1.70 GeV (baseline 1.78 GeV).',
         justification='Higher excited vector states and omitted K2star/tensor dynamics grow toward the top of the window.'),
    dict(name='babar_upper_energy_1p60', options=dict(babar_upper_energy_gev=1.60),
         definition='Drop BaBar bins with upper edge above 1.60 GeV.',
         justification='As above, with a stricter window.'),
    dict(name='regulator_2p0', options=dict(regulator_gev=2.0),
         definition='Regulator Lambda=2.0 GeV in the amplitude factor exp(-E^6/Lambda^6), E the e+e- energy, which multiplies every BaBar bin integrand squared (baseline 2.5 GeV).',
         justification='The regulator value is a model choice without a measurement.'),
    dict(name='regulator_3p0', options=dict(regulator_gev=3.0),
         definition='Regulator Lambda=3.0 GeV.', justification='As above.'),
    dict(name='tensor_products_zero', options={}, tensor_products_zero=True,
         definition='Fix g_T C_T and g_T D_T to zero (two fewer free coordinates).',
         justification='Tests whether the a2 tensor term is required by the data.'),
    dict(name='source_universal', options=dict(source_hypothesis='universal'),
         start_source=('universal', 'near_ideal'),
         definition='Universal source (z_V=0, one fewer free coordinate); domain F_V>0.',
         justification='Declared source hypothesis of the stage-5 specification.'),
    dict(name='source_single_trace_positive', options=dict(source_hypothesis='single_trace_positive'),
         start_source=('single_trace_positive', 'near_ideal'),
         definition='Single-trace source with positive light and strange residues.',
         justification='Declared source hypothesis of the stage-5 specification.'),
    dict(name='inverse_mass_second_branch', options=dict(inverse_mass_branch='second'),
         start_source=('single_trace_negative', 'second'),
         definition='Second (far from ideal) root of the inverse-mass equations, negative single-trace source; (Delta_V, kappa_V) recentred on the stage-3 root.',
         justification='Both mass branches are retained by the stage-5 specification.'),
    dict(name='phi_kk_correlation_zero', options=dict(phi_kk_correlation=False),
         definition='Set the PDG phi->K+K-/KLKS branching-fraction correlation -0.72 to zero.',
         justification='Covariance-treatment check.'),
    dict(name='phi_electronic_width_average', options=dict(phi_electronic_width='average'),
         definition='phi electronic width PDG OUR AVERAGE 1.251+-0.021 keV instead of OUR EVALUATION 1.27+-0.04 keV.',
         justification='Alternative phi electronic-width summary retained by the stage-5 specification.'),
    dict(name='babar_fractional_systematic_5_percent', options=dict(babar_fractional_systematic=.05),
         definition='Add 5% of the observed BaBar cross section in quadrature to each bin error, uncorrelated between bins.',
         justification='ASSUMPTION: no bin-to-bin systematic is given with the inclusive BaBar tables used here; size chosen a priori.'),
]

# The 36 +-1 sigma mass and width shifts of the fixed excited vector poles and of
# the a2(1320) pole (PDG 2024 uncertainties, data/vector_sector/pdg2024_excited_poles.json).
POLE_SHIFT_VARIATIONS = options_module.pole_shift_variations()
VARIATIONS += POLE_SHIFT_VARIATIONS

NOT_IMPLEMENTED = [
    dict(name='source_single_plus_double_trace',
         reason='deferred by user decision 2026-11-07; candidate model extension for follow-up bead'),
]


def accepted_context():
    """Accepted endpoints, optimizer scales and baseline local errors."""
    comparison = json.loads((ROOT/adequacy.COMPARISON).read_bytes())
    accepted = adequacy.load_accepted_endpoints(comparison)
    if accepted[0]['id'] != BEST:
        raise ValueError('best accepted endpoint changed')
    assessment = json.loads((ROOT/ADEQUACY).read_bytes())
    best_record = next(record for record in assessment['endpoints'] if record['id'] == BEST)
    return dict(accepted=[dict(id=item['id'], endpoint=item['endpoint']) for item in accepted],
                scales=comparison['optimizer_coordinate_scales'],
                local_errors=best_record['local_error_vector'],
                baseline_blocks=best_record['blocks'], baseline_chi2=best_record['chi2'])


def baseline_coordinates(model, values):
    """Express coordinates relative to the baseline even-sector reference point."""
    values = np.asarray(values, dtype=float).copy()
    physical = model.point+joint_fit.EVEN_SCALES*values[:6]
    if getattr(model, 'options', {}).get('source_hypothesis') == 'universal':
        physical[5] = 0.
    reference = np.array([model.preliminary[name] for name in joint_fit.EVEN_NAMES])
    values[:6] = (physical-reference)/joint_fit.EVEN_SCALES
    return values


def variation_start(model, endpoint, variation):
    """Warm start from an accepted endpoint, moved into the variation domain."""
    start = np.asarray(endpoint, dtype=float).copy()
    if 'start_source' in variation:
        hypothesis, branch = variation['start_source']
        if branch == 'second':
            start[2:4] = 0.
        vector_coupling, source_slope = model.stage3_source(hypothesis, branch)
        start[4] = (vector_coupling-model.point[4])/joint_fit.EVEN_SCALES[4]
        start[5] = (source_slope-model.point[5])/joint_fit.EVEN_SCALES[5]
    for index, value in model.fixed_coordinates(variation.get('tensor_products_zero', False)).items():
        start[index] = value
    return start


def near_singular_columns(jacobian, scales, free):
    """Directions of the scaled Jacobian with relative singular value < 1e-6.

    jacobian holds only the free columns, listed by index in free.
    """
    scaled = jacobian*np.asarray(scales)[free]
    norms = np.linalg.norm(scaled, axis=0)
    _, singular, right = np.linalg.svd(scaled, full_matrices=False)
    small = np.flatnonzero(singular < SINGULAR_RELATIVE_THRESHOLD*singular[0])
    directions = [dict(relative_singular_value=float(singular[index]/singular[0]),
                       coordinates={joint_fit.PARAMETER_NAMES[free[column]]: float(right[index, column])
                                    for column in np.flatnonzero(abs(right[index]) > SINGULAR_COMPONENT_THRESHOLD)})
                  for index in small]
    zero = [joint_fit.PARAMETER_NAMES[free[column]] for column in np.flatnonzero(norms == 0)]
    return dict(relative_threshold=SINGULAR_RELATIVE_THRESHOLD,
                component_threshold=SINGULAR_COMPONENT_THRESHOLD,
                relative_singular_values=(singular/singular[0]).tolist(),
                zero_columns=zero, near_singular_directions=directions,
                conditional=bool(zero or directions))


def refit(model, start, fixed, scales, settings, checkpoint_path=None):
    """Minimize over the free coordinates; keep failures.

    With checkpoint_path, every iteration is saved; a later call with the same
    start resumes from the saved point with the remaining evaluation budget.
    """
    lower, upper = joint_fit.parameter_bounds(model)
    free = np.array([index for index in range(32) if index not in fixed])
    start = np.asarray(start, dtype=float)
    record = dict(start=start.tolist(), fixed_coordinates={joint_fit.PARAMETER_NAMES[i]: v for i, v in fixed.items()},
                  free_coordinate_count=int(len(free)), prior_nfev=0, prior_seconds=0.)
    started = time.perf_counter()
    initial = start.copy()
    settings = dict(settings)
    checkpoint = Path(checkpoint_path) if checkpoint_path else None
    if checkpoint and checkpoint.exists():
        saved = json.loads(checkpoint.read_text())
        if saved['start'] == record['start']:
            initial = np.asarray(saved['parameter_values'])
            record.update(prior_nfev=saved['cumulative_nfev'], prior_seconds=saved['cumulative_seconds'],
                          resumed_from_checkpoint=True)
            settings['max_nfev'] = max(1, settings['max_nfev']-saved['cumulative_nfev'])
    try:
        multi.validate_start(model, start, lower, upper)
        multi.validate_start(model, initial, lower, upper)
        def full(free_values):
            values = start.copy()
            values[free] = free_values
            return values
        def callback(intermediate_result):
            if checkpoint:
                state = dict(start=record['start'], parameter_values=full(intermediate_result.x).tolist(),
                             cumulative_nfev=record['prior_nfev']+int(intermediate_result.nfev),
                             cumulative_seconds=record['prior_seconds']+time.perf_counter()-started)
                temporary = checkpoint.with_name(checkpoint.name+'.tmp')
                temporary.write_text(json.dumps(state))
                temporary.replace(checkpoint)
        result = least_squares(lambda x: model.residual(full(x)), initial[free],
                               jac=lambda x: model.jacobian(full(x))[:, free],
                               bounds=(lower[free], upper[free]), callback=callback,
                               x_scale=np.asarray(scales)[free], **settings)
        endpoint = full(result.x)
        record.update(state='finished', endpoint=endpoint.tolist(),
                      fit_order_chi2=float(result.fun @ result.fun), success=bool(result.success),
                      status=int(result.status), message=str(result.message), nfev=int(result.nfev),
                      cumulative_nfev=record['prior_nfev']+int(result.nfev),
                      njev=None if result.njev is None else int(result.njev),
                      optimality=float(result.optimality),
                      active_bounds={joint_fit.PARAMETER_NAMES[free[i]]: int(v)
                                     for i, v in enumerate(result.active_mask) if v},
                      penalty_endpoint=bool(np.all(result.fun == 1e8+np.linalg.norm(endpoint))))
        record['jacobian_conditioning'] = near_singular_columns(result.jac, scales, free)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        record.update(state='error', success=False, error_type=type(error).__name__, message=str(error))
    record['refit_seconds'] = record['prior_seconds']+time.perf_counter()-started
    return record


def high_order_evaluation(model, record, free_count, labels):
    """Re-evaluate a finished refit at the high order; blocks on retained rows."""
    started = time.perf_counter()
    values = np.asarray(record['endpoint'])
    try:
        model.unpack(values)
        full = model.full_residual(values)
        if not np.all(np.isfinite(full)) or np.all(full == 1e8+np.linalg.norm(values)):
            raise ValueError('invalid high-order residual')
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        return dict(valid=False, error_type=type(error).__name__, message=str(error),
                    seconds=time.perf_counter()-started)
    retained = full[model.row_mask]
    value = float(retained @ retained)
    dof = int(model.row_mask.sum())-free_count
    excluded_rows = np.flatnonzero(~model.row_mask)
    return dict(valid=True, chi2=value, residual_count=int(model.row_mask.sum()),
                degrees_of_freedom=dof,
                nominal_p_value=float(chi2_distribution.sf(value, dof)) if dof > 0 else None,
                blocks=adequacy.block_summary(full, mask=model.row_mask),
                removed_row_pulls=[dict(index=int(index), pull=float(full[index]), **labels[index])
                                   for index in excluded_rows],
                removed_rows_chi2=float(full[excluded_rows] @ full[excluded_rows]),
                full_residuals=full.tolist(), seconds=time.perf_counter()-started)


def add_constrained_degrees_of_freedom(refit_record, evaluation):
    """Alternative dof that does not count coordinates with zero Jacobian columns.

    When a removed block was the only one sensitive to some coordinates, those
    columns vanish and the coordinates are unconstrained; the default count
    (retained rows minus free coordinates) can then be zero or negative. This
    alternative count is an assumption and is reported next to the default.
    """
    if not evaluation.get('valid') or 'jacobian_conditioning' not in refit_record:
        return
    zero = len(refit_record['jacobian_conditioning']['zero_columns'])
    dof = evaluation['residual_count']-(refit_record['free_coordinate_count']-zero)
    evaluation.update(zero_column_count=zero, degrees_of_freedom_without_zero_columns=dof,
                      nominal_p_value_without_zero_columns=(
                          float(chi2_distribution.sf(evaluation['chi2'], dof)) if dof > 0 else None))


def make_models(options, orders_list):
    """Varied likelihoods at each requested order, four-bin blocks."""
    return [options_module.VariedLikelihood(*orders, bin_block_size=BLOCK_SIZE, **options)
            for orders in orders_list]


def shifts(model, endpoint, context):
    """Parameter shifts from the best endpoint in baseline local errors."""
    converted = baseline_coordinates(model, endpoint)
    best = np.asarray(context['accepted'][0]['endpoint'])
    shift = (converted-best)/np.asarray(context['local_errors'])
    distance = float(np.linalg.norm((converted-best)/np.asarray(context['scales'])))
    largest = int(np.argmax(abs(shift)))
    return dict(scaled_distance_to_best=distance,
                shifts_in_local_errors=dict(zip(joint_fit.PARAMETER_NAMES, shift.tolist())),
                maximum_absolute_shift=float(abs(shift[largest])),
                maximum_shift_parameter=joint_fit.PARAMETER_NAMES[largest])


def run_exclusion(item, context, models=None, settings=OPTIMIZER_SETTINGS, checkpoint_prefix=None):
    """Refit with one block removed; report remaining and removed rows."""
    started = time.perf_counter()
    low, high = models or make_models(dict(excluded_blocks=item['excluded_blocks']), (FIT_ORDERS, HIGH_ORDERS))
    labels = adequacy.row_labels(high)
    best = context['accepted'][0]
    record = dict(item, excluded_blocks=list(item['excluded_blocks']), start_endpoint=best['id'])
    record['refit'] = refit(low, best['endpoint'], {}, context['scales'], settings,
                            checkpoint_prefix and f'{checkpoint_prefix}.progress.json')
    if record['refit']['state'] == 'finished':
        record['high_order'] = high_order_evaluation(high, record['refit'], 32, labels)
        add_constrained_degrees_of_freedom(record['refit'], record['high_order'])
        record.update(shifts(low, record['refit']['endpoint'], context))
        record['conditional'] = record['refit']['jacobian_conditioning']['conditional']
    record['baseline_blocks'] = context['baseline_blocks']
    record['elapsed_seconds'] = time.perf_counter()-started
    return record


def audit(model, record):
    """Quadrature ladder and acceptance rule of the endpoint comparison."""
    attempt = dict(success=record.get('success'), endpoint=record.get('endpoint'),
                   endpoint_diagnostics=dict(valid=record.get('state') == 'finished'
                                             and not record.get('penalty_endpoint')))
    lower, upper = joint_fit.parameter_bounds(model)
    result = multi.refine_attempt(model, attempt, lower, upper)
    return result, bool(record.get('success') and comparison_module.audit_passed(result))


def run_variation_from(models, variation, start_endpoint, context, settings, with_audit, checkpoint_path=None):
    """One warm-started refit of a variation with evaluation and gates."""
    low, high, refinement = models
    fixed = low.fixed_coordinates(variation.get('tensor_products_zero', False))
    labels = adequacy.row_labels(high)
    run = dict(start_endpoint=start_endpoint['id'])
    run.update(refit(low, variation_start(low, start_endpoint['endpoint'], variation),
                     fixed, context['scales'], settings, checkpoint_path))
    if run['state'] != 'finished':
        run['accepted'] = False
        return run
    run['high_order'] = high_order_evaluation(high, run, run['free_coordinate_count'], labels)
    add_constrained_degrees_of_freedom(run, run['high_order'])
    run.update(shifts(low, run['endpoint'], context))
    if with_audit:
        started = time.perf_counter()
        run['quadrature_refinement'], run['accepted'] = audit(refinement, run)
        run['audit_seconds'] = time.perf_counter()-started
    else:
        run['accepted'] = False
    return run


def run_variation(variation, context, models=None, settings=OPTIMIZER_SETTINGS, with_audit=True,
                  checkpoint_prefix=None):
    """Refit a variation; rerun from the other accepted endpoints if it moved far."""
    started = time.perf_counter()
    models = models or make_models(variation['options'], (FIT_ORDERS, HIGH_ORDERS, REFINEMENT_ORDERS))
    record = dict(variation, implemented=True)
    def checkpoint(index):
        return checkpoint_prefix and f'{checkpoint_prefix}.{index}.progress.json'
    runs = [run_variation_from(models, variation, context['accepted'][0], context, settings, with_audit,
                               checkpoint(0))]
    if runs[0].get('scaled_distance_to_best', 0.) > BASIN_SWITCH_DISTANCE:
        runs += [run_variation_from(models, variation, endpoint, context, settings, with_audit, checkpoint(index))
                 for index, endpoint in enumerate(context['accepted'][1:], start=1)]
    return variation_record(variation, runs, time.perf_counter()-started)


def variation_record(variation, runs, elapsed_seconds):
    """Variation entry with the selected run: lowest high-order chi2, accepted runs first."""
    record = dict(variation, implemented=True, runs=runs)
    valid = [run for run in runs if run.get('high_order', {}).get('valid')]
    accepted = [run for run in valid if run['accepted']]
    pool = accepted or valid
    if pool:
        selected = min(pool, key=lambda run: run['high_order']['chi2'])
        record['selected_run'] = runs.index(selected)
        record['selected_accepted'] = selected['accepted']
    record['basin_switch_reruns'] = len(runs) > 1
    record['elapsed_seconds'] = elapsed_seconds
    return record


def run_single(arguments):
    """Worker entry: one run (start endpoint index) of a variation, saved as its own part."""
    name, index, parts = arguments
    started = time.perf_counter()
    context = accepted_context()
    variation = next(item for item in VARIATIONS if item['name'] == name)
    models = make_models(variation['options'], (FIT_ORDERS, HIGH_ORDERS, REFINEMENT_ORDERS))
    run = run_variation_from(models, variation, context['accepted'][index], context, OPTIMIZER_SETTINGS, True,
                             str(Path(parts)/f'{name}.{index}.progress.json'))
    run['run_seconds'] = time.perf_counter()-started
    multi.write_ledger(Path(parts)/f'{name}.run{index}.json', run)
    return name, index


def run_variations_by_run(names, parts, workers):
    """Schedule single runs: run 0 of every variation, then the basin-switch reruns.

    Same rule as run_variation, but every run is a separate task so all
    workers stay busy. The entry elapsed_seconds is the sum of the run times.
    """
    run_count = len(accepted_context()['accepted'])
    def run_part(name, index):
        return parts/f'{name}.run{index}.json'
    def needed_runs(name):
        first = run_part(name, 0)
        if not first.exists():
            return [0]
        distance = json.loads(first.read_text()).get('scaled_distance_to_best', 0.)
        return list(range(run_count)) if distance > BASIN_SWITCH_DISTANCE else [0]
    def finalize(name):
        indices = needed_runs(name)
        if not all(run_part(name, index).exists() for index in indices):
            return False
        runs = [json.loads(run_part(name, index).read_text()) for index in indices]
        variation = next(item for item in VARIATIONS if item['name'] == name)
        multi.write_ledger(parts/f'{name}.json',
                           variation_record(variation, runs, sum(run['run_seconds'] for run in runs)))
        print('finished', name, flush=True)
        return True
    def submit_missing(executor, name, futures):
        for index in needed_runs(name):
            if not run_part(name, index).exists() and (name, index) not in futures.values():
                futures[executor.submit(run_single, (name, index, str(parts)))] = (name, index)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for name in names:
            if not finalize(name):
                submit_missing(executor, name, futures)
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                name, index = futures.pop(future)
                future.result()
                print('finished run', name, index, flush=True)
                if not finalize(name):
                    submit_missing(executor, name, futures)


def run_named(arguments):
    """Worker entry: run one exclusion or variation and save its part."""
    study, name, part_path = arguments
    context = accepted_context()
    if study == 'exclusions':
        result = run_exclusion(next(item for item in EXCLUSIONS if item['name'] == name), context,
                               checkpoint_prefix=part_path.removesuffix('.json'))
    else:
        result = run_variation(next(item for item in VARIATIONS if item['name'] == name), context,
                               checkpoint_prefix=part_path.removesuffix('.json'))
    multi.write_ledger(part_path, result)
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--study', choices=('exclusions', 'variations'), required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--only', action='append', help='run only these names')
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--parallel-runs', action='store_true',
                        help='variations only: schedule every start endpoint as its own task')
    args = parser.parse_args()
    items = EXCLUSIONS if args.study == 'exclusions' else VARIATIONS
    known = [item['name'] for item in items]
    names = [name for name in args.only if name in known] if args.only else known
    parts = ROOT/(args.output+'.parts')
    parts.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    pending = [(args.study, name, str(parts/f'{name}.json')) for name in names
               if not (parts/f'{name}.json').exists()]
    if args.parallel_runs and args.study == 'variations':
        run_variations_by_run([name for _, name, _ in pending], parts, args.workers)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for name in executor.map(run_named, pending):
                print('finished', name, flush=True)
    context = accepted_context()
    entries = [json.loads((parts/f'{item["name"]}.json').read_text()) for item in items
               if (parts/f'{item["name"]}.json').exists()]
    for entry in entries:
        # Parts written before this field existed get it at assembly.
        if args.study == 'exclusions' and 'refit' in entry and 'high_order' in entry:
            add_constrained_degrees_of_freedom(entry['refit'], entry['high_order'])
        for run in entry.get('runs', []):
            if 'high_order' in run:
                add_constrained_degrees_of_freedom(run, run['high_order'])
    ledger = dict(study=args.study, best_accepted_endpoint=BEST,
                  accepted_endpoints=[item['id'] for item in context['accepted']],
                  baseline_chi2_high_order=context['baseline_chi2'],
                  baseline_blocks=context['baseline_blocks'],
                  fit_quadrature=dict(zip(('energy', 'pair_per_segment', 'angle_per_segment'), FIT_ORDERS),
                                      bin_block_size=BLOCK_SIZE),
                  evaluation_quadrature=dict(zip(('energy', 'pair_per_segment', 'angle_per_segment'), HIGH_ORDERS),
                                             bin_block_size=BLOCK_SIZE),
                  optimizer_settings=OPTIMIZER_SETTINGS,
                  optimizer_coordinate_scales=context['scales'],
                  degrees_of_freedom_counting='retained residual rows minus free coordinates; '
                                              'per-block p-values use the row count (assumption)',
                  shift_definition='(coordinates converted to the baseline reference point - best endpoint) / '
                                   'baseline (4,128,128) local errors',
                  near_singular_definition='SVD of the fit-order Jacobian times optimizer scales; relative singular '
                                           f'value below {SINGULAR_RELATIVE_THRESHOLD}; refit then labelled conditional',
                  entries=entries,
                  input_sha256={path: adequacy.file_hash(path) for path in
                                (adequacy.COMPARISON, ADEQUACY, *multi.INPUTS)},
                  source_sha256=multi.source_hashes(),
                  driver_elapsed_seconds=time.perf_counter()-started,
                  command=shlex.join([sys.executable, *sys.argv]))
    if args.study == 'variations':
        ledger['acceptance_rule'] = ('optimizer success and certified quadrature ladder '
                                     f'{multi.CERTIFICATION_ORDERS} with tolerances {joint_fit.QUADRATURE_TOLERANCES}')
        ledger['basin_switch_rule'] = (f'scaled distance > {BASIN_SWITCH_DISTANCE}: refit from the other three '
                                       'accepted endpoints; selected run = lowest high-order chi2 among accepted runs')
        ledger['not_implemented'] = NOT_IMPLEMENTED
    multi.write_ledger(ROOT/args.output, ledger)
    print(args.output)


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    main()
