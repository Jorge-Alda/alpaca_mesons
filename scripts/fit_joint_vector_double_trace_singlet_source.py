"""Joint vector fit with the double-trace singlet-source coordinate y_V switched on.

Dry run (no optimizer; builds the model, checks it and the source hashes,
evaluates the residual and Jacobian once at every start, prints a summary):
  .venv/bin/python scripts/fit_joint_vector_double_trace_singlet_source.py --dry-run
Fit (writes the ledger with bounds, starts and scales before any optimizer call):
  .venv/bin/python scripts/fit_joint_vector_double_trace_singlet_source.py --output data/vector_sector/joint_vector_double_trace_singlet_source_fit.json --workers 2

Model: DoubleTraceSingletSourceLikelihood of
joint_vector_double_trace_singlet_source_option.py with
double_trace_singlet_source=True: 33 coordinates (the 32 baseline coordinates
and y_V in GeV^-1, bounds [-0.5, 0.5] GeV^-1), 75 residual rows, 75 - 33 = 42
nominal degrees of freedom.

Starts (the plan record
memory/rcht/odd-parity/joint-vector-double-trace-singlet-source-fit-plan.md
states them and the bounds before any refit):
  * y_V = 0, y_V free, from each of the four accepted endpoints of alp-uzh.6
    (data/vector_sector/joint_vector_endpoint_comparison.json);
  * y_V = +0.022 and -0.022 GeV^-1, y_V free, from the best accepted endpoint
    (magnitude of the stage-3 near-ideal double-trace value, both signs);
  * zero_control: y_V = 0 held fixed from the best accepted endpoint; it must
    reproduce the baseline endpoint (scaled distance <= 0.01 and
    |chi2(4,128,128) - baseline| <= 0.01).
Each start is one refit (no basin-switch rerun: the four accepted endpoints
are all used as starts already).

Gates (same as refit_joint_vector_sensitivity.py and the alp-uzh.6 ledgers):
  * scipy least_squares, OPTIMIZER_SETTINGS (trf, max_nfev=5000, ftol=1e-8,
    xtol=1e-9, gtol=1e-7), fixed optimizer scales, quadrature (2,16,16) with
    four-bin blocks; status, message, nfev, njev, optimality, active bounds
    and near-singular directions recorded; failures kept; per-run checkpoints.
  * Re-evaluation at (4,128,128): chi2, dof, nominal p, block chi2.
  * Local errors sqrt(diag((J^T J)^-1)) of the free coordinates from the
    (4,128,128) Jacobian (linearized, not profile intervals).
  * Quadrature ladder of the endpoint comparison from (2,64,64); accepted =
    optimizer success and certified ladder.
  * Selected run: lowest (4,128,128) chi2 among accepted runs with y_V free.
  * Source-hash gate of fit_joint_vector_isospin_breaking_vpp.source_hash_gate.
  * Adequacy criterion: total nominal p > 0.01 and no block with nominal
    p < 1e-3 (row count as block dof), at (4,128,128).

Choices (stated, not results): optimizer scale of y_V = reciprocal Jacobian
column norm at the zero_control start at the fit order (same definition as the
baseline scales); scaled distances and shifts in baseline local errors use the
32 baseline coordinates only, y_V is reported in GeV^-1.

Units: GeV, nb; y_V in GeV^-1; residuals, chi2 and scaled distances
dimensionless. Conventions of Alda:2026nrd.
"""
import argparse
import json
import os
import shlex
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assess_joint_vector_fit_adequacy as adequacy  # noqa: E402
import compare_joint_vector_endpoints as comparison_module  # noqa: E402
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_joint_vector_isospin_breaking_vpp as vpp_driver  # noqa: E402
import fit_joint_vector_multistart as multi  # noqa: E402
import joint_vector_double_trace_singlet_source_option as option  # noqa: E402
import refit_joint_vector_sensitivity as sensitivity  # noqa: E402

ROOT = joint_fit.ROOT
OUTPUT = 'data/vector_sector/joint_vector_double_trace_singlet_source_fit.json'
BEST = sensitivity.BEST
FIT_ORDERS = sensitivity.FIT_ORDERS
HIGH_ORDERS = sensitivity.HIGH_ORDERS
REFINEMENT_ORDERS = sensitivity.REFINEMENT_ORDERS
BLOCK_SIZE = sensitivity.BLOCK_SIZE
OPTIMIZER_SETTINGS = sensitivity.OPTIMIZER_SETTINGS
MODEL_OPTIONS = dict(double_trace_singlet_source=True)
BASELINE_COUNT = option.BASELINE_COORDINATE_COUNT
COORDINATE_COUNT = BASELINE_COUNT+1
RESIDUAL_COUNT = 75
CONTROL = 'zero_control'


def slug(identifier):
    """Endpoint identifier as a file-name fragment."""
    return identifier.replace('/', '__')


def start_definitions(context):
    """Every start: endpoint identifier, y_V start in GeV^-1, whether y_V is held fixed."""
    accepted = [item['id'] for item in context['accepted']]
    if accepted[0] != BEST:
        raise ValueError('best accepted endpoint changed')
    starts = [dict(name=f'y_zero_from_{slug(identifier)}', start_endpoint=identifier, coefficient_gev_inverse=0.,
                   coefficient_fixed=False,
                   definition=f'y_V = 0 GeV^-1, free; other coordinates at the accepted endpoint {identifier}.')
              for identifier in accepted]
    for sign, label in ((1., 'plus'), (-1., 'minus')):
        value = sign*option.STAGE3_SCALE_GEV_INVERSE
        starts.append(dict(name=f'y_{label}_stage3_scale_from_best', start_endpoint=BEST,
                           coefficient_gev_inverse=value, coefficient_fixed=False,
                           definition=(f'y_V = {value:+g} GeV^-1 (magnitude of the stage-3 near-ideal double-trace '
                                       f'value), free; other coordinates at {BEST}.')))
    starts.append(dict(name=CONTROL, start_endpoint=BEST, coefficient_gev_inverse=0., coefficient_fixed=True,
                       definition=f'y_V = 0 held fixed; other coordinates at {BEST}. Must reproduce the baseline.'))
    return starts


def make_models(orders_list):
    """Option likelihoods (option on) at each order, four-bin blocks."""
    return [option.DoubleTraceSingletSourceLikelihood(*orders, bin_block_size=BLOCK_SIZE, **MODEL_OPTIONS)
            for orders in orders_list]


def check_model(model):
    """Check option, coordinate count and bounds; return a summary."""
    names = model.parameter_names
    if not model.double_trace_singlet_source or tuple(names[BASELINE_COUNT:]) != option.ADDED_COORDINATE_NAMES:
        raise ValueError('double-trace singlet-source option with one added coordinate is required')
    lower, upper = model.parameter_bounds()
    baseline_lower, baseline_upper = joint_fit.parameter_bounds(model)
    bound = option.COEFFICIENT_BOUND_GEV_INVERSE
    if not (np.array_equal(lower[:BASELINE_COUNT], baseline_lower)
            and np.array_equal(upper[:BASELINE_COUNT], baseline_upper)
            and lower[BASELINE_COUNT] == -bound and upper[BASELINE_COUNT] == bound):
        raise ValueError('bounds differ from the baseline bounds plus the y_V bounds')
    return dict(double_trace_singlet_source=True, coordinate_count=len(names), residual_count=RESIDUAL_COUNT,
                nominal_degrees_of_freedom=RESIDUAL_COUNT-len(names),
                added_coordinate_bounds_gev_inverse={names[BASELINE_COUNT]: [float(lower[BASELINE_COUNT]),
                                                                             float(upper[BASELINE_COUNT])]})


def validate_point(model, values, lower, upper):
    """33 finite coordinates inside the bounds, on the likelihood's physical branch."""
    values = np.asarray(values, dtype=float)
    if values.shape != (COORDINATE_COUNT,) or not np.all(np.isfinite(values)):
        raise ValueError(f'point must have {COORDINATE_COUNT} finite coordinates')
    if np.any(values < lower) or np.any(values > upper):
        raise ValueError('point is outside parameter bounds')
    model.unpack(values)
    return values


def endpoint_by_identifier(context, identifier):
    """32 baseline coordinates of an accepted endpoint."""
    return next(item['endpoint'] for item in context['accepted'] if item['id'] == identifier)


def start_from(model, context, definition):
    """(33-coordinate start, fixed coordinates) of a start definition."""
    values = model.extend(endpoint_by_identifier(context, definition['start_endpoint']),
                          definition['coefficient_gev_inverse'])
    validate_point(model, values, *model.parameter_bounds())
    fixed = {BASELINE_COUNT: float(definition['coefficient_gev_inverse'])} if definition['coefficient_fixed'] else {}
    return values, fixed


def optimizer_scales(context, control_jacobian):
    """Comparison-ledger scales plus the reciprocal y_V column norm at the zero_control start."""
    norm = float(np.linalg.norm(np.asarray(control_jacobian)[:, BASELINE_COUNT]))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError('zero or nonfinite y_V column at the zero_control start')
    return np.r_[np.asarray(context['scales'], dtype=float), 1/norm]


def refit(model, start, fixed, scales, settings, checkpoint_path=None):
    """33-coordinate version of refit_joint_vector_sensitivity.refit (same record keys)."""
    lower, upper = model.parameter_bounds()
    names = model.parameter_names
    start = np.asarray(start, dtype=float)
    free = np.array([index for index in range(len(start)) if index not in fixed])
    record = dict(start=start.tolist(), fixed_coordinates={names[i]: v for i, v in fixed.items()},
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
        validate_point(model, start, lower, upper)
        validate_point(model, initial, lower, upper)

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
                      active_bounds={names[free[i]]: int(v) for i, v in enumerate(result.active_mask) if v},
                      penalty_endpoint=bool(np.all(result.fun == 1e8+np.linalg.norm(endpoint))))
        with vpp_driver.extended_parameter_names(names):
            record['jacobian_conditioning'] = sensitivity.near_singular_columns(result.jac, scales, free)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        record.update(state='error', success=False, error_type=type(error).__name__, message=str(error))
    record['refit_seconds'] = record['prior_seconds']+time.perf_counter()-started
    return record


def shifts(model, endpoint, context):
    """Baseline-coordinate shifts in baseline local errors; y_V in GeV^-1."""
    converted = sensitivity.baseline_coordinates(model, endpoint)
    best = np.asarray(endpoint_by_identifier(context, BEST))
    shift = (converted[:BASELINE_COUNT]-best)/np.asarray(context['local_errors'])
    distance = float(np.linalg.norm((converted[:BASELINE_COUNT]-best)/np.asarray(context['scales'])))
    largest = int(np.argmax(abs(shift)))
    names = option.BASELINE_PARAMETER_NAMES
    return dict(scaled_distance_to_best=distance, shifts_in_local_errors=dict(zip(names, shift.tolist())),
                maximum_absolute_shift=float(abs(shift[largest])), maximum_shift_parameter=names[largest],
                coefficient_y_V_gev_inverse=float(np.asarray(endpoint)[BASELINE_COUNT]))


def audit(model, run):
    """Quadrature ladder of the endpoint comparison on 33 coordinates; acceptance flag."""
    if not run.get('success') or run.get('state') != 'finished' or run.get('penalty_endpoint'):
        return dict(status='not_evaluated', passed=False, reason='optimizer did not succeed at a valid endpoint'), False
    try:
        values = validate_point(model, run['endpoint'], *model.parameter_bounds())
        started = time.perf_counter()
        with vpp_driver.extended_parameter_names(model.parameter_names):
            result = joint_fit.audit_quadrature(model, values, multi.CERTIFICATION_ORDERS,
                                                parameter_scales=np.ones(COORDINATE_COUNT))
        result['elapsed_seconds'] = time.perf_counter()-started
        result['bin_block_size'] = model.bin_block_size
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        result = dict(status='unevaluated', passed=False, error_type=type(error).__name__, reason=str(error))
    return result, bool(run.get('success') and comparison_module.audit_passed(result))


def run_start(definition, context, scales, settings=OPTIMIZER_SETTINGS, models=None, checkpoint_path=None):
    """One refit of a start with all gates; failures are kept."""
    started = time.perf_counter()
    low, high, refinement = models or make_models((FIT_ORDERS, HIGH_ORDERS, REFINEMENT_ORDERS))
    values, fixed = start_from(low, context, definition)
    run = dict(definition, run_kind='refit')
    run.update(refit(low, values, fixed, scales, settings, checkpoint_path))
    if run['state'] != 'finished':
        run.update(accepted=False, elapsed_seconds=time.perf_counter()-started)
        return run
    run['high_order'] = sensitivity.high_order_evaluation(high, run, run['free_coordinate_count'],
                                                          adequacy.row_labels(high))
    sensitivity.add_constrained_degrees_of_freedom(run, run['high_order'])
    run['local_errors'] = vpp_driver.local_error_record(high, run)
    run.update(shifts(low, run['endpoint'], context))
    audit_started = time.perf_counter()
    run['quadrature_refinement'], run['accepted'] = audit(refinement, run)
    run['audit_seconds'] = time.perf_counter()-audit_started
    run['adequacy'] = vpp_driver.adequacy_verdict(run['high_order'])
    if definition['coefficient_fixed']:
        run['reproduces_baseline_endpoint'] = vpp_driver.control_reproduces_baseline(run, context)
    run['elapsed_seconds'] = time.perf_counter()-started
    return run


def selected_run(runs):
    """Index of the lowest (4,128,128) chi2 among accepted runs with y_V free (None if there is none)."""
    pool = [(index, run) for index, run in enumerate(runs)
            if not run['coefficient_fixed'] and run.get('accepted') and run.get('high_order', {}).get('valid')]
    if not pool:
        return None
    return min(pool, key=lambda item: item[1]['high_order']['chi2'])[0]


def prepare_ledger(model, context, scales, starts, gate, settings=OPTIMIZER_SETTINGS):
    """Ledger written before any optimizer call: model, bounds, starts, scales and gates."""
    lower, upper = model.parameter_bounds()
    quadrature = dict(zip(('energy', 'pair_per_segment', 'angle_per_segment'), FIT_ORDERS), bin_block_size=BLOCK_SIZE)
    return dict(
        study='double_trace_singlet_source', state='prepared', model_options=MODEL_OPTIONS,
        best_accepted_endpoint=BEST, accepted_endpoints=[item['id'] for item in context['accepted']],
        baseline_chi2_high_order=context['baseline_chi2'], baseline_blocks=context['baseline_blocks'],
        parameter_names=list(model.parameter_names),
        bounds=dict(lower=multi.json_bounds(lower), upper=multi.json_bounds(upper)),
        starts=starts, residual_count=RESIDUAL_COUNT, nominal_degrees_of_freedom=RESIDUAL_COUNT-COORDINATE_COUNT,
        fit_quadrature=quadrature, evaluation_quadrature=dict(zip(quadrature, HIGH_ORDERS), bin_block_size=BLOCK_SIZE),
        refinement_quadrature=dict(zip(quadrature, REFINEMENT_ORDERS), bin_block_size=BLOCK_SIZE,
                                   certification_orders=multi.CERTIFICATION_ORDERS),
        optimizer_settings=dict(settings), optimizer_coordinate_scales=list(map(float, scales)),
        added_coordinate_scale_definition=('reciprocal Jacobian column norm of y_V at the zero_control start, fit '
                                           'order (choice); first 32 from the endpoint comparison'),
        degrees_of_freedom_counting='retained residual rows minus free coordinates; per-block p-values use the '
                                    'row count (assumption)',
        shift_definition='(baseline coordinates - best endpoint) / baseline (4,128,128) local errors; y_V in GeV^-1',
        near_singular_definition=('SVD of the fit-order Jacobian times optimizer scales; relative singular value '
                                  f'below {sensitivity.SINGULAR_RELATIVE_THRESHOLD}; run then labelled conditional'),
        acceptance_rule=('optimizer success and certified quadrature ladder '
                         f'{multi.CERTIFICATION_ORDERS} with tolerances {joint_fit.QUADRATURE_TOLERANCES}'),
        selection_rule='lowest (4,128,128) chi2 among accepted runs with y_V free',
        adequacy_criterion=(f'total nominal p > {vpp_driver.ADEQUACY_TOTAL_P_VALUE} and no block with nominal p < '
                            f'{vpp_driver.ADEQUACY_BLOCK_P_VALUE} (row count as block dof), at (4,128,128)'),
        control_rule=(f'zero_control reproduces the baseline: scaled distance <= '
                      f'{comparison_module.SCALED_DISTANCE_TOLERANCE} and |chi2 - baseline| <= '
                      f'{comparison_module.CHI2_TOLERANCE}'),
        source_hash_gate=gate, runs=[],
        input_sha256={path: adequacy.file_hash(path) for path in
                      (adequacy.COMPARISON, sensitivity.ADEQUACY, vpp_driver.REFERENCE_LEDGER, *multi.INPUTS)},
        source_sha256=multi.source_hashes())


def evaluate_start(model, values, fixed_count):
    """One residual and Jacobian evaluation with finiteness and penalty checks; chi2 per block."""
    values = np.asarray(values, dtype=float)
    residual = model.residual(values)
    jacobian = model.jacobian(values)
    if (np.all(residual == 1e8+np.linalg.norm(values)) or not np.all(np.isfinite(residual))
            or not np.all(np.isfinite(jacobian))):
        raise ValueError('nonfinite or penalty residual or Jacobian at a start')
    free = COORDINATE_COUNT-fixed_count
    return residual, jacobian, dict(
        rows=int(residual.size), columns=int(jacobian.shape[1]), free_coordinate_count=free,
        degrees_of_freedom=int(residual.size)-free, chi2=float(residual @ residual),
        blocks={name: dict(chi2=block['chi2'], rows=block['count'])
                for name, block in adequacy.block_summary(residual).items()})


def dry_run(orders=FIT_ORDERS):
    """Everything up to the optimizer; returns the summary (nothing is written)."""
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    model = option.DoubleTraceSingletSourceLikelihood(*orders, bin_block_size=BLOCK_SIZE, **MODEL_OPTIONS)
    baseline = joint_fit.GlobalLikelihood(*orders, bin_block_size=BLOCK_SIZE)
    summary = dict(dry_run=True, optimizer_called=False, orders=list(orders), bin_block_size=BLOCK_SIZE,
                   optimizer_settings=dict(OPTIMIZER_SETTINGS), model=check_model(model),
                   source_hash_gate=vpp_driver.source_hash_gate(baseline), starts={})
    jacobians = {}
    for definition in start_definitions(context):
        values, fixed = start_from(model, context, definition)
        residual, jacobians[definition['name']], evaluation = evaluate_start(model, values, len(fixed))
        summary['starts'][definition['name']] = dict(evaluation, start_endpoint=definition['start_endpoint'],
                                                     coefficient_gev_inverse=definition['coefficient_gev_inverse'],
                                                     coefficient_fixed=definition['coefficient_fixed'])
        if definition['coefficient_gev_inverse'] == 0.:
            endpoint = np.asarray(endpoint_by_identifier(context, definition['start_endpoint']), dtype=float)
            summary['starts'][definition['name']]['equals_baseline_bitwise'] = dict(
                residual=bool(np.array_equal(residual, baseline.residual(endpoint))),
                jacobian_first_32_columns=bool(np.array_equal(jacobians[definition['name']][:, :BASELINE_COUNT],
                                                              baseline.jacobian(endpoint))))
            if not all(summary['starts'][definition['name']]['equals_baseline_bitwise'].values()):
                raise ValueError(f'start {definition["name"]} does not reproduce the baseline')
    scales = optimizer_scales(context, jacobians[CONTROL])
    summary.update(added_optimizer_scale_gev_inverse=float(scales[BASELINE_COUNT]),
                   baseline_chi2_high_order=context['baseline_chi2'], seconds=time.perf_counter()-started)
    return summary


def run_named(arguments):
    """Worker entry: run one start and save its part."""
    name, part_path, scales, settings = arguments
    context = sensitivity.accepted_context()
    definition = next(item for item in start_definitions(context) if item['name'] == name)
    result = run_start(definition, context, np.asarray(scales), settings,
                       checkpoint_path=part_path.removesuffix('.json')+'.progress.json')
    multi.write_ledger(part_path, result)
    return name


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', default=OUTPUT)
    parser.add_argument('--dry-run', action='store_true', help='check model, starts and source hashes; no optimizer')
    parser.add_argument('--workers', type=vpp_driver.positive_integer, default=2)
    args = parser.parse_args(argv)
    if args.dry_run:
        print(json.dumps(dry_run(), indent=1))
        return
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    low = make_models((FIT_ORDERS,))[0]
    check_model(low)
    gate = vpp_driver.source_hash_gate(joint_fit.GlobalLikelihood(*FIT_ORDERS, bin_block_size=BLOCK_SIZE))
    definitions = start_definitions(context)
    control = next(item for item in definitions if item['name'] == CONTROL)
    scales = optimizer_scales(context, low.jacobian(start_from(low, context, control)[0]))
    ledger = prepare_ledger(low, context, scales, definitions, gate)
    ledger['command'] = shlex.join([sys.executable, *sys.argv])
    multi.write_ledger(ROOT/args.output, ledger)
    parts = ROOT/(args.output+'.parts')
    parts.mkdir(parents=True, exist_ok=True)
    pending = [(item['name'], str(parts/f'{item["name"]}.json'), scales.tolist(), OPTIMIZER_SETTINGS)
               for item in definitions if not (parts/f'{item["name"]}.json').exists()]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for name in executor.map(run_named, pending):
            print('finished', name, flush=True)
    runs = [json.loads((parts/f'{item["name"]}.json').read_text()) for item in definitions]
    selected = selected_run(runs)
    ledger.update(runs=runs, selected_run=selected,
                  selected_adequacy=None if selected is None else runs[selected]['adequacy'],
                  control_reproduces_baseline=next(run.get('reproduces_baseline_endpoint') for run in runs
                                                   if run['name'] == CONTROL),
                  state='finished', driver_elapsed_seconds=time.perf_counter()-started)
    multi.write_ledger(ROOT/args.output, ledger)
    print(args.output)


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    main()
