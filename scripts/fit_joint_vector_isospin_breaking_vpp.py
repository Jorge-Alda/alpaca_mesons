"""Joint vector fit with the isospin-breaking VPP option switched on (prepared driver).

Dry run (no optimizer; builds the model, checks it, evaluates the residual and
Jacobian once at each start, checks the source-hash gate, prints a summary):
  .venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --dry-run
Fit (writes the ledger with bounds and starts before any optimizer call):
  .venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --output data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json

Model: IsospinBreakingVppLikelihood of joint_vector_isospin_breaking_vpp_option.py
with isospin_breaking_vpp=True and the default switches
vpp_kstar_excitation_scaling='ground_only', vpp_s_channel_factors=True. It has
35 coordinates: the 32 baseline coordinates and r1, r2, r3 in GeV^-2 with
bounds [-5, 5] GeV^-2. 75 residual rows; 75 - 35 = 40 nominal degrees of freedom.

Gates (same as scripts/refit_joint_vector_sensitivity.py, alp-uzh.6 ledgers):
  * scipy least_squares, OPTIMIZER_SETTINGS (trf, max_nfev=5000, ftol=1e-8,
    xtol=1e-9, gtol=1e-7), fixed optimizer scales, quadrature FIT_ORDERS
    (2,16,16) with four-bin blocks; optimizer status, message, nfev, njev,
    optimality, active bounds and near-singular directions recorded;
    failures kept in the ledger; per-run checkpoints allow resumption.
  * Re-evaluation at HIGH_ORDERS (4,128,128): chi2, dof, nominal p, block chi2.
  * Local errors sqrt(diag((J^T J)^-1)) of the free coordinates from the
    (4,128,128) Jacobian (assess_joint_vector_fit_adequacy.local_errors).
  * Quadrature ladder of the endpoint comparison at REFINEMENT_ORDERS
    (2,64,64) -> CERTIFICATION_ORDERS; accepted = optimizer success and
    certified ladder (compare_joint_vector_endpoints.audit_passed).
  * Multistart: every start is refit from the best accepted endpoint
    basin_probes/even_shared_widths; if its scaled distance to that endpoint
    exceeds 1 it is refit from the other three accepted endpoints with the
    same r start; the selected run has the lowest (4,128,128) chi2 among
    accepted runs (else among valid runs).
  * Source-hash gate: the six source ledgers of the endpoint comparison
    must match the hashes recorded in it and pass
    compare_joint_vector_endpoints.validate_ledgers; every input and source
    hash recorded in the alp-uzh.6 variation ledger REFERENCE_LEDGER must
    match the current files.
  * Adequacy criterion (openspec/changes/assess-joint-vector-fit-adequacy):
    total nominal p > 0.01 and no block with nominal p < 1e-3 (row count as
    block dof), at (4,128,128).

Starts (other 32 coordinates at the accepted endpoint basin_probes/even_shared_widths):
  scan_point: r = (-0.350, 1.855, 0.729) GeV^-2, r free.
  zero_control: r = (0, 0, 0) held fixed (choice: with r fixed the likelihood
      equals the baseline bitwise, so the refit must reproduce the baseline
      endpoint within the comparison tolerances: scaled distance <= 0.01 and
      |chi2(4,128,128) - baseline| <= 0.01).

Choices (stated, not results):
  * Optimizer scales: the 32 fixed scales of the endpoint comparison; for
    r1, r2, r3 the reciprocal Jacobian column norms at the zero_control start
    at the fit order (same definition as the baseline scales, which are
    reciprocal column norms at the stored start).
  * Scaled distance to the best endpoint (basin-switch rule) and parameter
    shifts in baseline local errors use the 32 baseline coordinates only; the
    best endpoint has no r coordinates. r values are reported in GeV^-2 with
    their own local errors.
  * The baseline quadrature audit labels Jacobian columns with
    joint_fit.PARAMETER_NAMES (32 names); during the audit and the
    near-singular analysis that module name is replaced by the 35 names
    (extended_parameter_names). FitLikelihood fixes its own names so the
    replacement does not change the likelihood.

Units: GeV, nb; r in GeV^-2; residuals, chi2 and scaled distances dimensionless.
Conventions of Alda:2026nrd.
"""
import argparse
import contextlib
import hashlib
import json
import os
import shlex
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_distribution

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assess_joint_vector_fit_adequacy as adequacy  # noqa: E402
import compare_joint_vector_endpoints as comparison_module  # noqa: E402
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_joint_vector_multistart as multi  # noqa: E402
import joint_vector_isospin_breaking_vpp_option as vpp_option  # noqa: E402
import refit_joint_vector_sensitivity as sensitivity  # noqa: E402

ROOT = joint_fit.ROOT
OUTPUT = 'data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json'
REFERENCE_LEDGER = 'data/vector_sector/joint_vector_sensitivity_variations.json'
BEST = sensitivity.BEST
FIT_ORDERS = sensitivity.FIT_ORDERS
HIGH_ORDERS = sensitivity.HIGH_ORDERS
REFINEMENT_ORDERS = sensitivity.REFINEMENT_ORDERS
BLOCK_SIZE = sensitivity.BLOCK_SIZE
OPTIMIZER_SETTINGS = sensitivity.OPTIMIZER_SETTINGS
BASIN_SWITCH_DISTANCE = sensitivity.BASIN_SWITCH_DISTANCE
MODEL_OPTIONS = dict(isospin_breaking_vpp=True, vpp_kstar_excitation_scaling='ground_only',
                     vpp_s_channel_factors=True)
BASELINE_COUNT = vpp_option.BASELINE_COORDINATE_COUNT
COORDINATE_COUNT = BASELINE_COUNT+len(vpp_option.ADDED_COORDINATE_NAMES)
RESIDUAL_COUNT = 75
ADEQUACY_TOTAL_P_VALUE = .01
ADEQUACY_BLOCK_P_VALUE = 1e-3
BASELINE_PARAMETER_NAMES = joint_fit.PARAMETER_NAMES
STARTS = (
    dict(name='scan_point', coefficients_gev_minus2=tuple(vpp_option.SCAN_POINT_GEV_MINUS2),
         coefficients_fixed=False,
         definition='r = (-0.350, 1.855, 0.729) GeV^-2 from the simultaneous scan, r free; '
                    'other coordinates at the accepted endpoint.'),
    dict(name='zero_control', coefficients_gev_minus2=(0., 0., 0.), coefficients_fixed=True,
         definition='r = (0, 0, 0) GeV^-2 held fixed; other coordinates at the accepted endpoint. '
                    'Must reproduce the baseline endpoint.'),
)


class FitLikelihood(vpp_option.IsospinBreakingVppLikelihood):
    """IsospinBreakingVppLikelihood whose coordinate names do not read joint_fit.PARAMETER_NAMES.

    Needed because extended_parameter_names temporarily replaces that module
    name by the 35 names during the quadrature audit.
    """

    @property
    def parameter_names(self):
        """The 32 baseline names followed by the r1, r2, r3 names."""
        return BASELINE_PARAMETER_NAMES+self.added_coordinate_names


@contextlib.contextmanager
def extended_parameter_names(names):
    """Let joint_fit.audit_quadrature and near_singular_columns label 35 columns."""
    original = joint_fit.PARAMETER_NAMES
    joint_fit.PARAMETER_NAMES = tuple(names)
    try:
        yield
    finally:
        joint_fit.PARAMETER_NAMES = original


def parse_arguments(argv=None):
    """Command-line options: --output, --dry-run, --workers."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', default=OUTPUT)
    parser.add_argument('--dry-run', action='store_true',
                        help='check model, starts and source hashes; no optimizer; print only')
    parser.add_argument('--workers', type=int, default=2)
    return parser.parse_args(argv)


def make_models(orders_list):
    """Option likelihoods (option on, default switches) at each order, four-bin blocks."""
    return [FitLikelihood(*orders, bin_block_size=BLOCK_SIZE, **MODEL_OPTIONS) for orders in orders_list]


def check_model(model):
    """Check option, switches, coordinate count and bounds; return a summary."""
    if not (model.isospin_breaking_vpp and model.s_channel and model.excitation_scaling == 'ground_only'):
        raise ValueError('isospin-breaking VPP option with default switches is required')
    names = model.parameter_names
    if len(names) != COORDINATE_COUNT or tuple(names[BASELINE_COUNT:]) != vpp_option.ADDED_COORDINATE_NAMES:
        raise ValueError(f'expected {COORDINATE_COUNT} coordinates')
    lower, upper = model.parameter_bounds()
    baseline_lower, baseline_upper = joint_fit.parameter_bounds(model)
    bound = vpp_option.COEFFICIENT_BOUND_GEV_MINUS2
    if not (np.array_equal(lower[:BASELINE_COUNT], baseline_lower)
            and np.array_equal(upper[:BASELINE_COUNT], baseline_upper)
            and np.all(lower[BASELINE_COUNT:] == -bound) and np.all(upper[BASELINE_COUNT:] == bound)):
        raise ValueError('bounds differ from the baseline bounds plus [-5, 5] GeV^-2')
    return dict(isospin_breaking_vpp=True, vpp_s_channel_factors=model.s_channel,
                vpp_kstar_excitation_scaling=model.excitation_scaling,
                coordinate_count=len(names), residual_count=RESIDUAL_COUNT,
                nominal_degrees_of_freedom=RESIDUAL_COUNT-len(names),
                added_coordinate_bounds_gev_minus2={name: [float(lower[index]), float(upper[index])]
                                                    for index, name in enumerate(names)
                                                    if index >= BASELINE_COUNT})


def validate_point(model, values, lower, upper):
    """35 finite coordinates inside the bounds, on the likelihood's physical branch."""
    values = np.asarray(values, dtype=float)
    if values.shape != (COORDINATE_COUNT,) or not np.all(np.isfinite(values)):
        raise ValueError(f'start must have {COORDINATE_COUNT} finite coordinates')
    if np.any(values < lower) or np.any(values > upper):
        raise ValueError('start is outside parameter bounds')
    model.unpack(values)
    return values


def start_from(model, endpoint, definition):
    """35-coordinate start: an accepted 32-coordinate endpoint extended by the start's r."""
    values = model.extend(endpoint, definition['coefficients_gev_minus2'])
    fixed = ({BASELINE_COUNT+index: float(value)
              for index, value in enumerate(definition['coefficients_gev_minus2'])}
             if definition['coefficients_fixed'] else {})
    validate_point(model, values, *model.parameter_bounds())
    return values, fixed


def start_points(model, context):
    """Every start of STARTS from the best accepted endpoint."""
    best = context['accepted'][0]
    if best['id'] != BEST:
        raise ValueError('best accepted endpoint changed')
    result = []
    for definition in STARTS:
        values, fixed = start_from(model, best['endpoint'], definition)
        result.append(dict(definition, start_endpoint=best['id'], start=values.tolist(),
                           fixed_coordinates={model.parameter_names[index]: value
                                              for index, value in fixed.items()},
                           free_coordinate_count=COORDINATE_COUNT-len(fixed)))
    return result


def optimizer_scales(context, control_jacobian):
    """Comparison-ledger scales plus reciprocal r column norms at the zero_control start."""
    norms = np.linalg.norm(np.asarray(control_jacobian)[:, BASELINE_COUNT:], axis=0)
    if not np.all(np.isfinite(norms)) or np.any(norms <= 0):
        raise ValueError('zero or nonfinite r column at the zero_control start')
    return np.r_[np.asarray(context['scales'], dtype=float), 1/norms]


def file_hash(path):
    """SHA-256 of a file relative to the repository root."""
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def source_hash_gate(baseline_model, reference_ledger=REFERENCE_LEDGER):
    """Raise ValueError unless the accepted endpoints and alp-uzh.6 sources are unchanged."""
    comparison = json.loads((ROOT/adequacy.COMPARISON).read_bytes())
    paths = {role: ROOT/entry['path'] for role, entry in comparison['sources'].items()}
    ledgers, hashes = comparison_module.read_ledgers(paths)
    mismatches = [role for role, entry in comparison['sources'].items()
                  if hashes[role]['sha256'] != entry['sha256']]
    if mismatches:
        raise ValueError(f'source ledgers differ from the endpoint comparison: {mismatches}')
    comparison_module.validate_ledgers(ledgers, hashes, baseline_model)
    reference = json.loads((ROOT/reference_ledger).read_bytes())
    recorded = {**reference['input_sha256'], **reference['source_sha256']}
    changed = [path for path, value in recorded.items()
               if not (ROOT/path).exists() or file_hash(path) != value]
    if changed:
        raise ValueError(f'files changed since {reference_ledger}: {changed}')
    return dict(passed=True, comparison_source_ledgers=sorted(hashes),
                reference_ledger=reference_ledger, reference_ledger_sha256=file_hash(reference_ledger),
                reference_files_checked=len(recorded))


def refit(model, start, fixed, scales, settings, checkpoint_path=None):
    """35-coordinate version of refit_joint_vector_sensitivity.refit (same record keys)."""
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
        with extended_parameter_names(names):
            record['jacobian_conditioning'] = sensitivity.near_singular_columns(result.jac, scales, free)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        record.update(state='error', success=False, error_type=type(error).__name__, message=str(error))
    record['refit_seconds'] = record['prior_seconds']+time.perf_counter()-started
    return record


def local_error_record(model, run):
    """Local errors of the free coordinates from the high-order Jacobian at the endpoint."""
    names = model.parameter_names
    free = [index for index, name in enumerate(names) if name not in run['fixed_coordinates']]
    try:
        jacobian = model.jacobian(np.asarray(run['endpoint']))[:, free]
        if not np.all(np.isfinite(jacobian)):
            raise ValueError('nonfinite high-order Jacobian')
        errors, singular = adequacy.local_errors(jacobian)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        return dict(valid=False, error_type=type(error).__name__, message=str(error))
    return dict(valid=True, local_errors={names[index]: float(value) for index, value in zip(free, errors)},
                column_normalized_singular_values=singular.tolist(),
                definition='sqrt(diag((J^T J)^-1)) of the free coordinates from the (4,128,128) Jacobian; '
                           'linearized, not profile intervals')


def shifts(model, endpoint, context):
    """Baseline-coordinate shifts in baseline local errors; r values in GeV^-2."""
    converted = sensitivity.baseline_coordinates(model, endpoint)
    best = np.asarray(context['accepted'][0]['endpoint'])
    shift = (converted[:BASELINE_COUNT]-best)/np.asarray(context['local_errors'])
    distance = float(np.linalg.norm((converted[:BASELINE_COUNT]-best)/np.asarray(context['scales'])))
    largest = int(np.argmax(abs(shift)))
    return dict(scaled_distance_to_best=distance,
                shifts_in_local_errors=dict(zip(BASELINE_PARAMETER_NAMES, shift.tolist())),
                maximum_absolute_shift=float(abs(shift[largest])),
                maximum_shift_parameter=BASELINE_PARAMETER_NAMES[largest],
                coefficients_gev_minus2=dict(zip(vpp_option.ADDED_COORDINATE_NAMES,
                                                 np.asarray(endpoint)[BASELINE_COUNT:].tolist())))


def audit(model, run):
    """Quadrature ladder of the endpoint comparison on 35 coordinates; acceptance flag."""
    if not run.get('success') or run.get('state') != 'finished' or run.get('penalty_endpoint'):
        result = dict(status='not_evaluated', passed=False,
                      reason='optimizer did not succeed at a valid endpoint')
        return result, False
    try:
        values = validate_point(model, run['endpoint'], *model.parameter_bounds())
        started = time.perf_counter()
        with extended_parameter_names(model.parameter_names):
            result = joint_fit.audit_quadrature(model, values, multi.CERTIFICATION_ORDERS,
                                                parameter_scales=np.ones(COORDINATE_COUNT))
        result['elapsed_seconds'] = time.perf_counter()-started
        result['bin_block_size'] = model.bin_block_size
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        result = dict(status='unevaluated', passed=False, error_type=type(error).__name__, reason=str(error))
    return result, bool(run.get('success') and comparison_module.audit_passed(result))


def adequacy_verdict(high_order):
    """Total nominal p > 0.01 and no block with nominal p < 1e-3 (row count as block dof)."""
    if not high_order.get('valid') or high_order.get('nominal_p_value') is None:
        return dict(evaluated=False, meets_criterion=False)
    blocks = {name: block['nominal_p_value_row_count_dof'] for name, block in high_order['blocks'].items()
              if block['nominal_p_value_row_count_dof'] is not None}
    failing = sorted(name for name, value in blocks.items() if value < ADEQUACY_BLOCK_P_VALUE)
    total = high_order['nominal_p_value']
    return dict(evaluated=True, total_nominal_p_value=total, failing_blocks=failing,
                meets_criterion=bool(total > ADEQUACY_TOTAL_P_VALUE and not failing))


def control_reproduces_baseline(run, context):
    """Scaled distance <= 0.01 and |chi2 - baseline| <= 0.01 at (4,128,128)."""
    if not run.get('high_order', {}).get('valid'):
        return False
    return bool(run['scaled_distance_to_best'] <= comparison_module.SCALED_DISTANCE_TOLERANCE
                and abs(run['high_order']['chi2']-context['baseline_chi2']) <= comparison_module.CHI2_TOLERANCE)


def run_from(models, definition, start_endpoint, context, scales, settings, checkpoint_path=None):
    """One refit of a start from one accepted endpoint, with all gates."""
    low, high, refinement = models
    values, fixed = start_from(low, start_endpoint['endpoint'], definition)
    run = dict(start_endpoint=start_endpoint['id'])
    run.update(refit(low, values, fixed, scales, settings, checkpoint_path))
    if run['state'] != 'finished':
        run['accepted'] = False
        return run
    run['high_order'] = sensitivity.high_order_evaluation(high, run, run['free_coordinate_count'],
                                                          adequacy.row_labels(high))
    sensitivity.add_constrained_degrees_of_freedom(run, run['high_order'])
    run['local_errors'] = local_error_record(high, run)
    run.update(shifts(low, run['endpoint'], context))
    started = time.perf_counter()
    run['quadrature_refinement'], run['accepted'] = audit(refinement, run)
    run['audit_seconds'] = time.perf_counter()-started
    run['adequacy'] = adequacy_verdict(run['high_order'])
    return run


def run_start(definition, context, scales, models=None, settings=OPTIMIZER_SETTINGS, checkpoint_prefix=None):
    """Refit a start from the best endpoint; from the other accepted endpoints if it moved far."""
    started = time.perf_counter()
    models = models or make_models((FIT_ORDERS, HIGH_ORDERS, REFINEMENT_ORDERS))

    def checkpoint(index):
        return checkpoint_prefix and f'{checkpoint_prefix}.{index}.progress.json'
    runs = [run_from(models, definition, context['accepted'][0], context, scales, settings, checkpoint(0))]
    if runs[0].get('scaled_distance_to_best', 0.) > BASIN_SWITCH_DISTANCE:
        runs += [run_from(models, definition, endpoint, context, scales, settings, checkpoint(index))
                 for index, endpoint in enumerate(context['accepted'][1:], start=1)]
    record = sensitivity.variation_record(definition, runs, time.perf_counter()-started)
    selected = runs[record['selected_run']] if 'selected_run' in record else {}
    record['adequacy'] = adequacy_verdict(selected.get('high_order', {}))
    if definition['coefficients_fixed']:
        record['reproduces_baseline_endpoint'] = control_reproduces_baseline(selected, context)
    return record


def prepare_ledger(model, context, scales, starts, gate):
    """Ledger written before any optimizer call: model, bounds, starts, scales and gates."""
    lower, upper = model.parameter_bounds()
    quadrature = dict(zip(('energy', 'pair_per_segment', 'angle_per_segment'), FIT_ORDERS), bin_block_size=BLOCK_SIZE)
    return dict(
        study='isospin_breaking_vpp', state='prepared', model_options=MODEL_OPTIONS,
        best_accepted_endpoint=BEST, accepted_endpoints=[item['id'] for item in context['accepted']],
        baseline_chi2_high_order=context['baseline_chi2'], baseline_blocks=context['baseline_blocks'],
        parameter_names=list(model.parameter_names),
        bounds=dict(lower=multi.json_bounds(lower), upper=multi.json_bounds(upper)),
        starts=starts, residual_count=RESIDUAL_COUNT,
        nominal_degrees_of_freedom=RESIDUAL_COUNT-COORDINATE_COUNT,
        fit_quadrature=quadrature,
        evaluation_quadrature=dict(zip(quadrature, HIGH_ORDERS), bin_block_size=BLOCK_SIZE),
        refinement_quadrature=dict(zip(quadrature, REFINEMENT_ORDERS), bin_block_size=BLOCK_SIZE,
                                   certification_orders=multi.CERTIFICATION_ORDERS),
        optimizer_settings=OPTIMIZER_SETTINGS, optimizer_coordinate_scales=list(map(float, scales)),
        added_coordinate_scale_definition='reciprocal Jacobian column norms of r1, r2, r3 at the zero_control '
                                          'start, fit order (assumption); first 32 from the endpoint comparison',
        degrees_of_freedom_counting='retained residual rows minus free coordinates; '
                                    'per-block p-values use the row count (assumption)',
        shift_definition='(baseline coordinates converted to the baseline reference point - best endpoint) / '
                         'baseline (4,128,128) local errors; r reported in GeV^-2',
        near_singular_definition='SVD of the fit-order Jacobian times optimizer scales; relative singular '
                                 f'value below {sensitivity.SINGULAR_RELATIVE_THRESHOLD}; run then labelled conditional',
        acceptance_rule=('optimizer success and certified quadrature ladder '
                         f'{multi.CERTIFICATION_ORDERS} with tolerances {joint_fit.QUADRATURE_TOLERANCES}'),
        basin_switch_rule=(f'scaled distance (32 baseline coordinates) > {BASIN_SWITCH_DISTANCE}: refit from the '
                           'other three accepted endpoints; selected run = lowest high-order chi2 among accepted runs'),
        adequacy_criterion=(f'total nominal p > {ADEQUACY_TOTAL_P_VALUE} and no block with nominal p < '
                            f'{ADEQUACY_BLOCK_P_VALUE} (row count as block dof), at (4,128,128)'),
        control_rule=(f'zero_control reproduces the baseline: scaled distance <= '
                      f'{comparison_module.SCALED_DISTANCE_TOLERANCE} and |chi2 - baseline| <= '
                      f'{comparison_module.CHI2_TOLERANCE}'),
        source_hash_gate=gate, entries=[],
        input_sha256={path: adequacy.file_hash(path) for path in
                      (adequacy.COMPARISON, sensitivity.ADEQUACY, REFERENCE_LEDGER, *multi.INPUTS)},
        source_sha256=multi.source_hashes())


def evaluate_start(model, values, fixed_count):
    """One residual and Jacobian evaluation with finiteness and penalty checks."""
    values = np.asarray(values, dtype=float)
    residual = model.residual(values)
    jacobian = model.jacobian(values)
    penalty = bool(np.all(residual == 1e8+np.linalg.norm(values)))
    if penalty or not np.all(np.isfinite(residual)) or not np.all(np.isfinite(jacobian)):
        raise ValueError('nonfinite or penalty residual or Jacobian at a start')
    free = COORDINATE_COUNT-fixed_count
    return residual, jacobian, dict(
        rows=int(residual.size), columns=int(jacobian.shape[1]), free_coordinate_count=free,
        degrees_of_freedom=int(residual.size)-free, chi2=float(residual @ residual),
        blocks={name: dict(chi2=block['chi2'], rows=block['count'])
                for name, block in adequacy.block_summary(residual).items()},
        residual_finite=True, jacobian_finite=True)


def dry_run(orders=FIT_ORDERS):
    """Everything up to the optimizer; returns the summary (nothing is written)."""
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    model = FitLikelihood(*orders, bin_block_size=BLOCK_SIZE, **MODEL_OPTIONS)
    baseline = joint_fit.GlobalLikelihood(*orders, bin_block_size=BLOCK_SIZE)
    summary = dict(dry_run=True, optimizer_called=False, orders=list(orders), bin_block_size=BLOCK_SIZE,
                   model=check_model(model), source_hash_gate=source_hash_gate(baseline))
    starts = start_points(model, context)
    best = np.asarray(context['accepted'][0]['endpoint'])
    evaluations, jacobians, residuals = {}, {}, {}
    for item in starts:
        residuals[item['name']], jacobians[item['name']], evaluations[item['name']] = evaluate_start(
            model, item['start'], len(item['fixed_coordinates']))
        evaluations[item['name']]['start_coefficients_gev_minus2'] = list(item['coefficients_gev_minus2'])
        evaluations[item['name']]['coefficients_fixed'] = item['coefficients_fixed']
    control_residual = residuals['zero_control']
    baseline_residual = baseline.residual(best)
    baseline_jacobian = baseline.jacobian(best)
    summary['starts'] = evaluations
    summary['zero_control_equals_baseline'] = dict(
        residual_bitwise=bool(np.array_equal(control_residual, baseline_residual)),
        jacobian_first_32_columns_bitwise=bool(np.array_equal(
            jacobians['zero_control'][:, :BASELINE_COUNT], baseline_jacobian)))
    if not all(summary['zero_control_equals_baseline'].values()):
        raise ValueError('zero_control start does not reproduce the baseline residual and Jacobian')
    scales = optimizer_scales(context, jacobians['zero_control'])
    summary['added_optimizer_scales'] = dict(zip(vpp_option.ADDED_COORDINATE_NAMES,
                                                 scales[BASELINE_COUNT:].tolist()))
    summary['baseline_chi2_high_order'] = context['baseline_chi2']
    summary['seconds'] = time.perf_counter()-started
    return summary


def run_named(arguments):
    """Worker entry: run one start and save its part."""
    name, part_path, scales = arguments
    definition = next(item for item in STARTS if item['name'] == name)
    result = run_start(definition, sensitivity.accepted_context(), np.asarray(scales),
                       checkpoint_prefix=part_path.removesuffix('.json'))
    multi.write_ledger(part_path, result)
    return name


def main(argv=None):
    args = parse_arguments(argv)
    if args.dry_run:
        print(json.dumps(dry_run(), indent=1))
        return
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    low = make_models((FIT_ORDERS,))[0]
    check_model(low)
    gate = source_hash_gate(joint_fit.GlobalLikelihood(*FIT_ORDERS, bin_block_size=BLOCK_SIZE))
    starts = start_points(low, context)
    control = next(item for item in starts if item['name'] == 'zero_control')
    scales = optimizer_scales(context, low.jacobian(np.asarray(control['start'])))
    ledger = prepare_ledger(low, context, scales, starts, gate)
    ledger['command'] = shlex.join([sys.executable, *sys.argv])
    multi.write_ledger(ROOT/args.output, ledger)
    parts = ROOT/(args.output+'.parts')
    parts.mkdir(parents=True, exist_ok=True)
    pending = [(item['name'], str(parts/f'{item["name"]}.json'), scales.tolist()) for item in STARTS
               if not (parts/f'{item["name"]}.json').exists()]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for name in executor.map(run_named, pending):
            print('finished', name, flush=True)
    ledger['entries'] = [json.loads((parts/f'{item["name"]}.json').read_text()) for item in STARTS]
    ledger.update(state='finished', driver_elapsed_seconds=time.perf_counter()-started)
    multi.write_ledger(ROOT/args.output, ledger)
    print(args.output)


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    main()
