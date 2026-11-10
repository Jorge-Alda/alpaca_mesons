"""Joint vector fit with the isospin-breaking VPP option switched on (prepared driver).

Dry run (no optimizer; builds the model, checks it, evaluates the residual and
Jacobian once at each start, checks the source-hash gate, prints a summary):
  .venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --dry-run
Fit (writes the ledger with bounds and starts before any optimizer call):
  .venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --output data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json
Restart from the lowest-(4,128,128)-chi2 run of entry 'scan_point' of an existing
ledger (single start, r free, same bounds; dry run shown, drop --dry-run to fit):
  .venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --dry-run --restart-from data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json
Option --max-nfev N (default 5000) sets the optimizer evaluation limit and is
recorded in the ledger. Options --ftol, --xtol, --gtol (defaults 1e-8, 1e-9,
1e-7, the OPTIMIZER_SETTINGS values) set the least_squares tolerances and are
recorded in the ledger's optimizer_settings.
Profile in r (needs --restart-from; dry run shown, drop --dry-run to run the
refits; several --profile-r flags give one 1-D scan per r):
  .venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --dry-run --restart-from data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json --profile-r 1=-1,-0.5,0,0.5,1

Model: IsospinBreakingVppLikelihood of joint_vector_isospin_breaking_vpp_option.py,
evaluated through IsospinBreakingVppLikelihoodWithReuse
(joint_vector_isospin_breaking_vpp_reuse.py; bitwise-equal residual and Jacobian,
fewer repeated calculations), with isospin_breaking_vpp=True and the default switches
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

Restart (--restart-from PATH): the 35-coordinate endpoint of the run with the
lowest (4,128,128) chi2 of entry 'scan_point' in the ledger PATH (accepted or
not) is the single start, with r free and the bounds above. The source ledger
path, its SHA-256 and the run index are stored in the start definition. The
driver refuses to start if that endpoint is outside the bounds or gives a
nonfinite or penalty residual at the fit order. No basin switch and no
zero_control run (choice: a restart continues one run; the basin-switch rule
belongs to the original multistart). Optimizer scales: same definition as above,
recomputed at the zero_control start. Default output in restart mode:
RESTART_OUTPUT (never the source ledger).

Profile (--profile-r INDEX=V1,V2,..., INDEX in 1, 2, 3 for r1, r2, r3, values
in GeV^-2 inside [-5, 5]; with --restart-from): for every grid point the chosen
r_i is held fixed at the grid value (same fixed-coordinate mechanism as
zero_control) and the 32 baseline coordinates and the other r are refit from the
restart endpoint with the fixed value substituted. Each --profile-r flag is a
separate 1-D scan; a point fixing several r at once exists only with
--profile-grid (Cartesian product of the flags). Per point: refit at
FIT_ORDERS, then the HIGH_ORDERS chi2 per block. Choices: no basin switch, no
zero_control run, no quadrature ladder (a profile point is a constrained refit,
not a new candidate endpoint); optimizer scales as in restart mode. Ledger
study PROFILE_STUDY with a 'profile' list; default output PROFILE_OUTPUT; an
existing output is refused unless --overwrite (finished points in
OUTPUT.parts/ are reused when their start and optimizer settings match).
With --dry-run: one residual and Jacobian at the first point of each 1-D scan
(or the first product point with --profile-grid), point count and a time
estimate under the assumption of PROFILE_MINUTES_PER_POINT minutes per point;
nothing is written.

Units: GeV, nb; r in GeV^-2; residuals, chi2 and scaled distances dimensionless.
Conventions of Alda:2026nrd.
"""
import argparse
import contextlib
import hashlib
import itertools
import json
import math
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
import joint_vector_isospin_breaking_vpp_reuse as vpp_reuse  # noqa: E402
import refit_joint_vector_sensitivity as sensitivity  # noqa: E402

ROOT = joint_fit.ROOT
OUTPUT = 'data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json'
RESTART_OUTPUT = 'data/vector_sector/joint_vector_isospin_breaking_vpp_restart_fit.json'
RESTART_ENTRY = 'scan_point'
PROFILE_OUTPUT = 'data/vector_sector/joint_vector_isospin_breaking_vpp_r_profile.json'
PROFILE_STUDY = 'isospin_breaking_vpp_r_profile'
# Assumption for the dry-run time estimate: wall time of one profile refit,
# the measured value supplied with the task (not measured by this driver).
PROFILE_MINUTES_PER_POINT = 9.
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


class FitLikelihood(vpp_reuse.IsospinBreakingVppLikelihoodWithReuse):
    """Option likelihood whose coordinate names do not read joint_fit.PARAMETER_NAMES.

    Needed because extended_parameter_names temporarily replaces that module
    name by the 35 names during the quadrature audit. The parent class reuses
    repeated calculations of IsospinBreakingVppLikelihood with bitwise-equal
    residuals and Jacobians (joint_vector_isospin_breaking_vpp_reuse.py and its
    tests); the original scan_point/zero_control fit ledger was produced with
    IsospinBreakingVppLikelihood itself.
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


def positive_integer(text):
    """argparse type: an integer >= 1."""
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError('must be a positive integer')
    return value


def positive_float(text):
    """argparse type: a finite float > 0."""
    value = float(text)
    if not (math.isfinite(value) and value > 0):
        raise argparse.ArgumentTypeError('must be a finite positive number')
    return value


def profile_specification(text):
    """argparse type: 'INDEX=V1,V2,...' -> (INDEX, (V1, V2, ...)); r_INDEX grid in GeV^-2.

    INDEX is 1, 2 or 3; values are finite, distinct and inside the r bounds
    [-5, 5] GeV^-2.
    """
    index_text, separator, values_text = text.partition('=')
    if not separator or index_text.strip() not in ('1', '2', '3'):
        raise argparse.ArgumentTypeError(f'{text!r}: expected INDEX=V1,V2,... with INDEX 1, 2 or 3')
    try:
        values = tuple(float(item) for item in values_text.split(','))
    except ValueError:
        raise argparse.ArgumentTypeError(f'{text!r}: grid values must be numbers') from None
    bound = vpp_option.COEFFICIENT_BOUND_GEV_MINUS2
    if not all(math.isfinite(value) and -bound <= value <= bound for value in values):
        raise argparse.ArgumentTypeError(f'{text!r}: grid values must lie in [-{bound:g}, {bound:g}] GeV^-2')
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError(f'{text!r}: repeated grid value')
    return int(index_text), values


def parse_arguments(argv=None):
    """Command-line options of the driver (see the module docstring).

    The default output is OUTPUT, RESTART_OUTPUT with --restart-from, or
    PROFILE_OUTPUT with --profile-r. --profile-r needs --restart-from and
    distinct INDEX values; --profile-grid needs --profile-r.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', default=None)
    parser.add_argument('--dry-run', action='store_true',
                        help='check model, starts and source hashes; no optimizer; print only')
    parser.add_argument('--workers', type=positive_integer, default=2)
    parser.add_argument('--restart-from', default=None, metavar='PATH',
                        help=f"ledger whose lowest high-order chi2 run of entry '{RESTART_ENTRY}' is the single start")
    parser.add_argument('--max-nfev', type=positive_integer, default=OPTIMIZER_SETTINGS['max_nfev'],
                        help='optimizer max_nfev (default %(default)s)')
    for name in ('ftol', 'xtol', 'gtol'):
        parser.add_argument(f'--{name}', type=positive_float, default=OPTIMIZER_SETTINGS[name],
                            help=f'least_squares {name} (default %(default)s)')
    parser.add_argument('--profile-r', type=profile_specification, action='append', default=None,
                        metavar='INDEX=V1,V2,...',
                        help='1-D profile of r_INDEX (1, 2 or 3) on the grid V1,V2,... GeV^-2; repeatable')
    parser.add_argument('--profile-grid', action='store_true',
                        help='profile on the Cartesian product of the --profile-r grids (several r fixed at once)')
    parser.add_argument('--overwrite', action='store_true',
                        help='profile only: allow replacing an existing output ledger')
    args = parser.parse_args(argv)
    if args.profile_r:
        if not args.restart_from:
            parser.error('--profile-r needs --restart-from (the profile starts from the restart endpoint)')
        indices = [index for index, _ in args.profile_r]
        if len(set(indices)) != len(indices):
            parser.error('--profile-r: each INDEX may be given once')
    elif args.profile_grid:
        parser.error('--profile-grid needs --profile-r')
    if args.output is None:
        args.output = PROFILE_OUTPUT if args.profile_r else RESTART_OUTPUT if args.restart_from else OUTPUT
    if args.restart_from and (ROOT/args.output).resolve() == (ROOT/args.restart_from).resolve():
        parser.error('--output must differ from the --restart-from ledger')
    return args


def optimizer_settings(max_nfev, ftol=OPTIMIZER_SETTINGS['ftol'], xtol=OPTIMIZER_SETTINGS['xtol'],
                       gtol=OPTIMIZER_SETTINGS['gtol']):
    """OPTIMIZER_SETTINGS with max_nfev and the three tolerances replaced (key order kept)."""
    return dict(OPTIMIZER_SETTINGS, max_nfev=int(max_nfev), ftol=float(ftol), xtol=float(xtol), gtol=float(gtol))


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


def restart_run(ledger, entry_name=RESTART_ENTRY):
    """(index, run) with the lowest valid (4,128,128) chi2 of an entry, accepted or not."""
    entry = next((item for item in ledger['entries'] if item['name'] == entry_name), None)
    if entry is None:
        raise ValueError(f'no entry {entry_name!r} in the restart ledger')
    candidates = [(index, run) for index, run in enumerate(entry['runs'])
                  if run.get('state') == 'finished' and run.get('high_order', {}).get('valid')
                  and 'endpoint' in run]
    if not candidates:
        raise ValueError(f'no finished run with a valid high-order chi2 in entry {entry_name!r}')
    return min(candidates, key=lambda item: item[1]['high_order']['chi2'])


def restart_definition(model, ledger_path, entry_name=RESTART_ENTRY):
    """Single restart start from a ledger; refuses endpoints outside bounds or with bad residuals.

    Returns (definition, start endpoint record) in the form used by run_from:
    the record holds the 32 baseline coordinates of the endpoint, the
    definition its r (free), so start_from rebuilds the 35-coordinate endpoint.
    """
    path = ROOT/ledger_path
    ledger = json.loads(path.read_bytes())
    if ledger.get('parameter_names') != list(model.parameter_names):
        raise ValueError('restart ledger has different coordinate names')
    index, run = restart_run(ledger, entry_name)
    endpoint = np.asarray(run['endpoint'], dtype=float)
    validate_point(model, endpoint, *model.parameter_bounds())
    residual = model.residual(endpoint)
    if not np.all(np.isfinite(residual)) or np.all(residual == 1e8+np.linalg.norm(endpoint)):
        raise ValueError('nonfinite or penalty residual at the restart endpoint')
    identifier = f'restart:{ledger_path}#{entry_name}.run{index}'
    definition = dict(
        name=f'{entry_name}_restart', coefficients_gev_minus2=tuple(endpoint[BASELINE_COUNT:].tolist()),
        coefficients_fixed=False,
        definition=(f"restart from run {index} of entry '{entry_name}' of {ledger_path} (lowest (4,128,128) "
                    'chi2, accepted or not): its 35-coordinate endpoint, r free, same bounds; single start, '
                    'no basin switch.'),
        restart_source=dict(ledger=str(ledger_path), ledger_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            entry=entry_name, run_index=index, run_start_endpoint=run.get('start_endpoint'),
                            run_high_order_chi2=run['high_order']['chi2'],
                            run_fit_order_chi2=run.get('fit_order_chi2'), run_accepted=run.get('accepted'),
                            run_message=run.get('message'), run_cumulative_nfev=run.get('cumulative_nfev')))
    return definition, dict(id=identifier, endpoint=endpoint[:BASELINE_COUNT].tolist())


def restart_start_record(model, definition, start_endpoint):
    """Ledger start item of a restart definition (same keys as start_points)."""
    values, fixed = start_from(model, start_endpoint['endpoint'], definition)
    return dict(definition, start_endpoint=start_endpoint['id'], start=values.tolist(),
                fixed_coordinates={model.parameter_names[index]: value for index, value in fixed.items()},
                free_coordinate_count=COORDINATE_COUNT-len(fixed))


def profile_points(specifications, cartesian_product=False):
    """Fixed r values of every profile point, as {coordinate index: value in GeV^-2}.

    specifications: (INDEX, grid values) pairs of --profile-r, INDEX 1, 2, 3
    for r1, r2, r3 (coordinate index BASELINE_COUNT+INDEX-1). Without
    cartesian_product the 1-D scans follow each other, one r per point; with it
    the points are the Cartesian product of the grids.
    """
    if cartesian_product:
        indices = [BASELINE_COUNT+index-1 for index, _ in specifications]
        return [dict(zip(indices, values)) for values in itertools.product(*(grid for _, grid in specifications))]
    return [{BASELINE_COUNT+index-1: value} for index, grid in specifications for value in grid]


def first_profile_points(specifications, cartesian_product=False):
    """First point of each 1-D scan, or the first Cartesian-product point."""
    if cartesian_product:
        return profile_points(specifications, True)[:1]
    return [{BASELINE_COUNT+index-1: grid[0]} for index, grid in specifications]


def profile_start(start, fixed, cartesian_product=False):
    """Copy of the 35-coordinate start with the fixed r values substituted.

    More than one fixed r is a Cartesian-product point and is refused unless
    cartesian_product is set.
    """
    if len(fixed) > 1 and not cartesian_product:
        raise ValueError('a profile point fixing more than one r needs --profile-grid')
    values = np.array(start, dtype=float)
    for index, value in fixed.items():
        values[index] = value
    return values


def profile_part_name(model, fixed):
    """Part-file name of a profile point, from the fixed coordinate names and values."""
    return '__'.join(f'{model.parameter_names[index]}={value!r}' for index, value in sorted(fixed.items()))+'.json'


def refuse_existing_output(path, overwrite):
    """Raise FileExistsError if the output ledger exists and overwrite is not set."""
    if Path(path).exists() and not overwrite:
        raise FileExistsError(f'{path} exists; give --overwrite to replace it')


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


def run_restart(definition, start_endpoint, context, scales, models=None, settings=OPTIMIZER_SETTINGS,
                checkpoint_prefix=None):
    """Single refit from a restart start (no basin switch), with all gates of run_from."""
    started = time.perf_counter()
    models = models or make_models((FIT_ORDERS, HIGH_ORDERS, REFINEMENT_ORDERS))
    checkpoint = checkpoint_prefix and f'{checkpoint_prefix}.0.progress.json'
    runs = [run_from(models, definition, start_endpoint, context, scales, settings, checkpoint)]
    record = sensitivity.variation_record(definition, runs, time.perf_counter()-started)
    selected = runs[record['selected_run']] if 'selected_run' in record else {}
    record['adequacy'] = adequacy_verdict(selected.get('high_order', {}))
    return record


def run_profile_point(models, start, fixed, scales, settings, checkpoint_path=None):
    """One profile refit: fixed r held at the grid values, then the (4,128,128) chi2 per block.

    models: (fit-order model, high-order model). No basin switch, no
    zero_control run and no quadrature ladder (choices, see the module docstring).
    """
    low, high = models
    names = low.parameter_names
    run = refit(low, start, fixed, scales, settings, checkpoint_path)
    record = dict(profile=True, fixed_values_gev_minus2={names[index]: float(value) for index, value in fixed.items()},
                  optimizer_settings=dict(settings), state=run['state'], status=run.get('status'),
                  success=run['success'], message=run.get('message'), nfev=run.get('nfev'), njev=run.get('njev'),
                  optimality=run.get('optimality'), fit_order_chi2=run.get('fit_order_chi2'),
                  coordinates=run.get('endpoint'))
    if run['state'] == 'finished':
        high_order = sensitivity.high_order_evaluation(high, run, run['free_coordinate_count'],
                                                       adequacy.row_labels(high))
        sensitivity.add_constrained_degrees_of_freedom(run, high_order)
        record.update(high_order_chi2=high_order.get('chi2'),
                      high_order_block_chi2={name: block['chi2'] for name, block in high_order.get('blocks', {}).items()},
                      high_order=high_order, adequacy=adequacy_verdict(high_order))
    record['refit'] = run
    return record


def profile_part_matches(path, start, settings):
    """True if a saved profile part exists with the same start and optimizer settings."""
    if not Path(path).exists():
        return False
    saved = json.loads(Path(path).read_text())
    return saved['refit']['start'] == list(map(float, start)) and saved['optimizer_settings'] == dict(settings)


def run_profile_named(arguments):
    """Worker entry: one profile point, saved as its own part."""
    part_path, start, fixed, scales, settings = arguments
    record = run_profile_point(make_models((FIT_ORDERS, HIGH_ORDERS)), start, fixed, np.asarray(scales), settings,
                               checkpoint_path=part_path.removesuffix('.json')+'.progress.json')
    multi.write_ledger(part_path, record)
    return part_path


def profile_ledger(model, context, scales, start, gate, settings, specifications, cartesian_product, point_count):
    """Profile ledger: the prepared-ledger fields without the multistart rules, plus the profile definition."""
    ledger = prepare_ledger(model, context, scales, [start], gate, settings)
    for key in ('entries', 'acceptance_rule', 'basin_switch_rule', 'control_rule'):
        ledger.pop(key)
    names = model.parameter_names
    ledger.update(
        study=PROFILE_STUDY, restart_source=start['restart_source'],
        profile_definition=dict(
            grids_gev_minus2=[{names[BASELINE_COUNT+index-1]: list(grid)} for index, grid in specifications],
            cartesian_product=bool(cartesian_product), point_count=point_count,
            rule=('each point: the listed r held fixed at the grid value; the 32 baseline coordinates and the other '
                  'r refit at the fit order from the restart endpoint with the fixed value substituted; '
                  'chi2 per block re-evaluated at the evaluation order'),
            not_run='basin switch, zero_control run and quadrature ladder (choices)'),
        profile=[])
    return ledger


def prepare_ledger(model, context, scales, starts, gate, settings=OPTIMIZER_SETTINGS):
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
        optimizer_settings=dict(settings), optimizer_coordinate_scales=list(map(float, scales)),
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


def restart_dry_run(model, context, restart_from, summary):
    """Restart part of the dry run: one residual and Jacobian at the restart point."""
    definition, start_endpoint = restart_definition(model, restart_from)
    start = restart_start_record(model, definition, start_endpoint)
    _, jacobian, evaluation = evaluate_start(model, start['start'], 0)
    control = start_from(model, context['accepted'][0]['endpoint'], STARTS[1])[0]
    scales = optimizer_scales(context, model.jacobian(control))
    source = json.loads((ROOT/restart_from).read_bytes())
    summary.update(
        restart_source=definition['restart_source'], restart_start_endpoint=start['start_endpoint'],
        restart_coefficients_gev_minus2=list(definition['coefficients_gev_minus2']),
        restart_evaluation=evaluation,
        added_optimizer_scales=dict(zip(vpp_option.ADDED_COORDINATE_NAMES, scales[BASELINE_COUNT:].tolist())),
        optimizer_scales_equal_source_ledger=bool(np.array_equal(
            scales, np.asarray(source['optimizer_coordinate_scales'], dtype=float))),
        baseline_chi2_high_order=context['baseline_chi2'])
    return summary


def profile_time_estimate(point_count, workers):
    """Wall-time estimate of a profile under the assumption PROFILE_MINUTES_PER_POINT per point."""
    return dict(assumption=(f'{PROFILE_MINUTES_PER_POINT:g} minutes per point (measured value supplied with the '
                            'task, not measured here); parallel estimate assumes no slowdown from sharing the CPU'),
                minutes_per_point=PROFILE_MINUTES_PER_POINT, point_count=point_count, workers=workers,
                serial_hours=point_count*PROFILE_MINUTES_PER_POINT/60,
                wall_hours_with_workers=math.ceil(point_count/workers)*PROFILE_MINUTES_PER_POINT/60)


def profile_dry_run(model, context, restart_from, specifications, cartesian_product, workers, output, summary):
    """Profile part of the dry run: one residual and Jacobian at the first point of each scan."""
    definition, start_endpoint = restart_definition(model, restart_from)
    start = restart_start_record(model, definition, start_endpoint)
    points = profile_points(specifications, cartesian_product)
    lower, upper = model.parameter_bounds()
    for fixed in points:
        validate_point(model, profile_start(start['start'], fixed, cartesian_product), lower, upper)
    names = model.parameter_names
    evaluations = []
    for fixed in first_profile_points(specifications, cartesian_product):
        _, _, evaluation = evaluate_start(model, profile_start(start['start'], fixed, cartesian_product), len(fixed))
        evaluations.append(dict(fixed_values_gev_minus2={names[index]: value for index, value in fixed.items()},
                                **evaluation))
    control = start_from(model, context['accepted'][0]['endpoint'], STARTS[1])[0]
    scales = optimizer_scales(context, model.jacobian(control))
    summary.update(
        restart_source=definition['restart_source'],
        restart_coefficients_gev_minus2=list(definition['coefficients_gev_minus2']),
        profile=dict(grids_gev_minus2=[{names[BASELINE_COUNT+index-1]: list(grid)} for index, grid in specifications],
                     cartesian_product=bool(cartesian_product), point_count=len(points),
                     first_point_evaluations=evaluations,
                     time_estimate=profile_time_estimate(len(points), workers)),
        added_optimizer_scales=dict(zip(vpp_option.ADDED_COORDINATE_NAMES, scales[BASELINE_COUNT:].tolist())),
        output=output, output_exists=(ROOT/output).exists(),
        baseline_chi2_high_order=context['baseline_chi2'])
    return summary


def dry_run(orders=FIT_ORDERS, restart_from=None, settings=OPTIMIZER_SETTINGS, profile_specifications=None,
            cartesian_product=False, workers=1, output=None):
    """Everything up to the optimizer; returns the summary (nothing is written).

    With restart_from, only the restart start is built and evaluated (residual
    and Jacobian once at the restart point, chi2 per block). With
    profile_specifications (needs restart_from) the evaluation is made instead
    at the first point of each profile scan.
    """
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    model = FitLikelihood(*orders, bin_block_size=BLOCK_SIZE, **MODEL_OPTIONS)
    baseline = joint_fit.GlobalLikelihood(*orders, bin_block_size=BLOCK_SIZE)
    summary = dict(dry_run=True, optimizer_called=False, orders=list(orders), bin_block_size=BLOCK_SIZE,
                   optimizer_settings=dict(settings),
                   model=check_model(model), source_hash_gate=source_hash_gate(baseline))
    if profile_specifications:
        if not restart_from:
            raise ValueError('a profile needs restart_from')
        profile_dry_run(model, context, restart_from, profile_specifications, cartesian_product, workers,
                        output, summary)
        summary['seconds'] = time.perf_counter()-started
        return summary
    if restart_from:
        restart_dry_run(model, context, restart_from, summary)
        summary['seconds'] = time.perf_counter()-started
        return summary
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
    name, part_path, scales, settings = arguments
    definition = next(item for item in STARTS if item['name'] == name)
    result = run_start(definition, sensitivity.accepted_context(), np.asarray(scales), settings=settings,
                       checkpoint_prefix=part_path.removesuffix('.json'))
    multi.write_ledger(part_path, result)
    return name


def main_restart(args, settings):
    """Restart fit: ledger with the single restart start, then one refit (in this process)."""
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    low = make_models((FIT_ORDERS,))[0]
    check_model(low)
    gate = source_hash_gate(joint_fit.GlobalLikelihood(*FIT_ORDERS, bin_block_size=BLOCK_SIZE))
    definition, start_endpoint = restart_definition(low, args.restart_from)
    start = restart_start_record(low, definition, start_endpoint)
    control = start_from(low, context['accepted'][0]['endpoint'], STARTS[1])[0]
    scales = optimizer_scales(context, low.jacobian(control))
    ledger = prepare_ledger(low, context, scales, [start], gate, settings)
    ledger.update(command=shlex.join([sys.executable, *sys.argv]), restart_source=definition['restart_source'])
    multi.write_ledger(ROOT/args.output, ledger)
    parts = ROOT/(args.output+'.parts')
    parts.mkdir(parents=True, exist_ok=True)
    part = parts/f'{definition["name"]}.json'
    if not part.exists():
        result = run_restart(definition, start_endpoint, context, scales, settings=settings,
                             checkpoint_prefix=str(part).removesuffix('.json'))
        multi.write_ledger(part, result)
    ledger['entries'] = [json.loads(part.read_text())]
    ledger.update(state='finished', driver_elapsed_seconds=time.perf_counter()-started)
    multi.write_ledger(ROOT/args.output, ledger)
    print(args.output)


def main_profile(args, settings):
    """Profile in r: ledger with the restart start and the grid, then one refit per point."""
    output = ROOT/args.output
    refuse_existing_output(output, args.overwrite)
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    low = make_models((FIT_ORDERS,))[0]
    check_model(low)
    gate = source_hash_gate(joint_fit.GlobalLikelihood(*FIT_ORDERS, bin_block_size=BLOCK_SIZE))
    definition, start_endpoint = restart_definition(low, args.restart_from)
    start = restart_start_record(low, definition, start_endpoint)
    control = start_from(low, context['accepted'][0]['endpoint'], STARTS[1])[0]
    scales = optimizer_scales(context, low.jacobian(control))
    points = profile_points(args.profile_r, args.profile_grid)
    lower, upper = low.parameter_bounds()
    starts = [validate_point(low, profile_start(start['start'], fixed, args.profile_grid), lower, upper)
              for fixed in points]
    print('profile points:', len(points), flush=True)
    ledger = profile_ledger(low, context, scales, start, gate, settings, args.profile_r, args.profile_grid,
                            len(points))
    ledger['command'] = shlex.join([sys.executable, *sys.argv])
    multi.write_ledger(output, ledger)
    parts = ROOT/(args.output+'.parts')
    parts.mkdir(parents=True, exist_ok=True)
    paths = [str(parts/profile_part_name(low, fixed)) for fixed in points]
    pending = [(path, values.tolist(), fixed, scales.tolist(), settings)
               for path, values, fixed in zip(paths, starts, points)
               if not profile_part_matches(path, values, settings)]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for path in executor.map(run_profile_named, pending):
            print('finished', Path(path).name, flush=True)
    ledger['profile'] = [json.loads(Path(path).read_text()) for path in paths]
    ledger.update(state='finished', driver_elapsed_seconds=time.perf_counter()-started)
    multi.write_ledger(output, ledger)
    print(args.output)


def main(argv=None):
    args = parse_arguments(argv)
    settings = optimizer_settings(args.max_nfev, args.ftol, args.xtol, args.gtol)
    if args.dry_run:
        print(json.dumps(dry_run(restart_from=args.restart_from, settings=settings,
                                 profile_specifications=args.profile_r, cartesian_product=args.profile_grid,
                                 workers=args.workers, output=args.output), indent=1))
        return
    if args.profile_r:
        main_profile(args, settings)
        return
    if args.restart_from:
        main_restart(args, settings)
        return
    started = time.perf_counter()
    context = sensitivity.accepted_context()
    low = make_models((FIT_ORDERS,))[0]
    check_model(low)
    gate = source_hash_gate(joint_fit.GlobalLikelihood(*FIT_ORDERS, bin_block_size=BLOCK_SIZE))
    starts = start_points(low, context)
    control = next(item for item in starts if item['name'] == 'zero_control')
    scales = optimizer_scales(context, low.jacobian(np.asarray(control['start'])))
    ledger = prepare_ledger(low, context, scales, starts, gate, settings)
    ledger['command'] = shlex.join([sys.executable, *sys.argv])
    multi.write_ledger(ROOT/args.output, ledger)
    parts = ROOT/(args.output+'.parts')
    parts.mkdir(parents=True, exist_ok=True)
    pending = [(item['name'], str(parts/f'{item["name"]}.json'), scales.tolist(), settings) for item in STARTS
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
