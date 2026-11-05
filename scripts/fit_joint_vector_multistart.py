"""Run deterministic starts of the existing joint vector-sector likelihood.

Example: .venv/bin/python scripts/fit_joint_vector_multistart.py --prepare-only --output /tmp/joint-starts.json
Masses and widths are in GeV, cross sections in nb, and residuals are dimensionless.
The recorded endpoints are not accepted minima until a separate refinement study.
"""
import argparse
import hashlib
import json
import platform
import shlex
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares

import fit_babar_kkpi_global_vector_sector as joint_fit


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = 'data/vector_sector/babar_kkpi_global_vector_sector_fit.json'
INPUTS = (REFERENCE, 'data/vector_sector/odd_preliminary_fit.json',
          'data/vector_sector/stage3_results.json', 'data/vector_sector/pdg2024.json',
          'data/vector_sector/alda_2026nrd_arxiv_v1.json', 'data/vector_sector/gv_measurements.json',
          'data/odd_parity/odd_low_energy_measurements.json',
          *(f'data/vector_sector/babar_2007_{channel}_inclusive_below_1p78.csv'
            for channel in joint_fit.CHANNELS))
BLOCKS = (('odd', 0, 12), ('masses_and_leptonic', 12, 19),
          ('shared_widths', 19, 26), ('phi_KK', 26, 28),
          ('BaBar', 28, 73), ('normalization_constraints', 73, 75))
CERTIFICATION_ORDERS = {'energy': [2, 4], 'pair': [64, 128], 'angle': [64, 128]}


def validate_start(model, values, lower, upper):
    """Check the bounds and the physical branches used by the likelihood."""
    values = np.asarray(values, dtype=float)
    if values.shape != (32,) or not np.all(np.isfinite(values)):
        raise ValueError('start must have 32 finite coordinates')
    if np.any(values < lower) or np.any(values > upper):
        raise ValueError('start is outside parameter bounds')
    model.unpack(values)
    return values


def named_starts(model, reference):
    """Use the stored point and reproducible sign variants on the same branch."""
    lower, upper = joint_fit.parameter_bounds(model)
    preliminary = model.initial()
    reference = np.asarray(reference, dtype=float)
    negative_c4 = reference.copy()
    negative_c4[16] *= -1
    negative_odd = reference.copy()
    negative_odd[12:24] *= -1
    starts = {'stored': reference, 'preliminary': preliminary,
              'stored_negative_c4': negative_c4, 'stored_negative_odd': negative_odd}
    return {name: validate_start(model, values, lower, upper).copy()
            for name, values in starts.items()}


def optimizer_scales(model, reference):
    """Compute one fixed coordinate scale from the stored-point Jacobian."""
    column_norms = np.linalg.norm(model.jacobian(reference), axis=0)
    if not np.all(np.isfinite(column_norms)) or np.any(column_norms <= 0):
        raise ValueError('stored-point Jacobian has a zero or nonfinite column norm')
    return 1 / column_norms


def json_bounds(values):
    """Represent unbounded endpoints by JSON null."""
    return [float(value) if np.isfinite(value) else None for value in values]


def input_hashes():
    """Identify the direct data inputs, starting point, and fit code."""
    paths = (*INPUTS, 'scripts/fit_joint_vector_multistart.py',
             'scripts/fit_babar_kkpi_global_vector_sector.py',
             'scripts/fit_babar_kkpi_g_coefficients.py', 'scripts/calibrate_odd_low_energy.py',
             'scripts/fit_vector_gv.py')
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in paths}


def source_hashes():
    """Hash project Python modules loaded by the likelihood implementation."""
    paths = {Path(module.__file__).resolve() for module in tuple(sys.modules.values())
             if getattr(module, '__file__', None) and Path(module.__file__).suffix == '.py'
             and Path(module.__file__).resolve().is_relative_to(ROOT)
             and not Path(module.__file__).resolve().is_relative_to(ROOT / '.venv')}
    paths.add(Path(__file__).resolve())
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def prepare_ledger(model, orders, block_size, settings, optimizer_coordinate_scales=None):
    reference = json.loads((ROOT / REFERENCE).read_text())['parameter_values']
    starts = named_starts(model, reference)
    lower, upper = joint_fit.parameter_bounds(model)
    scales = (optimizer_scales(model, starts['stored']) if optimizer_coordinate_scales is None
              else np.asarray(optimizer_coordinate_scales, dtype=float))
    return {'parameter_names': list(joint_fit.PARAMETER_NAMES),
            'bounds': {'lower': json_bounds(lower), 'upper': json_bounds(upper)},
            'physical_branch': {'source': 'single_trace_negative',
                'inverse_mass_branch_sign': float(model.branch_sign),
                'conditions': ['light source positive', 'strange source negative',
                               'G_V positive', 'total widths positive', 'absolute excited residues at most 5']},
            'even_physical_scales': joint_fit.EVEN_SCALES.tolist(),
            'optimizer_coordinate_scales': scales.tolist(),
            'scale_definition': 'reciprocal Jacobian column norms at stored start, fixed for all starts',
            'quadrature': {'energy': orders[0], 'pair_per_segment': orders[1],
                           'angle_per_segment': orders[2], 'bin_block_size': block_size},
            'optimizer_settings': settings, 'input_sha256': input_hashes(),
            'source_sha256': source_hashes(),
            'software': {'python': sys.version, 'numpy': np.__version__,
                         'scipy': scipy.__version__, 'platform': platform.platform()},
            'attempts': [{'name': name, 'start': values.tolist(), 'state': 'pending',
                          'accepted': False} for name, values in starts.items()]}


def write_ledger(path, ledger):
    """Persist every attempt without dropping unsuccessful or distinct endpoints."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(ledger, indent=2, allow_nan=False) + '\n')


def endpoint_diagnostics(model, values, lower, upper):
    """Evaluate domain validity, residual blocks and BaBar predictions."""
    try:
        values = validate_start(model, values, lower, upper)
        residual = model.residual(values)
        if residual.shape != (75,) or not np.all(np.isfinite(residual)):
            raise ValueError('endpoint residuals are not 75 finite values')
        if np.all(residual == 1e8 + np.linalg.norm(values)):
            raise ValueError('likelihood returned its invalid-point penalty')
        _, predictions = model.babar_blocks(values, model.unpack(values))
        if not all(np.all(np.isfinite(prediction)) for prediction in predictions.values()):
            raise ValueError('endpoint predictions are nonfinite')
        return {'valid': True, 'residuals': residual.tolist(),
                'chi2': float(residual @ residual),
                'blocks': {name: {'count': end - start,
                                  'chi2': float(residual[start:end] @ residual[start:end])}
                           for name, start, end in BLOCKS},
                'babar_predictions': {name: values.tolist() for name, values in predictions.items()}}
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        return {'valid': False, 'error_type': type(error).__name__, 'message': str(error)}


def run_attempt(model, attempt, lower, upper, scales, settings, callback=None):
    """Minimize once and retain optimizer diagnostics, including failures."""
    start = validate_start(model, attempt['start'], lower, upper)
    start_time = time.perf_counter()
    try:
        options = dict(settings)
        if callback is not None:
            options['callback'] = callback
        result = least_squares(model.residual, start, jac=model.jacobian,
                               bounds=(lower, upper), x_scale=scales, **options)
        attempt.update(state='finished', endpoint=result.x.tolist(),
                       residuals=result.fun.tolist(), chi2=float(result.fun @ result.fun),
                       success=bool(result.success), status=int(result.status), message=str(result.message),
                       nfev=int(result.nfev), njev=None if result.njev is None else int(result.njev),
                       optimality=float(result.optimality), gradient=result.grad.tolist(),
                       active_bounds=result.active_mask.tolist())
        attempt['endpoint_diagnostics'] = endpoint_diagnostics(model, result.x, lower, upper)
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        attempt.update(state='error', success=False, error_type=type(error).__name__, message=str(error))
    attempt['elapsed_seconds'] = time.perf_counter() - start_time


def prepare_continuation(model, source_path, start_name, orders, block_size, settings):
    """Load an endpoint with its original fixed optimizer scale and provenance."""
    source_path = Path(source_path)
    source_bytes = source_path.read_bytes()
    source = json.loads(source_bytes)
    if source['parameter_names'] != list(joint_fit.PARAMETER_NAMES):
        raise ValueError('source parameter coordinates differ from the current fit')
    # The continuation driver is expected to change when resume support is
    # introduced. The likelihood implementation and data must remain identical.
    changed_driver = str(Path(__file__).resolve().relative_to(ROOT))
    for recorded_hashes in (source['input_sha256'], source['source_sha256']):
        for path, expected in recorded_hashes.items():
            if path != changed_driver and hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
                raise ValueError(f'source data or likelihood code changed: {path}')
    matches = [item for item in source['attempts'] if item['name'] == start_name]
    if len(matches) != 1:
        raise ValueError('source ledger must contain exactly one named attempt')
    prior = matches[0]
    if prior.get('endpoint_diagnostics', {}).get('valid') is not True:
        raise ValueError('source endpoint was not marked physically valid')
    lower, upper = joint_fit.parameter_bounds(model)
    if source['bounds'] != {'lower': json_bounds(lower), 'upper': json_bounds(upper)}:
        raise ValueError('current parameter bounds differ from source ledger')
    if source['physical_branch']['inverse_mass_branch_sign'] != float(model.branch_sign):
        raise ValueError('current inverse-mass branch differs from source ledger')
    values = validate_start(model, prior['endpoint'], lower, upper)
    scales = np.asarray(source['optimizer_coordinate_scales'], dtype=float)
    if scales.shape != (32,) or not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError('source optimizer scales must contain 32 positive finite values')
    prior_nfev = prior.get('cumulative_nfev', prior.get('nfev'))
    if type(prior_nfev) is not int or prior_nfev < 0:
        raise ValueError('source attempt has no valid function evaluation count')
    ledger = prepare_ledger(model, orders, block_size, settings, scales)
    ledger['optimizer_coordinate_scales'] = source['optimizer_coordinate_scales']
    ledger['scale_definition'] = source['scale_definition']
    ledger['resume_source'] = {'path': str(source_path),
                               'sha256': hashlib.sha256(source_bytes).hexdigest(),
                               'attempt': start_name, 'source_nfev': prior_nfev,
                               'source_success': prior['success'],
                               'source_driver_sha256': source['input_sha256'].get(changed_driver),
                               'current_driver_sha256': ledger['input_sha256'][changed_driver]}
    ledger['attempts'] = [{'name': start_name, 'start': values.tolist(),
                           'state': 'pending', 'accepted': False,
                           'prior_nfev': prior_nfev}]
    if ledger['attempts'][0]['start'] != prior['endpoint']:
        raise ValueError('continuation start differs from source endpoint')
    return ledger


def write_checkpoint(path, intermediate_result, source_nfev):
    """Save the current SciPy iteration without modifying the source ledger."""
    residual = np.asarray(intermediate_result.fun, dtype=float)
    values = np.asarray(intermediate_result.x, dtype=float)
    nfev = int(intermediate_result.nfev)
    available_optimality = getattr(intermediate_result, 'optimality', None)
    checkpoint = {'parameter_values': values.tolist(), 'chi2': float(residual @ residual),
                  'nfev': nfev, 'prior_nfev': source_nfev,
                  'cumulative_nfev': source_nfev + nfev,
                  'optimality': None if available_optimality is None else float(available_optimality),
                  'iteration': int(intermediate_result.nit)}
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(checkpoint, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def refine_attempt(model, attempt, lower, upper):
    """Apply the historical high-order quadrature audit to a successful endpoint."""
    if not attempt.get('success') or not attempt.get('endpoint_diagnostics', {}).get('valid'):
        return {'status': 'not_evaluated', 'passed': False,
                'reason': 'optimizer did not succeed at a valid endpoint'}
    try:
        values = validate_start(model, attempt['endpoint'], lower, upper)
        started = time.perf_counter()
        refinement = joint_fit.audit_quadrature(model, values, CERTIFICATION_ORDERS,
                                                parameter_scales=np.ones(32))
        refinement['elapsed_seconds'] = time.perf_counter() - started
        refinement['bin_block_size'] = model.bin_block_size
        return refinement
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        return {'status': 'unevaluated', 'passed': False,
                'error_type': type(error).__name__, 'reason': str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--refine-from', help='Refine successful endpoints in this saved ledger')
    parser.add_argument('--resume-from', help='Continue a named endpoint from this saved ledger')
    parser.add_argument('--resume-start', help='Name of the saved endpoint to continue')
    parser.add_argument('--starts', nargs='+', help='Names of starts to run or refine')
    parser.add_argument('--max-nfev', type=int, default=1500)
    parser.add_argument('--ftol', type=float, default=1e-8)
    parser.add_argument('--xtol', type=float, default=1e-9)
    parser.add_argument('--gtol', type=float, default=1e-7)
    parser.add_argument('--energy-order', type=int, default=2)
    parser.add_argument('--pair-order', type=int, default=16)
    parser.add_argument('--angle-order', type=int, default=16)
    parser.add_argument('--bin-block-size', type=int, default=4)
    args = parser.parse_args()
    if bool(args.resume_from) != bool(args.resume_start):
        parser.error('--resume-from and --resume-start must be supplied together')
    if args.resume_from and (args.refine_from or args.starts):
        parser.error('--resume-from cannot be combined with --refine-from or --starts')
    if args.refine_from:
        if args.prepare_only:
            parser.error('--prepare-only cannot be combined with --refine-from')
        ledger = json.loads(Path(args.refine_from).read_text())
        if ledger['parameter_names'] != list(joint_fit.PARAMETER_NAMES):
            parser.error('saved parameter coordinates differ from the current fit')
        for path, expected in ledger['input_sha256'].items():
            if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
                parser.error(f'input or fit code changed: {path}')
        for path, expected in ledger['source_sha256'].items():
            if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
                parser.error(f'loaded source changed: {path}')
        if Path(args.output).resolve() == Path(args.refine_from).resolve():
            parser.error('refinement output must differ from its input ledger')
        model = joint_fit.GlobalLikelihood(2, 64, 64, bin_block_size=4)
        lower, upper = joint_fit.parameter_bounds(model)
        selected = set(args.starts or [attempt['name'] for attempt in ledger['attempts']])
        if selected - {attempt['name'] for attempt in ledger['attempts']}:
            parser.error('unknown start name')
        ledger['refinement_command'] = shlex.join([sys.executable, *sys.argv])
        for attempt in ledger['attempts']:
            if attempt['name'] in selected:
                attempt['quadrature_refinement'] = refine_attempt(model, attempt, lower, upper)
                write_ledger(args.output, ledger)
        print(args.output)
        return
    if min(args.energy_order, args.pair_order, args.angle_order, args.bin_block_size, args.max_nfev) <= 0:
        parser.error('orders, block size and max-nfev must be positive')
    orders = (args.energy_order, args.pair_order, args.angle_order)
    model = joint_fit.GlobalLikelihood(*orders, bin_block_size=args.bin_block_size)
    settings = {'method': 'trf', 'max_nfev': args.max_nfev, 'ftol': args.ftol,
                'xtol': args.xtol, 'gtol': args.gtol}
    if args.resume_from:
        checkpoint_path = Path(args.output).with_suffix('.progress.json')
        source_path = Path(args.resume_from).resolve()
        if Path(args.output).resolve() == source_path or checkpoint_path.resolve() == source_path:
            parser.error('output and checkpoint must not overwrite the source ledger')
        ledger = prepare_continuation(model, args.resume_from, args.resume_start,
                                      orders, args.bin_block_size, settings)
    else:
        ledger = prepare_ledger(model, orders, args.bin_block_size, settings)
    ledger['command'] = shlex.join([sys.executable, *sys.argv])
    ledger['working_directory'] = str(Path.cwd())
    selected = set(args.starts or [attempt['name'] for attempt in ledger['attempts']])
    if selected - {attempt['name'] for attempt in ledger['attempts']}:
        parser.error('unknown start name')
    write_ledger(args.output, ledger)
    if not args.prepare_only:
        lower, upper = joint_fit.parameter_bounds(model)
        for attempt in ledger['attempts']:
            if attempt['name'] not in selected:
                continue
            attempt['state'] = 'running'
            write_ledger(args.output, ledger)
            callback = None
            if args.resume_from:
                def callback(intermediate_result):
                    write_checkpoint(checkpoint_path, intermediate_result, attempt['prior_nfev'])
            run_attempt(model, attempt, lower, upper, np.array(ledger['optimizer_coordinate_scales']), settings, callback)
            if args.resume_from and 'nfev' in attempt:
                attempt['cumulative_nfev'] = attempt['prior_nfev'] + attempt['nfev']
            write_ledger(args.output, ledger)
    print(args.output)


if __name__ == '__main__':
    main()
