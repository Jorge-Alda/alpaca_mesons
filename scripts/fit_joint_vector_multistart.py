"""Run deterministic starts of the existing joint vector-sector likelihood.

Example: .venv/bin/python scripts/fit_joint_vector_multistart.py --prepare-only --output /tmp/joint-starts.json
Masses and widths are in GeV, cross sections in nb, and residuals are dimensionless.
The recorded endpoints are not accepted minima until a separate refinement study.
"""
import argparse
import hashlib
import json
import platform
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


def prepare_ledger(model, orders, block_size, settings):
    reference = json.loads((ROOT / REFERENCE).read_text())['parameter_values']
    starts = named_starts(model, reference)
    lower, upper = joint_fit.parameter_bounds(model)
    scales = optimizer_scales(model, starts['stored'])
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
            'software': {'python': sys.version, 'numpy': np.__version__,
                         'scipy': scipy.__version__, 'platform': platform.platform()},
            'attempts': [{'name': name, 'start': values.tolist(), 'state': 'pending',
                          'accepted': False} for name, values in starts.items()]}


def write_ledger(path, ledger):
    """Persist every attempt without dropping unsuccessful or distinct endpoints."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(ledger, indent=2, allow_nan=False) + '\n')


def run_attempt(model, attempt, lower, upper, scales, settings):
    """Minimize once and retain optimizer diagnostics, including failures."""
    start = validate_start(model, attempt['start'], lower, upper)
    start_time = time.perf_counter()
    try:
        result = least_squares(model.residual, start, jac=model.jacobian,
                               bounds=(lower, upper), x_scale=scales, **settings)
        attempt.update(state='finished', endpoint=result.x.tolist(),
                       residuals=result.fun.tolist(), chi2=float(result.fun @ result.fun),
                       success=bool(result.success), status=int(result.status), message=str(result.message),
                       nfev=int(result.nfev), njev=None if result.njev is None else int(result.njev),
                       optimality=float(result.optimality), gradient=result.grad.tolist(),
                       active_bounds=result.active_mask.tolist())
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        attempt.update(state='error', success=False, error_type=type(error).__name__, message=str(error))
    attempt['elapsed_seconds'] = time.perf_counter() - start_time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--max-nfev', type=int, default=1500)
    parser.add_argument('--ftol', type=float, default=1e-8)
    parser.add_argument('--xtol', type=float, default=1e-9)
    parser.add_argument('--gtol', type=float, default=1e-7)
    parser.add_argument('--energy-order', type=int, default=2)
    parser.add_argument('--pair-order', type=int, default=16)
    parser.add_argument('--angle-order', type=int, default=16)
    parser.add_argument('--bin-block-size', type=int, default=4)
    args = parser.parse_args()
    if min(args.energy_order, args.pair_order, args.angle_order, args.bin_block_size, args.max_nfev) <= 0:
        parser.error('orders, block size and max-nfev must be positive')
    orders = (args.energy_order, args.pair_order, args.angle_order)
    model = joint_fit.GlobalLikelihood(*orders, bin_block_size=args.bin_block_size)
    settings = {'method': 'trf', 'max_nfev': args.max_nfev, 'ftol': args.ftol,
                'xtol': args.xtol, 'gtol': args.gtol}
    ledger = prepare_ledger(model, orders, args.bin_block_size, settings)
    write_ledger(args.output, ledger)
    if not args.prepare_only:
        lower, upper = joint_fit.parameter_bounds(model)
        for attempt in ledger['attempts']:
            attempt['state'] = 'running'
            write_ledger(args.output, ledger)
            run_attempt(model, attempt, lower, upper, np.array(ledger['optimizer_coordinate_scales']), settings)
            write_ledger(args.output, ledger)
    print(args.output)


if __name__ == '__main__':
    main()
