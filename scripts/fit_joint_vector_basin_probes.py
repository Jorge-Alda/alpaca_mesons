"""Fit three deterministic starts near the high-order joint-fit endpoint.

Run: .venv/bin/python scripts/fit_joint_vector_basin_probes.py --source data/vector_sector/joint_vector_continuation_e4p128a128.json --output /tmp/joint-vector-probes.json --refine-successful
Coordinates use the original fit units; residuals and scaled distances are
dimensionless. The nearby starts probe the search and are not confidence regions.
"""
import argparse
import json
import shlex
import sys
from pathlib import Path

import numpy as np

import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_joint_vector_multistart as multi


BLOCKS = {'even_shared_widths': tuple(range(0, 12)),
          'odd': tuple(range(12, 24)),
          'babar': tuple(range(24, 32))}
SCALED_RADIUS = 0.25
MAXIMUM_HALVINGS = 40


def basin_starts(model, seed, scales):
    """Displace each block by a signed vector of scaled length 0.25.

    Alternate signs within each block and halve the common displacement until
    the existing bounds and physical branch are satisfied. Displacements are
    numerical search choices and are not uncertainty estimates.
    """
    lower, upper = joint_fit.parameter_bounds(model)
    seed = multi.validate_start(model, seed, lower, upper)
    scales = np.asarray(scales, dtype=float)
    if scales.shape != (32,) or not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError('expected 32 positive finite optimizer scales')
    starts = []
    for name, indices in BLOCKS.items():
        direction = np.zeros(32)
        for offset, index in enumerate(indices):
            direction[index] = (1 if offset % 2 == 0 else -1) / np.sqrt(len(indices))
        for halvings in range(MAXIMUM_HALVINGS + 1):
            values = seed + SCALED_RADIUS * (0.5 ** halvings) * scales * direction
            if np.array_equal(values, seed):
                raise ValueError(f'{name}: perturbation is not representable')
            try:
                multi.validate_start(model, values, lower, upper)
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                continue
            starts.append({'name': name, 'start': values.tolist(), 'state': 'pending',
                           'accepted': False, 'perturbed_indices': list(indices),
                           'initial_scaled_radius': SCALED_RADIUS,
                           'halvings': halvings,
                           'realized_scaled_distance': float(np.linalg.norm((values-seed)/scales))})
            break
        else:
            raise ValueError(f'{name}: no valid perturbation after {MAXIMUM_HALVINGS} halvings')
    return starts


def prepare_probes(model, source_path, source_name, max_nfev):
    """Create a ledger using the high-order endpoint and its recorded scales."""
    source = json.loads(Path(source_path).read_bytes())
    quadrature = source['quadrature']
    if (quadrature['energy'], quadrature['pair_per_segment'], quadrature['angle_per_segment'],
            quadrature['bin_block_size']) != (4, 128, 128, 4):
        raise ValueError('source must be the four-bin (4,128,128) continuation')
    matches = [attempt for attempt in source['attempts'] if attempt['name'] == source_name]
    if len(matches) != 1 or matches[0].get('success') is not True:
        raise ValueError('source must contain one successful named endpoint')
    settings = {'method': 'trf', 'max_nfev': max_nfev, 'ftol': 1e-8,
                'xtol': 1e-9, 'gtol': 1e-7}
    ledger = multi.prepare_continuation(model, source_path, source_name, (2, 16, 16), 4, settings)
    seed = ledger['attempts'][0]['start']
    ledger['probe_seed'] = {'name': source_name, 'endpoint': seed,
                            'chi2': float(matches[0]['chi2']), 'success': True,
                            'endpoint_valid': True}
    ledger['attempts'] = basin_starts(model, seed, ledger['optimizer_coordinate_scales'])
    ledger['perturbation'] = {'scaled_radius': SCALED_RADIUS,
                              'direction': 'alternating signs, equal scaled magnitude per coordinate',
                              'domain_repair': 'halve the displacement until valid',
                              'maximum_halvings': MAXIMUM_HALVINGS,
                              'interpretation': 'deterministic search, not a confidence region'}
    return ledger


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--source-name', default='stored_negative_c4')
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-nfev', type=int, default=1500)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--refine-successful', action='store_true',
                        help='Run the certified quadrature audit on every successful valid probe')
    args = parser.parse_args()
    if args.max_nfev <= 0:
        parser.error('max-nfev must be positive')
    if Path(args.source).resolve() == Path(args.output).resolve():
        parser.error('output must differ from source')
    model = joint_fit.GlobalLikelihood(2, 16, 16, bin_block_size=4)
    ledger = prepare_probes(model, args.source, args.source_name, args.max_nfev)
    ledger['command'] = shlex.join([sys.executable, *sys.argv])
    ledger['working_directory'] = str(Path.cwd())
    multi.write_ledger(args.output, ledger)
    if not args.prepare_only:
        lower, upper = joint_fit.parameter_bounds(model)
        scales = np.asarray(ledger['optimizer_coordinate_scales'])
        for attempt in ledger['attempts']:
            attempt['state'] = 'running'
            multi.write_ledger(args.output, ledger)
            checkpoint_path = Path(args.output).with_name(
                Path(args.output).stem + '.' + attempt['name'] + '.progress.json')
            def callback(intermediate_result):
                multi.write_checkpoint(checkpoint_path, intermediate_result, 0)
            multi.run_attempt(model, attempt, lower, upper, scales, ledger['optimizer_settings'], callback)
            multi.write_ledger(args.output, ledger)
        if args.refine_successful:
            refinement_model = joint_fit.GlobalLikelihood(2, 64, 64, bin_block_size=4)
            for attempt in ledger['attempts']:
                attempt['quadrature_refinement'] = multi.refine_attempt(refinement_model, attempt, lower, upper)
                multi.write_ledger(args.output, ledger)
    print(args.output)


if __name__ == '__main__':
    main()
