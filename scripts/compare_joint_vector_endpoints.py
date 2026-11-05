"""Compare saved joint-fit endpoints at common high quadrature order.

Run: .venv/bin/python scripts/compare_joint_vector_endpoints.py --output /tmp/joint-vector-comparison.json
Masses and widths are in GeV, cross sections in nb, and scaled distances and
chi-squared values are dimensionless. Pairwise agreement does not prove that
the search found a global minimum.
"""
import argparse
import hashlib
import json
import shlex
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_joint_vector_multistart as multi


ROOT = joint_fit.ROOT
DEFAULT_INPUTS = {
    'multistart': 'data/vector_sector/joint_vector_multistart.json',
    'continuation': 'data/vector_sector/joint_vector_continuation.json',
    'continuation_26464': 'data/vector_sector/joint_vector_continuation_e26464.json',
    'continuation_4128128': 'data/vector_sector/joint_vector_continuation_e4p128a128.json',
    'continuation_audit': 'data/vector_sector/joint_vector_continuation_e4p128a128_refinement.json',
    'basin_probes': 'data/vector_sector/joint_vector_basin_probes.json',
}
EXPECTED_NAMES = {
    'multistart': {'stored', 'preliminary', 'stored_negative_c4', 'stored_negative_odd'},
    'continuation': {'stored_negative_c4'}, 'continuation_26464': {'stored_negative_c4'},
    'continuation_4128128': {'stored_negative_c4'},
    'continuation_audit': {'stored_negative_c4'},
    'basin_probes': {'even_shared_widths', 'odd', 'babar'},
}
EXPECTED_ORDERS = {'multistart': (2, 16, 16), 'continuation': (2, 16, 16),
                   'continuation_26464': (2, 64, 64), 'continuation_4128128': (4, 128, 128),
                   'continuation_audit': (4, 128, 128), 'basin_probes': (2, 16, 16)}
SCALED_DISTANCE_TOLERANCE = .01
CHI2_TOLERANCE = .01


def read_ledgers(paths):
    """Read exact source bytes and retain their SHA-256 values."""
    ledgers, hashes = {}, {}
    for role, path in paths.items():
        data = Path(path).read_bytes()
        ledgers[role] = json.loads(data)
        hashes[role] = {'path': str(path), 'sha256': hashlib.sha256(data).hexdigest()}
    return ledgers, hashes


def audit_passed(audit):
    """Require the complete previously certified four-comparison ladder."""
    if not isinstance(audit, dict):
        return False
    if (audit.get('status') != 'certified' or audit.get('passed') is not True
            or audit.get('orders') != multi.CERTIFICATION_ORDERS
            or audit.get('bin_block_size') != 4
            or audit.get('tolerances') != joint_fit.QUADRATURE_TOLERANCES):
        return False
    refinements = audit.get('refinements', {})
    return all(refinements.get(name, {}).get('passed') is True
               for name in ('energy', 'pair', 'angle', 'combined'))


def validate_ledgers(ledgers, hashes, model):
    """Check source identity, fit conventions, recorded scales and domains."""
    if set(ledgers) != set(DEFAULT_INPUTS) or set(hashes) != set(DEFAULT_INPUTS):
        raise ValueError('all six named source ledgers are required')
    reference = ledgers['multistart']
    reference_scales = np.asarray(reference['optimizer_coordinate_scales'], dtype=float)
    if (reference_scales.shape != (32,) or not np.all(np.isfinite(reference_scales))
            or np.any(reference_scales <= 0)):
        raise ValueError('invalid optimizer scales')
    lower, upper = joint_fit.parameter_bounds(model)
    expected_bounds = {'lower': multi.json_bounds(lower), 'upper': multi.json_bounds(upper)}
    reference_inputs = {path: value for path, value in reference['input_sha256'].items()
                        if path.startswith('data/')}
    for path, expected in reference_inputs.items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
            raise ValueError(f'likelihood input changed: {path}')
    for role, ledger in ledgers.items():
        quadrature = ledger['quadrature']
        if ((quadrature['energy'], quadrature['pair_per_segment'], quadrature['angle_per_segment'])
                != EXPECTED_ORDERS[role] or quadrature['bin_block_size'] != 4):
            raise ValueError(f'{role}: unexpected fit quadrature')
        if ledger['parameter_names'] != list(joint_fit.PARAMETER_NAMES):
            raise ValueError(f'{role}: parameter names do not match')
        if ledger['optimizer_coordinate_scales'] != reference['optimizer_coordinate_scales']:
            raise ValueError(f'{role}: optimizer scales do not match')
        if ledger['bounds'] != expected_bounds:
            raise ValueError(f'{role}: optimizer bounds do not match')
        if ledger['physical_branch']['inverse_mass_branch_sign'] != float(model.branch_sign):
            raise ValueError(f'{role}: inverse-mass branch does not match')
        inputs = {path: value for path, value in ledger['input_sha256'].items() if path.startswith('data/')}
        if inputs != reference_inputs:
            raise ValueError(f'{role}: recorded input hashes differ')
        # The continuation driver differs between campaigns; all other
        # recorded modules must still match their source hashes.
        for path, expected in ledger['source_sha256'].items():
            if path == 'scripts/fit_joint_vector_multistart.py':
                continue
            if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
                raise ValueError(f'{role}: likelihood source changed: {path}')
        names = [attempt['name'] for attempt in ledger['attempts']]
        if len(names) != len(EXPECTED_NAMES[role]) or set(names) != EXPECTED_NAMES[role]:
            raise ValueError(f'{role}: unexpected attempt names')
        for attempt in ledger['attempts']:
            if attempt['state'] != 'finished' or type(attempt.get('success')) is not bool:
                raise ValueError(f'{role}/{attempt["name"]}: optimizer status is missing')
            if (type(attempt.get('status')) is not int or not isinstance(attempt.get('message'), str)
                    or type(attempt.get('nfev')) is not int or attempt['nfev'] < 1
                    or (attempt.get('njev') is not None and
                        (type(attempt['njev']) is not int or attempt['njev'] < 0))
                    or not np.isfinite(attempt.get('chi2', np.nan))
                    or not np.isfinite(attempt.get('optimality', np.nan))):
                raise ValueError(f'{role}/{attempt["name"]}: incomplete optimizer diagnostics')
            gradient = np.asarray(attempt.get('gradient', []), dtype=float)
            active = np.asarray(attempt.get('active_bounds', []))
            if (bool(attempt['status'] > 0) != attempt['success'] or gradient.shape != (32,)
                    or not np.all(np.isfinite(gradient)) or active.shape != (32,)
                    or not np.all(np.isin(active, (-1, 0, 1)))):
                raise ValueError(f'{role}/{attempt["name"]}: inconsistent optimizer diagnostics')
            if attempt.get('endpoint_diagnostics', {}).get('valid') is not True:
                raise ValueError(f'{role}/{attempt["name"]}: recorded endpoint is invalid')
            multi.validate_start(model, attempt['endpoint'], lower, upper)
    for role, previous in (('continuation', 'multistart'), ('continuation_26464', 'continuation'),
                           ('continuation_4128128', 'continuation_26464'),
                           ('basin_probes', 'continuation_4128128')):
        if ledgers[role]['resume_source']['sha256'] != hashes[previous]['sha256']:
            raise ValueError(f'{role}: source-ledger hash does not match {previous}')
    high = ledgers['continuation_4128128']['attempts'][0]
    audited = ledgers['continuation_audit']['attempts'][0]
    if (high['endpoint'] != audited['endpoint'] or high['success'] != audited['success']
            or high['chi2'] != audited['chi2']):
        raise ValueError('high-order audit is not for the saved endpoint')
    return reference_scales


def compare(ledgers, hashes, model):
    """Evaluate every endpoint at (4,128,128), rank, and compare each pair."""
    scales = validate_ledgers(ledgers, hashes, model)
    if model.bin_block_size != 4 or model.phase.orders != (4, 128, 128):
        raise ValueError('common objective requires (4,128,128) with four-bin blocks')
    candidates = []
    for role in ('multistart', 'continuation', 'continuation_26464',
                 'continuation_4128128', 'basin_probes'):
        for attempt in ledgers[role]['attempts']:
            values = np.asarray(attempt['endpoint'], dtype=float)
            residual = model.residual(values)
            if (residual.shape != (75,) or not np.all(np.isfinite(residual))
                    or np.all(residual == 1e8 + np.linalg.norm(values))):
                raise ValueError(f'{role}/{attempt["name"]}: invalid common-order residual')
            audit = (ledgers['continuation_audit']['attempts'][0].get('quadrature_refinement')
                     if role == 'continuation_4128128' else attempt.get('quadrature_refinement'))
            certified = audit_passed(audit)
            candidates.append({'id': f'{role}/{attempt["name"]}', 'source': role,
                               'name': attempt['name'], 'endpoint': values.tolist(),
                               'source_chi2': float(attempt['chi2']),
                               'common_chi2': float(residual @ residual),
                               'optimizer_success': attempt['success'],
                               'optimizer_status': attempt['status'],
                               'optimizer_message': attempt['message'],
                               'optimizer_optimality': attempt['optimality'],
                               'endpoint_valid': True, 'quadrature_certified': certified,
                               'accepted': bool(attempt['success'] and certified)})
    ranked = sorted(candidates, key=lambda candidate: (candidate['common_chi2'], candidate['id']))
    pairs = []
    for first, second in combinations(candidates, 2):
        scaled_distance = float(np.linalg.norm((np.asarray(first['endpoint']) - np.asarray(second['endpoint'])) / scales))
        chi2_difference = abs(first['common_chi2'] - second['common_chi2'])
        pairs.append({'first': first['id'], 'second': second['id'], 'scaled_distance': scaled_distance,
                      'absolute_chi2_difference': chi2_difference,
                      'within_both_tolerances': bool(scaled_distance <= SCALED_DISTANCE_TOLERANCE
                                                     and chi2_difference <= CHI2_TOLERANCE)})
    return {'sources': hashes, 'quadrature': {'energy': 4, 'pair_per_segment': 128,
                                              'angle_per_segment': 128, 'bin_block_size': 4},
            'parameter_names': list(joint_fit.PARAMETER_NAMES),
            'optimizer_coordinate_scales': scales.tolist(),
            'comparison_tolerances': {'scaled_distance': SCALED_DISTANCE_TOLERANCE,
                                      'absolute_chi2_difference': CHI2_TOLERANCE},
            'distance_definition': 'Euclidean norm of coordinate differences divided by optimizer scales',
            'candidates': candidates, 'ranking': [candidate['id'] for candidate in ranked],
            'pairwise_comparisons': pairs,
            'scope': 'Finite saved endpoints; common-order comparison does not prove a global minimum.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    for role, default in DEFAULT_INPUTS.items():
        parser.add_argument('--' + role.replace('_', '-'), default=default)
    args = parser.parse_args()
    paths = {role: getattr(args, role) for role in DEFAULT_INPUTS}
    if Path(args.output).resolve() in {Path(path).resolve() for path in paths.values()}:
        parser.error('output must not overwrite an input ledger')
    ledgers, hashes = read_ledgers(paths)
    model = joint_fit.GlobalLikelihood(4, 128, 128, bin_block_size=4)
    result = compare(ledgers, hashes, model)
    result['command'] = shlex.join([sys.executable, *sys.argv])
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['runtime_source_sha256'] = multi.source_hashes()
    multi.write_ledger(args.output, result)
    print(args.output)


if __name__ == '__main__':
    main()
