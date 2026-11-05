"""Compare saved joint-fit endpoints and choose a continuation point.

Run: .venv/bin/python scripts/select_joint_vector_continuation.py --ledger data/vector_sector/joint_vector_multistart.json --output /tmp/joint-vector-continuation.json
Coordinates use the units of the original fit; chi-squared and scaled distances
are dimensionless. This comparison makes no optimizer convergence claim.
"""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np


EXPECTED_NAMES = frozenset(('stored', 'preliminary', 'stored_negative_c4', 'stored_negative_odd'))
CONTINUATION_NAME = 'stored_negative_c4'


def compare_endpoints(ledger, source_path, source_sha256):
    """Rank the four physically valid endpoints and compute pairwise distances."""
    scales = np.asarray(ledger['optimizer_coordinate_scales'], dtype=float)
    names = ledger['parameter_names']
    if scales.shape != (32,) or len(names) != 32 or not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError('expected 32 named coordinates with finite positive scales')
    attempts = ledger['attempts']
    attempt_names = [attempt['name'] for attempt in attempts]
    if len(attempt_names) != 4 or set(attempt_names) != EXPECTED_NAMES:
        raise ValueError('expected exactly the four named starts')
    endpoints = {}
    for attempt in attempts:
        name = attempt['name']
        if attempt.get('endpoint_diagnostics', {}).get('valid') is not True:
            raise ValueError(f'{name}: endpoint is not physically valid')
        if type(attempt.get('success')) is not bool:
            raise ValueError(f'{name}: optimizer success is missing or invalid')
        values = np.asarray(attempt['endpoint'], dtype=float)
        chi2 = float(attempt['chi2'])
        if values.shape != (32,) or not np.all(np.isfinite(values)) or not np.isfinite(chi2):
            raise ValueError(f'{name}: endpoint or chi-squared is not finite')
        endpoints[name] = (attempt, values, chi2)
    ranked = sorted(endpoints, key=lambda name: (endpoints[name][2], name))
    if ranked[0] != CONTINUATION_NAME:
        raise ValueError(f'{CONTINUATION_NAME} is not the lowest-objective endpoint')
    selected = endpoints[ranked[0]][0]
    if selected['success']:
        raise ValueError('expected continuation point to retain unsuccessful optimizer status')
    comparisons = []
    for first, second in combinations(sorted(endpoints), 2):
        difference = (endpoints[first][1] - endpoints[second][1]) / scales
        distance = float(np.linalg.norm(difference))
        if not np.isfinite(distance):
            raise ValueError(f'{first}, {second}: nonfinite scaled distance')
        comparisons.append({'first': first, 'second': second, 'scaled_distance': distance})
    selected_endpoint = list(selected['endpoint'])
    if selected_endpoint != selected['endpoint']:
        raise ValueError('selected vector differs from source endpoint')
    return {'source_ledger': {'path': str(source_path), 'sha256': source_sha256},
            'parameter_names': names, 'optimizer_coordinate_scales': scales.tolist(),
            'distance_definition': 'Euclidean norm of coordinate differences divided by optimizer_coordinate_scales',
            'objective_ranking': [{'rank': index + 1, 'name': name,
                                   'chi2': endpoints[name][2],
                                   'success': endpoints[name][0]['success'],
                                   'endpoint_valid': True}
                                  for index, name in enumerate(ranked)],
            'pairwise_distances': comparisons,
            'continuation_seed': {'name': selected['name'], 'endpoint': selected_endpoint,
                                  'chi2': float(selected['chi2']),
                                  'success': selected['success'],
                                  'endpoint_valid': True,
                                  'accepted': False,
                                  'selection': 'lowest objective among the four physically valid endpoints'}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ledger', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    source = Path(args.ledger)
    destination = Path(args.output)
    if source.resolve() == destination.resolve():
        parser.error('output must not overwrite the source ledger')
    source_bytes = source.read_bytes()
    result = compare_endpoints(json.loads(source_bytes), str(source), hashlib.sha256(source_bytes).hexdigest())
    result['calculation_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    # Verify the emitted JSON reproduces the source vector exactly.
    if json.loads(destination.read_text())['continuation_seed']['endpoint'] != next(
            attempt['endpoint'] for attempt in json.loads(source_bytes)['attempts']
            if attempt['name'] == CONTINUATION_NAME):
        raise ValueError('written continuation vector differs from source endpoint')
    print(destination)


if __name__ == '__main__':
    main()
